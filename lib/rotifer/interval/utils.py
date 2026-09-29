#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from collections import OrderedDict

import rotifer
from rotifer.core.functions import loadConfig
from rotifer.devel.beta import sequence as rdbs
logger = rotifer.logging.getLogger(__name__)
config = loadConfig(__name__, defaults = {
    'default': {
        'reference': ['query'],
        'collapse': 'hit',
        'start': 'qstart',
        'end': 'qend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'sum_probs':False, 'probability':False, 'score':False, 'region_length':False}),
    },
    'archtsv': {
        'reference': ['ID'],
        'collapse': 'domain',
        'start': 'start',
        'end': 'end',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'evalue':True, 'region_length':False}),
    },
    'hhsuite': {
        'reference': ['query'],
        'collapse': 'hit',
        'start': 'qstart',
        'end': 'qend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'sum_probs':False, 'probability':False, 'score':False, 'region_length':False}),
    },
    'hmmer': {
        'reference': ['sequence'],
        'collapse': 'model',
        'start': 'estart',
        'end': 'eend',
        'maximum_overlap': 0.4,
        'criteria': OrderedDict({'evalue':True, 'region_length':False}),
    },
})

def architecture(df,
        reference       = config['default']['reference'],
        collapse        = config['default']['collapse'],
        start           = config['default']['start'],
        end             = config['default']['end'],
        maximum_overlap = config['default']['maximum_overlap'],
        criteria        = config['default']['criteria'],
    ):
    '''
    This function builds a compact representation of all non-overlapping
    regions for a given reference column.

    It is often used to produce a compact protein domain architecture table

    The input Pandas DataFrame must contain columns corresponding
    to unidimensional interval coordinates (start and end columns)
    and at least one reference column (e.g. pid or sequence).

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    reference : (list of) strings, default ['query']
        Name of the column storing sequence identifiers
    collapse : (list of) strings, default 'hits'
        Name of the column storing region identifiers
    start : string, default estart
        Name of the column with the first coordinate of each interval
    end : string, default eend
        Name of the column with the last coordinate of each interval
    maximum_overlap : integer, default 0.1 (10% of the shortest)
        Maximum number of overlapping residues to ignore.
        If maximum_overlap < 1, the parameter is treated as a percentage 
    criteria : ordered dictionary
        Criteria for best interval selection.

        Keys in this dictionary must match column names in the input
        Pandas DataFrame. Values should be booleans (True or False).

        'region_length' is an exception: this column is automatically
        calculated from the start and end coordinates of each region.

        A True value for a given column will give precedence for rows
        with the smallest value in that column. A False value
        will select rows with the largest values.

    Examples
    --------
    HH-suite results based on the largest probability and/or region length:

    >>> from rotifer.interval.utils import architecture
    >>> df = pd.read_csv("results.hhsuite2table.tsv", sep="\t")
    >>> dfc = architecture(df, start='qstart', end='qend', criteria=c)
    '''
    arch = filter_nonoverlapping_regions(df, reference=reference, start=start, end=end, maximum_overlap=maximum_overlap, criteria=criteria)
    arch.sort_values([*reference, start, end], ascending=True, inplace=True)
    arch = arch.groupby(reference).agg(architecture=(collapse,lambda x: "+".join([ str(y) for y in x ]))).reset_index()
    return arch

def filter_nonoverlapping_regions(df,
        reference       = config['default']['reference'],
        start           = config['default']['start'],
        end             = config['default']['end'],
        maximum_overlap = config['default']['maximum_overlap'],
        criteria        = config['default']['criteria'],
        *args, **kwargs
    ):
    '''
    This function selects the best set of non-overlapping intervals.

    The input Pandas DataFrame must contain columns corresponding
    to unidimensional interval coordinates (start and end columns)
    and at least one reference column (e.g. pid or sequence).

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    reference : (list of) strings, default ['sequence']
        Name of the column storing sequence identifiers
    start : string, default estart
        Name of the column with the first coordinate of each interval
    end : string, default eend
        Name of the column with the last coordinate of each interval
    maximum_overlap : integer, default 0.1 (10% of the shortest)
        Maximum number of overlapping residues to ignore.
        If maximum_overlap < 1, the parameter is treated as a percentage 
    criteria : ordered dictionary
        Criteria for best interval selection.

        Keys in this dictionary must match column names in the input
        Pandas DataFrame. Values should be booleans (True or False).

        'region_length' is an exception: this column is automatically
        calculated from the start and end coordinates of each region.

        A True value for a given column will give precedence for rows
        with the smallest value in that column. A False value
        will select rows with the largest values.

    Examples
    --------
    HMMER results, parsed by hmmer2table, using the environment limits:

    >>> df = pd.read_csv("results.hmmer2table.tsv", sep="\t")
    >>> dfc = filter_overlapping_regions(df)

    HH-suite results based on the largest probability and/or region length:

    >>> df = pd.read_csv("results.hhsuite2table.tsv", sep="\t")
    >>> c = {'probability':False, 'region_length':False}
    >>> dfc = filter_overlapping_regions(df, start='qstart', end='qend', criteria=c)
    '''
    from copy import deepcopy

    # Process arguments and select columns for analysis
    if not isinstance(criteria,OrderedDict):
        criteria = OrderedDict(criteria)
    if not isinstance(reference,list):
        reference = [ reference ]

    # Copy and sort the dataframe
    # Sort row priority: best rows must appear first
    crit = list(criteria.keys())
    clean = deepcopy(df)
    remove = False
    if 'region_length' in crit and 'region_length' not in clean.columns:
        clean['region_length'] = clean[end] - clean[start] + 1
        remove = True
    clean.sort_values(reference + crit, ascending=len(reference)*[True] + [ criteria[x] for x in crit ], inplace=True)
    clean.reset_index(drop=True, inplace=True)

    # Remove non-optimal layers
    keep = pd.Series([ pd.NA ] * len(clean), index=clean.index.copy(), dtype="boolean")
    while keep.isna().any():
        logger.info(f'Running next layer with {keep.notna().sum()} processed rows out of {len(clean)} rows')
        best = clean.loc[keep.isna()].drop_duplicates(reference)
        keep.loc[best.index] = True
        best = clean.loc[keep == True]
        worst = clean.loc[keep.isna()].reset_index(drop=False).rename({'index':'_oid'}, axis=1)
        worst = best.merge(worst, on=reference, how='inner')
        overlap  = np.minimum(worst[f'{end}_x'],   worst[f'{end}_y'])
        overlap -= np.maximum(worst[f'{start}_x'], worst[f'{start}_y'])
        overlap += 1
        maxoverlap = maximum_overlap
        if maxoverlap < 1:
            maxoverlap = np.floor(maxoverlap * np.minimum(worst.region_length_x, worst.region_length_y))
        overlap = (overlap > maxoverlap)
        overlap = overlap | ((worst[f'{start}_x'] <= worst[f'{start}_y']) & (worst[f'{end}_x'] >= worst[f'{end}_y']))
        worst = worst[overlap]._oid.unique().tolist()
        keep.loc[worst] = False

    # Remove region_length and return
    if remove:
        clean.drop('region_length', axis=1, inplace=True)
    clean = clean.loc[keep].copy()
    clean.sort_values(reference + [start,end], ascending=(len(reference)+2)*[True], inplace=True)
    return clean

def complement(df, include=True, exclude=False, add=False, reference='ID', start='start', end='end'):
    
    '''
    Given a DataFrame describing regions in one-dimensional
    reference systems, return a similar DataFrame describing
    the set of regions required to fully cover all the length
    of the references.

    The input table is expected to include one or more columns
    describing the reference object and a pair of integer columns
    corresponding to the start and end positions of each region.

    If a column with the length of each interval is supplied, the
    last interval and the length of each fragment will be included
    in the returned dataframe.

    Parameters
    ----------
    df : Pandas dataframe
       The dataframe must have at least a reference column and coordinates
    include: list
    exclude: list
    add: dict
    reference : (list of) strings, default ['sequence']
        Name of the column(s) storing sequence identifiers
    start : string, default 'start'
        Name of the column with the first coordinate of each interval
    end : string, default 'end'
        Name of the column with the last coordinate of each interval
    '''
    if isinstance(reference, str):
        reference = [ reference ]

    regions = df.sort_values([*reference, start, end])
    if 'length' not in regions.columns:
        regions['length'] = regions.ID.map(regions.set_index(reference).end.to_dict())
        regions.loc[regions.length.isna(), 'length'] = regions.loc[regions.length.isna(), 'end']

    first = pd.DataFrame({
        'ID':     regions.ID.tolist(),
        'start':  np.where(regions.ID != regions.ID.shift(1),1,regions.end.shift(1) + 1).tolist(),
        'end':    (regions.start - 1).tolist(),
    }).astype({'start':int}).eval('length = end - start + 1').query('end > 0')

    uncovered = (regions.groupby(['ID','length']).end.max()+1)
    uncovered = uncovered.reset_index()
    uncovered = uncovered.query('end < length')
    uncovered = uncovered.rename({'end':'start'}, axis=1)
    uncovered = uncovered.eval('end = length').eval('length = end - start + 1')
    uncovered = pd.concat([ first, uncovered ]).sort_values(['ID','start','end'])

    if isinstance(include, list):
        uncovered = uncovered.filter(include).reindex(df.columns, axis=1)
    elif include:
        uncovered = uncovered.reindex(df.columns, axis=1)

    if isinstance(exclude, list):
        uncovered = uncovered.drop(exclude, axis=1)

    if isinstance(add, dict):
        for keys in values.keys():
            uncovered[keys] = values[keys]

    return uncovered

