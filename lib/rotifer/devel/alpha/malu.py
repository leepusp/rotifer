# Import rotifer modules
import rotifer

from rotifer.core.functions import loadConfig
config = loadConfig(__name__, defaults = {
    'read_tables': {
        'suffixes':[".scan.arch","_cluster.tsv"],
        'order':['c100i100','c80i70','c80i0','pfam','aravind'],
        'colnames':{'c100i100':['c100i100','pid'], 'c80i70':['c80i70','c100i100'], 'c80i0':['c80i0','c80i70']},
        'rename':{'pfam':{'ID':'c100i100','architecture':'pfam'}, 'aravind':{'ID':'c100i100','architecture':'aravind'}},
        'filter_columns':['pid','c100i100','c80i70','c80i0','pfam','aravind'],
    }})
logger = rotifer.logging.getLogger(__name__)

def findFilesByPrefixSuffix(
        prefix,
        suffixes=[".scan.arch","_cluster.tsv"],
        path=".",
        order=['c100i100','c80i70','c80i0','pfam','aravind']):
    import os
    from glob import glob
    import pandas as pd

    # Finding...
    tables = []
    for fname in glob(f'''{path}/{prefix}.*'''):
        name = os.path.basename(fname)
        name = name[len(prefix):]
        if name[0] == ".":
            name = name[1:]
        for pattern in suffixes:
            if name[-len(pattern):] == pattern:
                name = name[:len(name)-len(pattern)]
                if name[-1] == ".":
                    name = name[:-1]
                tables.append([ name, prefix, pattern, fname ])
    tables = pd.DataFrame(tables, columns=['target','prefix','suffix','path'])

    # Sorting
    if order:
        myorder = { order[i]: i for i in range(0,len(order)) }
    else:
        tables['idx'] = tables.index.tolist()
    tables['_idx'] = tables['target'].map(myorder)
    tables.sort_values(['_idx','target'], inplace=True)
    tables.reset_index(inplace=True, drop=True)
    tables.drop('_idx', inplace=True, axis=1)

    return tables

def read_tables(prefix,
                suffixes=config['read_tables']['suffixes'],
                path=".",
                order=config['read_tables']['order'],
                source=False,
                merge=True,
                colnames=config['read_tables']['colnames'],
                rename=config['read_tables']['rename'],
                filter_columns=config['read_tables']['filter_columns'],
                sep="\t",
                ):
    """
    This function loads, process and either join or concatenate a set of tables.

    The list of tables is derived iby matching prefixes and suffixes to a target directory.
    """
    import pandas as pd

    # Loading...
    df = pd.DataFrame()
    tables = findFilesByPrefixSuffix(prefix, suffixes, path, order)
    for idx, row in tables.iterrows():
        name = row['target']

        # Set column names
        columns = None
        if name in colnames:
            columns = colnames[name]

        # Read
        tmp = pd.read_csv(row['path'], sep=sep, names=columns)

        # Adjust
        if name in rename:
            tmp.rename(rename[name], inplace=True, axis=1)
        if filter_columns:
            tmp = tmp.filter(filter_columns)
        if source:
            tmp['source'] = name

        # Operate
        if df.empty:
            df = tmp.copy()
        elif merge:
            df = df.merge(tmp, how='left')
        else:
            df = pd.concat([ df, tmp ],ignore_index=True)

    # return dataframe
    return df

def extract_by_hmm(seqobj1, seqobj2, suffixes=('','_YyYyYy')):
    '''
    Use a HMM derived from seqobj2 to find similar regions
    in the sequences in seqobj1.

    If more than one match is found in a protein, two 
    different sequences will be returned with the same
    identifier and can only be distinguished using the
    start and end columns.

    Returns
    -------
    Unaligned sequence set object
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import sys
    import pandas as pd
    import rotifer
    #from rotifer.core import loadpath
    
    hmmer2table = f'{rotifer.config["base"]}/bin/hmmer2table'
    domain2architecture = f'{rotifer.config["base"]}/bin/domain2architecture'
    architecture2table = f'{rotifer.config["base"]}/bin/architecture2table'

    sharedColumns = set()
    nseqobj = seqobj1.copy()
    with tempfile.TemporaryDirectory() as tmpdirname:
        nseqobj.df.sequence = nseqobj.df.sequence.str.replace('-','')
        nseqobj.to_file(f'{tmpdirname}/seqfile') 

        if isinstance(seqobj2, rotifer.devel.beta.sequence.sequence):
            seqobj2.to_file(f'{tmpdirname}/model')
            Popen(f'hmmbuild {tmpdirname}/model.hmm {tmpdirname}/model',
                  stdout=PIPE,
                  shell=True).communicate()
            model = f'{tmpdirname}/model.hmm'

        Popen(f'hmmsearch {model} {tmpdirname}/seqfile > {tmpdirname}/hhsearch_r',
              stdout=PIPE,
              shell=True).communicate()
        exe  = f'{hmmer2table} {tmpdirname}/hhsearch_r'
        exe += f' | {domain2architecture} | {architecture2table} > {tmpdirname}/hhsearch_table'
        Popen(exe, stdout=PIPE, shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/hhsearch_table', sep="\t")
        t.rename({'ID':'id'}, axis=1, inplace=True)
        sharedColumns = set(t.columns).intersection(set(seqobj1.df.columns))

    if "start" in sharedColumns:
        start = f'start{suffixes[1]}'
    else:
        start = "start"
    if "end" in sharedColumns:
        end = f'end{suffixes[1]}'
    else:
        end = "end"

    nseqobj.df = nseqobj.df.merge(t, on='id', how='left', suffixes=suffixes)
    nseqobj.df[end] = nseqobj.df[end].fillna(nseqobj.df.length)
    nseqobj.df[start] = nseqobj.df[start].fillna(1)
    nseqobj.df.sequence = nseqobj.df.apply(lambda x: x.sequence[(int(x[start])-1):int(x[end])], axis=1)
    nseqobj = nseqobj.filter('domain == "model"')
    return (nseqobj)

def positions_to_coordinates(seqobj, df, seqid='ID', annotation='fixed', start='start', end='end', maxoverlap=10, how='merge'):
    coord = df.filter([seqid,annotation,start,end]).values
    coord = [ (*seqobj.position_to_column([s,e],i),f,i) for i, f, s, e in coord  ]
    coord = pd.DataFrame(coord, columns=["start","end","annotation","seqid"])
    coord.sort_values(['start','end'], inplace=True)
    coord['same'] = ((coord.annotation != coord.annotation.shift(1)) | (coord.start > maxoverlap+coord.end.shift(1)))
    coord['same'] = coord.same.cumsum()
    if how == "merge":
        coord = coord.groupby('same').agg({'start':'min', 'end':'max', 'annotation':'first'}).values
    elif how == "longest":
        coord['length'] = coord.end - coord.start + 1
        coord.sort_values(['same','length'], ascending=[True,False], inplace=True)
        coord = coord.drop_duplicates('same')
    coord = [ tuple(x) for x in coord.filter(['start','end','annotation']).values ] #substitui ':' por ',' nessa linha
    return coord

def annotate_columns_from_sequence_coordinates(seqobj, df, **kwargs):
    from rotifer.devel.alpha import gian_func as gf
    return gf.annotation(seqobj, positions_to_coordinates(seqobj, df, **kwargs))

def merge_ranges(df, identifier=None, merge_distance=10, method=None):
    """
    Merge overlapping or nearby (within merge_distance) regions for each sequence ID.
    Takes a df with an identifier colum, a merge distance that will be used to fuse sequence lengths
    Sometimes the hits of the domain are distant enough that they won't merge, currently the only selection method is sliced length

    df : A df with start and end columns and a sequence identifier
    identifier: the column which contains the sequence identifiers
    merge_distance : how far the start and end position of the domain hits can be from each other
    method: 'length', if the sequences are not mergeable, chooses the largest fragment
    it may be nice to select by evalue! Not implemented yet

    """
    import pandas as pd

    merged = []
    df = df.copy()

    if identifier:
        df.rename({identifier: 'id'}, axis=1, inplace=True)

    for _, group in df.groupby("id"):
        # Sort ranges for each id
        ranges = group[["start", "end"]].sort_values("start").values.tolist()
        merged_ranges = []
        for rng in ranges:
            if not merged_ranges:
                merged_ranges.append(rng)
            else:
                last = merged_ranges[-1]
                if rng[0] <= last[1] + merge_distance:
                    # Merge overlapping/nearby ranges
                    last[1] = max(last[1], rng[1])
                else:
                    merged_ranges.append(rng)

        for start, end in merged_ranges:
            merged.append({"id": group["id"].iloc[0], "start": start, "end": end})

    result = pd.DataFrame(merged)
    result["length"] = result["end"] - result["start"]

    if method == 'length':
        result = result.sort_values(['id', 'length'], ascending=[True, False])
        result = result.groupby('id', as_index=False).first()

    return result
                    


def extract_by_table(seqobj, df, suffixes=('','_YyYyYy'),method=None, merge_distance=10, downstream=0,upstream=0,):
    """

    Extracts the domains by the provided df, **needs to be for only one domain**, currently seems to add +1 residue
    good enough for now, but needs checking later

    seqobj: your sequence object
    df : the df with the hits, must contain at least the columns ID, start and end.
    method: see merge_ranges
    merge_distance: see merge_ranges
    downstream : number of bases downstream
    upstream: number of bases upstream
    """

    import pandas as pd

    t = df.copy()
    t.rename({'ID': 'id'}, axis=1, inplace=True)
    nseqobj=seqobj.copy()
    if {'start', 'end'}.issubset(t.columns):
        t = merge_ranges(t, merge_distance=merge_distance,method=method)

    sharedColumns = set(t.columns).intersection(set(nseqobj.df.columns))

    if "start" in sharedColumns:
        start = f'start{suffixes[1]}'
    else:
        start = "start"
    if "end" in sharedColumns:
        end = f'end{suffixes[1]}'
    else:
        end = "end"

    nseqobj.df = nseqobj.df.merge(t, on='id', how='left', suffixes=suffixes)
    nseqobj.df[end] = nseqobj.df[end].fillna(nseqobj.df.length)
    nseqobj.df[start] = nseqobj.df[start].fillna(1)                                                                                                                                                              
    nseqobj.df.sequence = nseqobj.df.apply(lambda x: x.sequence[max(0, int(x[start]) - upstream - 1):min(int(x[end]) + downstream, int(x['length']))], axis=1)
    nseqobj.df['length_slice'] = nseqobj.df.sequence.str.len()
    return nseqobj



def make_hist(df, column, bins=10):
      """
      This method displays an histogram for a given int column

      Parameters
      ----------
      column: string
          the name of the column for which the histogram is desired
      bins : int
          Number of histogram bins to discretize the distribution

      Returns
      -------
      None

      Examples
      --------
      >>> make_hist(df,'column', 20)
      """
      from ascii_graph import Pyasciigraph
      import collections
      import pandas as pd
      # Raise error if the column given is not numeric
      if not pd.api.types.is_numeric_dtype(df[column]):
          raise TypeError(f"Column '{column}' must be numeric.")

      collections.Iterable = collections.abc.Iterable
      print(f'Number of lines analyzed: {len(df[column])}')
      a = df[column].value_counts().to_frame().reset_index()
      a.columns = [column, 'count']
      a = a.sort_values(column)
      a['raw_bin'] = pd.cut(a[f'{column}'],bins,precision=0)
      a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
      test = a.groupby('bin').agg({'count':'sum'}).reset_index().apply(tuple, axis=1)
      graph = Pyasciigraph()
      for line in  graph.graph('count \t ' f'{column}', test):
          print(line)


### Igem 2026


def filter_fimo(fimoraw, gentab, repdist=2, gensize=5):
    """
    Apply the repeat-distance filter and mark repeats that overlap genes
    on the same nucleotide and strand.

    Parameters
    ----------
    fimoraw : pandas.DataFrame
        Must contain sequence_name, strand, start, stop, and distance.

    gentab : pandas.DataFrame
        Must contain nucleotide, strand, start, end, plen, and type.

    repdist : int or float, default=2
        Minimum repeat distance retained.

    gensize : int or float, default=5
        Minimum protein length used in the overlap check.

    Returns
    -------
    pandas.DataFrame
        Distance-filtered repeats with an ``intragenic`` column.
    """
    import numpy as np
    import pandas as pd

    # Doing .copy() to avoid pandas' warnings

    distfilt = fimoraw.loc[fimoraw["distance"] >= repdist].copy()
    
    # I'm  only interested in genes, so I grab only the CDS and while I am at it, I will ignore proteins that are smaller than 5 AA
    
    genfil = gentab.loc[(gentab["type"] == "CDS") & (gentab["plen"] >= gensize)].copy()

    # I need to compare repeats with genes on the same strand, but FIMO has "+" and "-" representation, so I'll change that on the fly
    
    distfilt["strand"] = distfilt["strand"].map({"+": 1,"-": -1})
    
    # Crash if there are unexpected values in strand
    
    if distfilt["strand"].isna().any():
        raise ValueError("FIMO strand contains values other than '+' and '-'")

    # Now, the operation I intend to do wants the df to be ordered

    distfilt = distfilt.sort_values(["sequence_name", "start", "stop","strand"])
    genfil = genfil.sort_values(["nucleotide", "start", "end","strand"])

    # Create gene groups where only one strand is present
    
    gene_groups = {key: group for key, group in genfil.groupby(["nucleotide", "strand"],sort=False,dropna=False)}

    # Initialize  series to store results defaulting to False
    
    result = pd.Series(False,index=distfilt.index,dtype=bool)

    # Iterate for every matching key 
    
    for key, repeats in distfilt.groupby(["sequence_name", "strand"],sort=False,dropna=False):
        
        # Grab the matching key
        genes = gene_groups.get(key)
        
        # Handling if it continues
        if genes is None or genes.empty:
            continue

        # NP arrays from each group
        
        gene_starts = genes["start"].to_numpy()
        gene_ends = genes["end"].to_numpy()

        repeat_starts = repeats["start"].to_numpy()
        repeat_stops = repeats["stop"].to_numpy()

        # Largest gene end observed at or before each gene position.
        max_end_so_far = np.maximum.accumulate(gene_ends)

        # Last gene whose start is <= repeat stop.
        pos = np.searchsorted(gene_starts, repeat_stops,side="right") - 1

        valid = pos >= 0
        overlap = np.zeros(len(repeats), dtype=bool)

        overlap[valid] = (
            max_end_so_far[pos[valid]]
            >= repeat_starts[valid]
        )

        # Assign the value
        result.loc[repeats.index] = overlap

    distfilt["intragenic"] = result

    return distfilt