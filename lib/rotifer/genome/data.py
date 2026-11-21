#!/usr/bin/env python3

# Data here is a child of Pandas DataFrame

from collections import OrderedDict
import itertools
import numpy as np
import pandas as pd
import rotifer
import rotifer.pipeline
from rotifer.core.functions import not_kwargs
from rotifer.taxonomy.utils import lineage as rtlineage
from rotifer.genome.utils import seqrecords_to_dataframe
logger = rotifer.logging.getLogger(__name__)

_column_dict = OrderedDict({
  'nucleotide'     :'nucleotide',
  'start'          :'location',
  'end'            :'location',
  'strand'         :'location',
  'nlen'           :'nucleotide',
  'block_id'       :'block',
  'query'          :'block',
  'pid'            :'feature',
  'type'           :'feature',
  'plen'           :'feature',
  'locus'          :'feature',
  'seq_type'       :'nucleotide',
  'assembly'       :'nucleotide',
  'gene'           :'feature',
  'origin'         :'feature',
  'topology'       :'nucleotide',
  'product'        :'feature',
  'taxid'          :'taxonomy',
  'organism'       :'taxonomy',
  'lineage'        :'taxonomy',
  'classification' :'taxonomy',
  'feature_order'  :'feature',
  'internal_id'    :'feature',
  'is_fragment'    :'nucleotide',
  'replaced'       :'feature',
})

def _check_choice(name, value, allowed):
    """Raise ValueError unless ``value`` is one of ``allowed``.

    The message names the argument, the value given and the values
    accepted.
    """
    ok = (value is None and None in allowed) or (isinstance(value, str) and value in allowed)
    if not ok:
        raise ValueError(
            f"{name}={value!r} is not supported; use one of "
            + ", ".join(repr(a) for a in allowed))


class NeighborhoodDF(pd.DataFrame, rotifer.pipeline.Annotatable):
    _metadata = ['filterby','NDFProperties']
    NDFProperties = {
        "columns": list(_column_dict.keys()),
        "required": [ 'assembly','nucleotide','start','end','strand','type' ],
        "order": [ 'assembly','nucleotide','start','end','strand' ],
    }

    def __init__(self, *args, **kwargs):
        self.filterby  = kwargs.pop('filterby', None)
        update_lineage = kwargs.pop('update_lineage',False)
        preferred_taxa = kwargs.pop('preferred_taxa',None)
        if len(args) == 0 and 'columns' not in kwargs:
            kwargs['columns'] = __class__.NDFProperties["columns"]

        # Call super().__init__()
        super(NeighborhoodDF, self).__init__(*args, **kwargs)

        # Update or add lineage column
        cols = self.columns.to_list()
        if isinstance(update_lineage,str) & (update_lineage in cols):
            if "lineage" in cols:
                self.lineage = rtlineage(self[update_lineage], preferred_taxa=preferred_taxa)
            else:
                taxonomy = [ i for i in range(0,len(cols)) if cols[i] == update_lineage ][0]
                self.insert(taxonomy, "lineage", rtlineage(self[update_lineage], preferred_taxa=preferred_taxa))

    # Pandas contructor
    @property
    def _constructor(self):
        return NeighborhoodDF



###### New functions added 12/08/2019

    def group_strand(self):
        """
        A function to add the column strand_groups in the acc2operon dataframe.
        """
        dataframe = self.copy()
        dataframe['strand_group'] = np.cumsum(
            np.where(dataframe.block_id == dataframe.block_id.shift(1),
                     np.where(
                         dataframe.strand.shift(1) == dataframe.strand, 0, 1
                     ), 1)
        )
        return dataframe

    def same_strand(self,
                    query='query'):
        """
        Using the strand group and a query it will filter the self with
        strands that present one  query at least.
        """
        dataframe = self.group_strand().copy()
        return dataframe[dataframe.strand_group.isin(
            dataframe.groupby('strand_group')[query].sum().where(
                lambda x: x > 0
            ).dropna().index)]

    def value_counts_full(self, column):
        """
        Calculate absolute and relative frequencies for a given column.
        If dropna is True, it will ignore the empty values to count
        the relative frequency.
        """
        x = self[self['type']=="CDS"]
        lpid = x['pid'].nunique()
        lbid = x['block_id'].nunique()
        g_df = x.groupby(column).agg({'block_id':'nunique','pid':'count'}).sort_values('block_id', ascending=False).reset_index()
        g_df['relative_block_id'] = g_df.block_id/lbid *100
        g_df['relative_pid'] = g_df.pid/lpid * 100
        g_df.rename({'block_id':'block_id_count', 'pid':'pid_count'}, axis=1)
        logger.info(f'Total of block_id in selection:{lbid}')
        logger.info(f'Total of pid in selection :{lpid}')

        return g_df

    def select_neighbors(self,
                         column='query',
                         value=1,
                         new_query=True,
                         **kwargs):
        """
        Slice the self in x neighbors of a given value query.
        self : The Dataframe to be sliced.
        Column : The column that will be filtered.
        value: The value to be filtered.
        after and before: The numbers of neighbors filtered.
        strand:Filter for neighbors that are at the same strand of the query
        stats: If set it will give the stats of the column
        Examples:
        get_x_neighbors(DF_ATPases, column='Pfam', value='AAA_9',before=3,after=3).
        In this example the Dataframe DF_ATPases will be filtered by only blocks
        that present the AAA_9 ATPase and will return 3 before and 3 after neighbors
        get_x_neighbors(DF_ATPases, column='cluster', value='sc_29',
        before=3,after=3, stats='pfam).
        The output of this example will be the stats of pfam column
        in the neighborhood of cluster sc_29.
        """

        after = not_kwargs(kwargs, 'after', 1)
        before = not_kwargs(kwargs, 'before', 1)
        strand = not_kwargs(kwargs, 'strand', False)
        stats = not_kwargs(kwargs, 'stats', False)

        query_df = self[self[column] == value].copy()
        full_query_df = \
                self[self.block_id.isin(query_df.block_id)].copy()
        min_max = full_query_df.reset_index().groupby(
            'block_id'
        ).agg({'index': ['min', 'max']})
        min_max.columns = min_max.columns.droplevel()
        query_min_max = \
                query_df.reset_index()[['block_id', 'index']].set_index('block_id')
        query_min_max = query_min_max.join(min_max)
        query_min_max['before_limit'] = 1 + (
            np.where(query_min_max['index'] + before >= query_min_max['max'],
                     query_min_max['max'],
                     query_min_max['index'] + before)
        )
        query_min_max['after_limit'] = np.where(
            query_min_max['index'] -after <= query_min_max['min'],
            query_min_max['min'],
            query_min_max['index'] - after
        )
        query_min_max['Diference'] = \
                query_min_max.before_limit - query_min_max.after_limit
        index_list = list(
            itertools.chain.from_iterable(
                [list(range(x[0], x[1])) for x in list(
                    zip(query_min_max.after_limit, query_min_max.before_limit))]))
        selected_df = self.iloc[index_list].copy()
        selected_df.loc[:, 'new_query'] = 0
        selected_df.loc[selected_df[column] == value, 'new_query'] = 1
        if strand:
            if stats:
                return selected_df.same_strand(query='new_query').value_counts_full(stats)
            if new_query:
                return selected_df.same_strand(query='new_query')
            return selected_df.same_strand(query='new_query').drop(
                'new_query', axis=1)

        if stats:
            return selected_df.value_counts_full(stats)
        if new_query:
            return selected_df
        return selected_df.drop('new_query', axis=1)

    def compact_neighborhood(self,
                             annotation='pfam',
                             query='query',
                             strand=False,
                             fill_unk=False):
        from rotifer.view import functions
        from rotifer.devel.alpha import gian_func as gf
        jupyter = functions.is_running_in_jupyter()
        if jupyter:
            to_display = gf.operon_fig2(self,domain_column = annotation, to_string=True) 
            functions.display_html_popup_from_file(to_display)
            return
        """
        A tabular and compact representation of  selected neighborhoods block.
        Given a acc2operons Dataframe it will print the block ID as first column
        and a compact form of all locus within the block_id as compact form using
        a given
        column to annotate it.
        If stand is True, it will display only the loci that are at the same
        strand of the query, and also the ID displayed is a new calculated Strand ID.
        """

        columns = ['compact', 'block_id', 'lineage', 'organism', 'seq_type',
                   'classification', 'pid_compact']
        carried = columns[1:-1]

        df = self[self['block_id'].notna()]
        if df.empty:
            return pd.DataFrame(columns=columns).set_index('block_id')

        # Strand runs: a new run starts where a gene's strand differs from
        # the previous gene of the same block. Comparing within the block,
        # not with the row above, keeps a block's runs intact when its rows
        # are not contiguous. Runs are numbered in order of appearance.
        previous_strand = df['strand'].groupby(df['block_id']).shift(1)
        new_run = (df['strand'] != previous_strand) | previous_strand.isna()
        run_in_block = new_run.groupby(df['block_id']).cumsum()
        run = pd.Series(
            pd.factorize(pd.MultiIndex.from_arrays([df['block_id'], run_in_block]))[0],
            index=df.index)

        # One output row per group: a block, or with strand=True a run of
        # same-strand genes that contains a query. A row is a query when
        # query == 1. The key is named so that grouping does not clash with
        # the block_id column.
        if strand:
            keep = run.isin(run[df[query] == 1])
            df, run = df[keep], run[keep]
            if df.empty:
                return pd.DataFrame(columns=columns).set_index('block_id')
            key = run.rename('_group')
        else:
            key = df['block_id'].rename('_group')

        # Each group is read in the direction of its first query row (in
        # table order): forward if that row is on the plus strand, backward
        # otherwise. Backward groups are joined in reverse, so rows are never
        # sorted.
        is_query = df[query] == 1
        first_query = pd.DataFrame({
            'key': key[is_query],
            'strand': df.loc[is_query, 'strand'],
        }).drop_duplicates('key')
        no_query = ~key.isin(first_query['key'])
        if no_query.any():
            blocks = pd.unique(df.loc[no_query, 'block_id'])
            listed = ', '.join(str(b) for b in blocks[:10])
            more = f' and {len(blocks) - 10} more' if len(blocks) > 10 else ''
            raise NameError(f'Query_Missing: no query in block_id(s) {listed}{more}')
        forward_by_key = pd.Series((first_query['strand'] == 1).to_numpy(),
                                   index=first_query['key'].to_numpy())
        forward = key.map(forward_by_key).to_numpy(dtype=bool)

        # '|' goes after a gene whose next gene, in reading order, belongs to
        # another strand run. The last gene of a group gets none.
        next_run = run.groupby(key).shift(-1)
        previous_run = run.groupby(key).shift(1)
        bar = np.where(forward,
                       next_run.notna() & (run != next_run),
                       previous_run.notna() & (run != previous_run))

        # Text for each gene: '?' when unannotated, the feature type for
        # non-CDS rows, and optionally the fill_unk column for '?'. A value
        # still missing after that (no type, or an empty fill_unk value) is
        # shown as '?' too. The annotation pass runs first and the pid pass
        # sees its result, as fill_unk may name either column.
        work = df[list(dict.fromkeys(
            [annotation, 'pid', 'type'] + ([fill_unk] if fill_unk else [])
        ))].copy()

        def text(column):
            work.fillna({column: '?'}, inplace=True)
            work[column] = np.where(work['type'] != 'CDS', work['type'], work[column])
            if fill_unk:
                work[column] = np.where(work[column] == '?', work[fill_unk], work[column])
            work[column] = work[column].fillna('?').astype(str)
            return work[column]

        # Strands are flipped for backward groups, so the query always points
        # forward.
        plus = ((df['strand'] * np.where(forward, 1, -1)) == 1).to_numpy()
        star = is_query.to_numpy()

        def tokens(values):
            token = np.where(plus,
                             np.where(star, '*' + values + '->', values + '->'),
                             np.where(star, '<-' + values + '*', '<-' + values))
            return np.where(bar, token + '|', token)

        annotation_tokens = tokens(text(annotation))
        pid_tokens = tokens(text('pid'))
        token_frame = pd.DataFrame({'compact': annotation_tokens,
                                    'pid_compact': pid_tokens}, index=df.index)
        joined = pd.concat([
            token_frame[forward].groupby(key[forward]).agg(''.join),
            token_frame[~forward].groupby(key[~forward]).agg(lambda x: ''.join(x.iloc[::-1])),
        ]).sort_index()

        # Carried columns take the first non-missing value in reading order,
        # which for backward groups is the last one in the table.
        by_key = df[carried].groupby(key)
        first, last = by_key.first(), by_key.last()
        backward_keys = forward_by_key.index[~forward_by_key.to_numpy()]
        first.loc[first.index.isin(backward_keys)] = last.loc[last.index.isin(backward_keys)]

        result = pd.concat([joined, first], axis=1)[columns]
        return result.reset_index(drop=True).set_index('block_id')

    #Alias for compact neighborhood
    compact = compact_neighborhood

    # Find internal_id (iimin,iimax) and feature_order (fomin,fomax) boundaries
    def boundaries(self, nucleotide=None):
        """
        Find maximum and minimum values of block_id, feature_order and internal_id.
        Columns feature_order and internal_id are evaluated per nucleotide entry.

        Returns a Pandas dataframe.
        """
        df = self
        if isinstance(nucleotide,list) or isinstance(nucleotide,pd.Series):
            df = self[self.nucleotide.isin(nucleotide)]
        dflim = df.groupby(['assembly','nucleotide','block_id','type']).agg({'feature_order':['min','max']})
        dflim.reset_index(inplace=True)
        dflim.columns = ['assembly','nucleotide','block','type','fomin','fomax']
        iidlim = df.groupby(['assembly','nucleotide','block_id']).agg({'internal_id':['min','max']})
        iidlim.reset_index(inplace=True)
        iidlim.columns = ['assembly','nucleotide','block','iidmin','iidmax']
        bidlim = df.groupby(['assembly']).agg({'block_id':[lambda x: 1,'nunique']})
        bidlim.reset_index(inplace=True)
        bidlim.columns = ['assembly','bidmin','bidmax']
        dflim = dflim.merge(iidlim, on=['assembly','nucleotide','block'], how='left')
        dflim = dflim.merge(bidlim, on=['assembly'], how='left')
        return dflim

    def vicinity(self, targets=['query == 1'], before=3, after=3, min_block_distance=0, fttype='same'):
        """
        Locate genomic regions that contain features selected by some criteria.

        Each region will be described by assembly, nucleotide and a pair of minimum and
        maximum internal IDs.

        The user may choose a set of rows (targets) as anchors, whose neighbors will be
        evaluated by user-defined parameters, such as feature type, strand or distance.

        Parameters
        ----------
        targets: string, (list of) booleans or pd.Series
          Rules used to identify target rows.

          Rules are based on Pandas's eval method or its results.

          If all rules return True for a given row, the row is a
          target and its neighbors will be identified.
        before: int
          Keep at most this number of features *before* each target
        after: int
          Keep at most this number of features *after* each target
        min_block_distance: int
          Minimum distance between two blocks.

          Blocks that are separated by smaller number of features
          will be merged into a single block,
        fttype: str
          How to count neighbors by feature type.
          Supported values:
          ```same```
            Consider only features of the same type as the target.
            This is the only supported value; any other raises ValueError.

        Raises
        ------
        ValueError
            If ``fttype`` is not ``'same'`` or a string rule in ``targets``
            cannot be evaluated.
        TypeError
            If an item of ``targets`` is neither a string nor a boolean Series.
        """
        _check_choice('fttype', fttype, ['same'])
        import sys
        import numpy as np
        from copy import deepcopy
        blks = deepcopy(self)

        # Make sure data is sorted and has compatible internal ids
        blks.sort_values(['assembly','nucleotide','block_id','internal_id'], ascending=[True,True,True,True], inplace=True)

        # Build boolean pandas.Series to mark targets
        select = True
        if not isinstance(targets, list):
            targets = [targets]
        for code in targets:
            if isinstance(code,str):
                try:
                    select &= blks.eval(code)
                except Exception as exc:
                    raise ValueError(f'Rule {code!r} failed: {exc}') from exc
            elif isinstance(code,pd.Series):
                select &= code
            else:
                raise TypeError(f'targets must be strings or boolean Series, not {type(code).__name__}')

        # Find targets and their neighborhood
        if (not isinstance(select,pd.Series)) or (select.sum() == 0):
            logger.error(f'No anchors to define neighborhoods! Targets: {targets}')
            return pd.DataFrame()

        # Initialize dataframe for each region (blocks)
        cols = ['assembly','nucleotide','topology','type']
        dflim = blks.boundaries()

        # Process features with type
        if (fttype == 'same'):
            blksCols = [ x for x in [*cols,'block_id','feature_order','internal_id','is_fragment'] if x in blks.columns ]
            blks = blks[select].filter(blksCols).copy()
            if 'is_fragment' not in blks.columns:
                blks['is_fragment'] = False
            blks['foup']   = blks.feature_order - before
            blks['fodown'] = blks.feature_order + after

            # Identify blocks by merging neighborhoods
            blks.sort_values([*cols,'feature_order'], inplace=True)
            bid = (blks.nucleotide != blks.nucleotide.shift(1))
            bid = bid | (blks.is_fragment & (blks.block_id != blks.block_id.shift(1)))
            bid = bid | (blks.type != blks.type.shift(1))
            bid = bid | ((blks.foup - blks.fodown.shift(1) - 1) > min_block_distance)
            blks.rename({'block_id':'block'}, axis=1, inplace=True)
            blks['block_id'] = bid.cumsum()
            blks = blks.groupby([*cols,'block','block_id']).agg({
                'feature_order': lambda x: ", ".join(sorted(set([ str(y) for y in x]))),
                'internal_id': lambda x: ", ".join(sorted(set([ str(y) for y in x]))),
                'foup': 'min',
                'fodown': 'max',
                'is_fragment': 'all'
            }).reset_index()
            blks = blks.merge(dflim, on=['assembly','nucleotide','block','type'], how='left')
            blks.rename({'feature_order':'targets','internal_id':'tiids'}, axis=1, inplace=True)
            blks['origin'] = 0

            # Initiate analysis of circular replicons
            circular = blks.query('topology == "circular" and not is_fragment')
            if len(circular) > 0:
                circular = circular.groupby(['assembly','nucleotide','type'])
                circular = circular.agg({'block_id':['min','max'],'foup':'min','fodown':'max','fomin':'min','fomax':'max'})
                circular.reset_index(inplace=True)
                circular.columns = ['assembly','nucleotide','type','bidmin','bidmax','foup','fodown','fomin','fomax']

                # Merge first and last blocks of circular replicons, if too close
                tooclose = circular.eval(f'(fomax - fodown) + (foup - fomin) <= {min_block_distance}')
                if tooclose.any():
                    bidmaxIndex = blks[blks.block_id.isin(circular[tooclose].bidmax)].index
                    bidminIndex = blks[blks.block_id.isin(circular[tooclose].bidmin)].index
                    #logger.debug("First and last are too close:\n"+blks.iloc[bidminIndex.append(bidmaxIndex)].to_string()+"\n")
                    blks.loc[bidminIndex,'foup']     = blks.loc[bidminIndex,'fomin'] # Set foup   = fomin
                    blks.loc[bidmaxIndex,'fodown']   = blks.loc[bidmaxIndex,'fomax'] # Set fodown = fomax
                    blks.loc[bidmaxIndex,'block_id'] = blks.loc[bidminIndex,'block_id'].to_list() # Set block_ids to be the same
                    blks.loc[bidminIndex,'origin']   = 1
                    blks.loc[bidmaxIndex,'origin']   = 1
                    circular = circular[~tooclose]

                # Neighborhood runs through the entire replicon
                overrun = (circular.foup < circular.fomin) & (circular.fodown > circular.fomax)
                if overrun.any():
                    overIndex = blks[blks.block_id.isin(circular[overrun].bidmin)].index
                    #logger.debug("Overrun:\n"+blks.iloc[overIndex].to_string()+"\n")
                    blks.loc[overIndex,'foup']   = blks.loc[overIndex,'fomin'] # Set foup   = fomin
                    blks.loc[overIndex,'fodown'] = blks.loc[overIndex,'fomax'] # Set fodown = fomax
                    blks.loc[overIndex,'origin'] = 2
                    circular = circular[~overrun]

                # Neighborhood starts before the origin of replication
                startsBefore = (circular.foup < circular.fomin)
                if startsBefore.any():
                    beforeOriginIndex = blks[blks.block_id.isin(circular[startsBefore].bidmin)].index
                    blks.loc[beforeOriginIndex,'origin'] = 1
                    beforeOrigin = blks.iloc[beforeOriginIndex].copy()
                    #logger.debug("Starts before origin:\n"+beforeOrigin.to_string()+"\n")
                    beforeOrigin.foup   = beforeOrigin.fomax + beforeOrigin.foup - beforeOrigin.fomin + 1
                    beforeOrigin.fodown = beforeOrigin.fomax
                    blks.loc[beforeOriginIndex,'foup'] = blks.loc[beforeOriginIndex,'fomin'] # Set foup = fomin
                    blks = pd.concat([blks, beforeOrigin]).reset_index(drop=True)
                    circular = circular[~startsBefore]

                # Neighborhood ends after the origin of replication
                endsAfter = (circular.fodown > circular.fomax)
                if endsAfter.any():
                    afterOriginIndex = blks[blks.block_id.isin(circular[endsAfter].bidmax)].index
                    blks.loc[afterOriginIndex,'origin'] = 1
                    afterOrigin = blks.iloc[afterOriginIndex].copy()
                    #logger.debug("Ends after origin:\n"+afterOrigin.to_string()+"\n")
                    afterOrigin.foup   = afterOrigin.fomin
                    afterOrigin.fodown = afterOrigin.fomin + afterOrigin.fodown - afterOrigin.fomax - 1
                    blks.loc[afterOriginIndex,'fodown'] = blks.loc[afterOriginIndex,'fomax'] # Set fodown = fomax
                    blks = pd.concat([blks, afterOrigin]).reset_index(drop=True)

            # Truncate feature_id boundaries
            blks.foup   = np.where(blks.foup < blks.fomin, blks.fomin, blks.foup)
            blks.fodown = np.where(blks.fodown > blks.fomax, blks.fomax, blks.fodown)

            # Sort again and add row ID, i.e. within block row number
            blks.sort_values([*cols,'block_id','foup'], ascending=[True,True,True,True,True,False], inplace=True)
            blks.reset_index(inplace=True, drop=True)
            blks['rid'] = blks.index.to_series() - (blks.index.to_series() * (blks.block_id != blks.block_id.shift(1))).cummax()

            # Add first and last internal_id per block
            blks = blks.merge(
                    self[['assembly','nucleotide','type','feature_order','internal_id']],
                    left_on=['assembly','nucleotide','type','foup'],
                    right_on=['assembly','nucleotide','type','feature_order'],
                    how='left').rename({'internal_id':'up'}, axis=1).drop('feature_order', axis=1)
            blks = blks.merge(
                    self[['assembly','nucleotide','type','feature_order','internal_id']],
                    left_on=['assembly','nucleotide','type','fodown'],
                    right_on=['assembly','nucleotide','type','feature_order'],
                    how='left').rename({'internal_id':'down'}, axis=1).drop('feature_order', axis=1)

            # Adjust internal_id boundaries for blocks with all features of the same type as the query
            blks.up   = np.where(blks.foup   == blks.fomin, blks.iidmin, blks.up)
            blks.down = np.where(blks.fodown == blks.fomax, blks.iidmax, blks.down)
            blks.is_fragment = blks.is_fragment | (blks.up > blks.iidmin) | (blks.down < blks.iidmax)

        # If type is to be ignored, all analysis should focus on internal_ids
        elif (fttype == 'any'):
            logger.error(f'Unsuported fttype {fttype}')
            return pd.DataFrame()

        # Unkown fttype
        else:
            logger.error(f'Unknown fttype {fttype}')
            return pd.DataFrame()

        # Return a summary of all regions
        blks.drop_duplicates(inplace=True)
        return blks

    def neighbors(self, targets=['query == 1'], before=3, after=3, min_block_distance=0, strand='any', fttype='same'):
        """
        Find sets of rows, representing genomic regions, that are located near a set of targets.

        The user may choose a set of rows (targets) as anchors, whose neighbors will be
        evaluated by user-defined parameters, such as feature type, strand or distance.

        The returned value is copy of a slice of the original NeighborhoodDF
        object with updated block_id, origin, is_fragment and query columns.

        Arguments:
         - targets : (list of) boolean pd.Series or rules to select targets.

            Rules will used identify target rows using Pandas's eval method.

            If all rules return True for a given row, the row is a target and
            its neighbors will be identified.

         - before : keep at most this number of features, of the same type as the target,
                    before each target
         - after  : keep at most this number of features, of the same type as the target,
                    after each target

         - min_block_distance : minimum distance between two consecutive blocks
                                Blocks separated by more features 

         - strand : how to evaluate rows concerning the value of the strand column
                    Possible values for this option are:
                     - None or 'any' : ignore strand
                     - 'same' : keep only runs of same-strand features that
                                contain a target

         - fttype : how to process feature types of neighbors
                    Supported value:
                     - 'same' : consider only features of the same type as the target

        Raises ValueError for any other strand or fttype, or when a string
        rule in targets cannot be evaluated; TypeError when a target is
        neither a string nor a boolean Series.
        """
        _check_choice('strand', strand, [None, 'any', 'same'])
        _check_choice('fttype', fttype, ['same'])
        import sys
        import numpy as np
        from copy import deepcopy
        df = deepcopy(self)

        # Make sure data is sorted and has compatible internal ids
        df.sort_values(['assembly','nucleotide','block_id','internal_id'], ascending=[True,True,True,True], inplace=True)

        # Build boolean pandas.Series to mark targets
        select = True
        if not isinstance(targets, list):
            targets = [targets]
        for code in targets:
            if isinstance(code,str):
                try:
                    select &= df.eval(code)
                except Exception as exc:
                    raise ValueError(f'Rule {code!r} failed: {exc}') from exc
            elif isinstance(code,pd.Series):
                select &= code
            else:
                raise TypeError(f'targets must be strings or boolean Series, not {type(code).__name__}')

        # Find targets and their neighborhood
        if (not isinstance(select,pd.Series)) or (select.sum() == 0):
            logger.error(f'No anchors to search for neighbors were found! Revise your list of targets!')
            return seqrecords_to_dataframe([])

        # Initialize dataframe for each region (blocks)
        select = df.filter(['assembly','nucleotide','internal_id']).assign(query=select)
        blks = df.vicinity(select['query'], before, after, min_block_distance, fttype)
        if blks.empty:
            return seqrecords_to_dataframe([])

        # Expand list of features within each neighborhood, tag queries
        # IMPORTANT: row number may increase if min_block_distance < 0!
        blks['internal_id'] = list(map(lambda x: list(range(x[0], x[1] + 1)), blks[["up", "down"]].values.tolist()))
        blks = blks.filter(['assembly','nucleotide','internal_id','block_id','rid','origin','is_fragment']).explode('internal_id')
        blks = blks.merge(select.drop_duplicates(), on=['assembly','nucleotide','internal_id'], how="left")

        # Replace columns in self with those from blks
        dropList = ['query','block_id','origin','rid','is_fragment']
        cols = df.columns
        dropList = [ x for x in dropList if x in cols ]
        copy = df.drop(dropList, axis=1).sort_values(['assembly','nucleotide','start','end'])
        copy = copy.merge(blks, left_on=['assembly','nucleotide','internal_id'], right_on=['assembly','nucleotide','internal_id'], how='inner')
        copy.sort_values(['assembly','nucleotide','block_id','rid','internal_id'], ascending=[True,True,True,True,True], inplace=True)
        copy['query'] = np.where((~copy['query'].isna()) & copy['query'],1,0).astype(int)
        copy = copy[cols].copy()

        # Evaluate other criteria for defining blocks
        if strand == 'same':
            # Locate new blocks
            tmp  = (copy.assembly   == copy.shift(1).assembly)
            tmp &= (copy.nucleotide == copy.shift(1).nucleotide)
            tmp &= (copy.block_id   == copy.shift(1).block_id)
            tmp &= (copy.strand     == copy.shift(1).strand)
            copy['block_id'] = (~tmp).cumsum()

            # Cleanup
            valid = set(copy.groupby('block_id').agg({'query':'sum'}).query('query > 0').index)
            copy = copy.query('block_id in @valid').copy()
            copy['block_id'] = (~(copy.block_id == copy.shift(1).block_id)).cumsum()

        # Set within block row number
        copy.reset_index(inplace=True, drop=True)
        copy['internal_id'] = copy.index.to_series()

        # Convert block_ids to universal identifiers: nucleotide:start-end
        bid = copy.groupby('block_id').agg({'nucleotide':'first','start':'first','end':'last'})
        bid = bid.nucleotide + ":" + bid.start.astype(str) + "-" + bid.end.astype(str)
        bid = bid.to_dict()
        copy.block_id = copy.block_id.map(bid)

        # Return
        return copy

    def jaccard(self, column='c80e3', *, min_members=3, min_jaccard=0.3, seed=0):
        """
        Group neighbourhoods that share protein clusters.

        Each block is described by the set of distinct values of ``column``
        found in it, by default its ``c80e3`` protein clusters. Blocks are
        joined when the Jaccard index of their sets reaches ``min_jaccard``,
        and the graph made by those joins is split into Louvain communities
        [1]_ and into connected components.

        Parameters
        ----------
        column : str, default 'c80e3'
            Column whose distinct values describe a block. Missing values
            are ignored.
        min_members : int, default 3
            Blocks with fewer distinct values than this are kept in the
            output but left out of the graph, so their ``community`` and
            ``component`` are missing.
        min_jaccard : float, default 0.3
            Lowest Jaccard index, from 0 to 1, at which two blocks are
            joined. A pair exactly at the threshold is joined.
        seed : int or None, default 0
            Random seed for Louvain, so a rerun gives the same communities.
            None gives an unseeded run.

        Returns
        -------
        pandas.DataFrame
            One row per block, in the order in which blocks first appear in
            the table, with the columns:

            - ``block_id``: the block.
            - ``members``: the distinct values of ``column`` in the block,
              as a sorted tuple (empty when all are missing).
            - ``n_members``: the length of ``members``.
            - ``community``: the block's Louvain community (Int64).
            - ``component``: the block's connected component (Int64).

            Communities and components are numbered by size, 0 being the
            largest; ties go to the one whose earliest block comes first in
            the table. Both are missing (``<NA>``) for a block that is not in
            the graph: it has fewer than ``min_members`` values, or no other
            block reaches ``min_jaccard`` with it. ``self`` is not modified.

        Raises
        ------
        ImportError
            If :mod:`rotifer.graph.network` cannot be imported. Install
            the ``graph`` extra with ``pip install rotifer[graph]``.
        KeyError
            If ``block_id`` or ``column`` is missing.
        ValueError
            If ``min_members`` is lower than 1 or ``min_jaccard`` is not
            between 0 and 1.

        See Also
        --------
        rotifer.graph.network.overlap_edges : Computes the Jaccard index
            between blocks.
        rotifer.graph.network.get_network_community : Louvain and Leiden
            communities of a domain network.

        Notes
        -----
        Blocks with the same set of values are scored once, as one node,
        and all get that node's community and component. This keeps the
        number of pairs down when many blocks are identical, and it never
        changes the result, since identical blocks have a Jaccard index
        of 1.

        Louvain maximises modularity with a random node order, so another
        ``seed`` can give another partition on graphs with weak community
        structure. Components do not depend on the seed.

        References
        ----------
        .. [1] Blondel, V. D., Guillaume, J.-L., Lambiotte, R. & Lefebvre, E.
           "Fast unfolding of communities in large networks." J. Stat. Mech.
           P10008 (2008). https://doi.org/10.1088/1742-5468/2008/10/P10008

        Examples
        --------
        >>> ndf = NeighborhoodDF(pd.DataFrame({
        ...     'block_id': ['b1'] * 3 + ['b2'] * 3 + ['b3'] * 3 + ['b4'],
        ...     'c80e3': ['A', 'B', 'C', 'A', 'B', 'D', 'X', 'Y', 'Z', 'A']}))
        >>> ndf.jaccard()
          block_id    members  n_members  community  component
        0       b1  (A, B, C)          3          0          0
        1       b2  (A, B, D)          3          0          0
        2       b3  (X, Y, Z)          3       <NA>       <NA>
        3       b4       (A,)          1       <NA>       <NA>
        """
        try:
            from rotifer.graph.network import overlap_edges
            import networkx as nx
        except ImportError as exc:
            raise ImportError(
                "NeighborhoodDF.jaccard needs the graph extra: "
                "pip install rotifer[graph]") from exc

        missing = [c for c in ('block_id', column) if c not in self.columns]
        if missing:
            raise KeyError(f"columns missing from the table: {missing}")
        if min_members < 1:
            raise ValueError(f"min_members must be at least 1, got {min_members}")
        if not 0 <= min_jaccard <= 1:
            raise ValueError(f"min_jaccard must be between 0 and 1, got {min_jaccard}")

        cols = ['block_id', 'members', 'n_members', 'community', 'component']
        sets = {}
        for bid, values in self[['block_id', column]].groupby('block_id', sort=False)[column]:
            sets[bid] = tuple(sorted(values.dropna().unique()))
        out = pd.DataFrame({
            'block_id': list(sets),
            'members': pd.Series(list(sets.values()), dtype=object),
            'n_members': pd.Series([len(m) for m in sets.values()], dtype='int64'),
        })
        out['community'] = pd.array([pd.NA] * len(out), dtype='Int64')
        out['component'] = pd.array([pd.NA] * len(out), dtype='Int64')
        if out.empty:
            return out[cols]

        eligible = out.n_members >= min_members
        set_id = {}
        block_set = [None] * len(out)
        for i, (m, ok) in enumerate(zip(out.members, eligible)):
            if ok:
                block_set[i] = set_id.setdefault(m, len(set_id))
        if not set_id:
            return out[cols]
        long = pd.DataFrame(
            [(sid, v) for m, sid in set_id.items() for v in m],
            columns=['set_id', 'member'])
        edges = overlap_edges(long, 'set_id', 'member', metric='jaccard',
                              min_weight=min_jaccard)
        G = nx.Graph()
        G.add_weighted_edges_from(
            zip(edges.source, edges.target, edges.weight), weight='weight')
        # Blocks with identical sets have a Jaccard index of 1 with each
        # other, so a set held by two or more blocks is a node even when no
        # other set is joined to it.
        counts = pd.Series([s for s in block_set if s is not None]).value_counts()
        G.add_nodes_from(counts.index[counts > 1])
        if G.number_of_nodes() == 0:
            return out[cols]

        n_blocks, first = {}, {}
        for pos, sid in enumerate(block_set):
            if sid is not None and sid in G:
                n_blocks[sid] = n_blocks.get(sid, 0) + 1
                first.setdefault(sid, pos)

        def number(groups):
            keyed = []
            for g in groups:
                g = set(g)
                keyed.append((-sum(n_blocks[s] for s in g),
                              min(first[s] for s in g), g))
            keyed.sort(key=lambda t: t[:2])
            return {s: rank for rank, (_, _, g) in enumerate(keyed) for s in g}

        comm = number(nx.community.louvain_communities(G, weight='weight', seed=seed))
        comp = number(nx.connected_components(G))
        out['community'] = pd.array([comm.get(s, pd.NA) if s is not None else pd.NA
                                     for s in block_set], dtype='Int64')
        out['component'] = pd.array([comp.get(s, pd.NA) if s is not None else pd.NA
                                     for s in block_set], dtype='Int64')
        return out[cols]

    def to_blocks(self, **kwargs):
        """
        Render the neighborhood as text blocks, one per c80e3 cluster.

        Applies :data:`rotifer.genome.blocks.NEIGHBORHOOD_BLOCKS` to
        :func:`rotifer.pandas.to_blocks`. The preset groups by ``c80e3``,
        shows the columns ``pid``, ``c80i70``, ``c100i100``, ``plen``,
        ``pfam``, ``aravind``, ``product``, ``organism``, ``lineage`` and
        ``classification``, sorts the rows by ``lineage``,
        ``classification`` and ``organism``, and summarises each block by
        the number of distinct ``lineage``, ``block_id``, ``c80i70`` and
        ``pid`` values.

        Parameters
        ----------
        **kwargs
            Any argument of :func:`rotifer.pandas.to_blocks`. A keyword
            that is also a preset entry (``groupby``, ``columns``,
            ``stats``, ``sortrows``) replaces that entry.

        Returns
        -------
        str or None
            The text, or None when ``buf`` is given.

        See Also
        --------
        to_yaml : Same blocks as a YAML document of markdown tables.
        rotifer.pandas.to_blocks : The function this method calls.
        """
        import rotifer.pandas
        from rotifer.genome.blocks import NEIGHBORHOOD_BLOCKS
        return rotifer.pandas.to_blocks(self, **{**NEIGHBORHOOD_BLOCKS, **kwargs})

    def to_yaml(self, **kwargs):
        """
        Write the neighborhood as a YAML document of markdown tables.

        Applies :data:`rotifer.genome.blocks.NEIGHBORHOOD_YAML` to
        :func:`rotifer.pandas.yaml_blocks.to_yaml`. The preset groups by
        ``c80e3`` and shows the same columns as :meth:`to_blocks`. It sorts
        rows by ``superkingdom``, ``phylum``, ``lineage``,
        ``classification`` and ``c80i70`` in ascending order, and sorts
        blocks by ``phylum`` and ``c80i70`` in descending order. Each block gets an
        annotation table holding the most frequent ``pfam``, a ``stats``
        table with the number of distinct ``superkingdom``, ``lineage``,
        ``phylum``, ``taxid``, ``c80i70`` and ``pid`` values, and a
        ``dist`` table with the number of distinct ``c80i70`` per ``pfam``.

        Parameters
        ----------
        **kwargs
            Any argument of :func:`rotifer.pandas.yaml_blocks.to_yaml`. A
            keyword that is also a preset entry (``groupby``,
            ``columns``, ``sortrows``, ``sortblocks``, ``annotation``,
            ``stats``, ``dist``) replaces that entry.

        Returns
        -------
        str or None
            The text, or None when ``buf`` is given.

        See Also
        --------
        to_blocks : Plain text blocks.
        rotifer.pandas.yaml_blocks.read_yaml : Reads the document back into a frame.
        """
        import rotifer.pandas.yaml_blocks
        from rotifer.genome.blocks import NEIGHBORHOOD_YAML
        return rotifer.pandas.yaml_blocks.to_yaml(self, **{**NEIGHBORHOOD_YAML, **kwargs})

    def to_network(self, target='pfam', *, strand=True, ftype='CDS', interaction=True,
                   ignore=None, separator='+', separator_regex=False, replace=None,
                   replace_regex=False, stats=None):
        """
        List the edges between features that are adjacent in each block.

        Calls :func:`rotifer.graph.network.adjacency_edges` on the
        neighborhood table. With ``strand=True`` the table is first reduced
        to the query rows and the neighbours on the same strand, ordered so
        that every block reads in the direction of its query (blocks of a
        minus-strand query are reversed). With ``strand=False`` the table is
        used in its current row order, and neighbour edges are made
        undirected (``orient=True``).

        Parameters
        ----------
        target : str or list of str, default 'pfam'
            Column or columns that label each feature.
        strand : bool, default True
            Keep only features on the strand of the query and orient every
            block along it. With False, the frame is used in its current
            row order and neighbor edges are made undirected (the smaller
            label is ``source``); fusion edges keep the N to C order.
        ftype : str or list of str, default 'CDS'
            Values of the ``type`` column to keep.
        interaction : bool, default True
            Tell fusions from neighbours and count each kind separately.
        ignore : iterable of str, optional
            Labels to remove before the edges are made. A single string is
            one label.
        separator : str or None, default '+'
            Splits each label into one label per domain. None keeps each
            label whole.
        separator_regex : bool, default False
            Read ``separator`` as a regular expression.
        replace : mapping, optional
            Renames labels, after splitting.
        replace_regex : bool, default False
            Read the keys of ``replace`` as regular expressions.
        stats : dict, optional
            Named aggregations, as ``{name: (column, function)}``. Defaults
            to the number of rows (``weight``) and of distinct blocks
            (``blocks``) with each edge.

        Returns
        -------
        pandas.DataFrame
            One row per distinct edge, with the columns described in
            :func:`rotifer.graph.network.adjacency_edges`. ``self`` is not
            modified.

        Raises
        ------
        ImportError
            If :mod:`rotifer.graph.network` cannot be imported. Install
            the ``graph`` extra with ``pip install rotifer[graph]``.
        KeyError
            If ``block_id``, ``pid``, ``type`` or a ``target`` column is
            missing.

        See Also
        --------
        rotifer.graph.network.adjacency_edges : The function this method calls.
        rotifer.graph.network.get_relations : Relates every pair of domains in a
            block instead of adjacent ones.

        Notes
        -----
        The graph module is imported when the method is called, so
        :mod:`rotifer.genome.data` can be used without the ``graph`` extra.
        """
        try:
            from rotifer.graph.network import adjacency_edges
        except ImportError as exc:
            raise ImportError(
                "NeighborhoodDF.to_network needs the graph extra: pip install rotifer[graph]"
            ) from exc
        if strand:
            w = self.neighbors(before=len(self), after=len(self), strand='same')
            w['rid'] = range(1, len(w) + 1)
            w['rid'] = w['rid'] * w['strand']
            w = w.sort_values(['rid'])
        else:
            w = self
        return adjacency_edges(w, target, ftype=ftype, interaction=interaction,
                               ignore=ignore, separator=separator,
                               separator_regex=separator_regex, replace=replace,
                               replace_regex=replace_regex, orient=not strand,
                               stats=stats)
