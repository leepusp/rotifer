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

    def add_features(self, data, *args, inplace=False, **kwargs):
        # Prepare to copy data from self
        nucdata = ['nucleotide','taxonomy'] 
        nucdata = [ x for x in self.NDFProperties["columns"] if _column_dict[x] in nucdata ]
        nucdata = self.filter(nucdata).drop_duplicates()

        # Process input
        if not isinstance(data,list):
            data = [data]
        stack = []
        for datum in data:
            # Prepare input
            if not isinstance(datum, pd.DataFrame):
                datum = pd.DataFrame(datum, *args, **kwargs)
            else:
                datum = datum.copy()

            # Verify whether required columns are present and set
            cols = datum.columns.to_list()
            if not set(cols).issuperset(self.NDFProperties["required"]):
                missing = [ x for x in self.NDFProperties["required"] if x not in cols ]
                logger.error(f'DataFrame must contain all required columns: {self.NDFProperties["required"]}. Missing columns: {missing}')
                continue
            if self[self.NDFProperties["required"]].isna().any().any():
                missing = self[self.NDFProperties["required"]].isna().any().where(lambda x: x).dropna().index.to_list()
                logger.error(f'The following columns are not allowed to contain null values: {missing}')
                continue

            # Add default columns to input data
            tobeadded = nucdata.columns.difference(cols)
            if not tobeadded.empty:
                tobeadded = nucdata[tobeadded.insert(0,'nucleotide')].drop_duplicates()
                datum = datum.merge(tobeadded, on="nucleotide", how="left")
                cols = datum.columns.to_list()
            missing = [ x for x in self.NDFProperties["columns"] if x not in cols ]
            if missing:
                datum[missing] = np.nan
            if cols[0:len(self.NDFProperties["columns"])] != self.NDFProperties["columns"]:
                cols = self.NDFProperties["columns"] + [ x for x in cols if x not in self.NDFProperties["columns"] ]
                datum = datum[cols].copy()
            stack.append(datum)

        # Reset internal_id and feature_id
        if stack:
            stack = pd.concat([self] + stack, ignore_index=True)
            stack.sort_values(self.NDFProperties["order"], inplace=True)
            stack.reset_index(drop=True, inplace=True)
            stack['query'] = stack['query'].fillna(0).astype(int)
            stack['is_fragment'] = stack['is_fragment'].fillna(False).astype(bool)
            stack['block_id'] = np.where(stack.block_id.isna(), stack.nucleotide, stack.block_id)

            # Update internal_id
            stack.internal_id = (stack.assembly != stack.assembly.shift(1))
            stack.internal_id = stack.internal_id | (stack.nucleotide != stack.nucleotide.shift(1))
            stack.internal_id = (stack.internal_id * stack.index.to_series()).cummax()
            stack.internal_id = stack.index.to_series() - stack.internal_id

            # Update internal_id
            fomin = self.boundaries().filter(['assembly','nucleotide','type','fomin'])
            stack.sort_values(['type'] + self.NDFProperties["order"], inplace=True)
            stack.reset_index(drop=True, inplace=True)
            stack.feature_order = (stack.type != stack.type.shift(1))
            stack.feature_order = stack.feature_order | (stack.assembly != stack.assembly.shift(1))
            stack.feature_order = stack.feature_order | (stack.nucleotide != stack.nucleotide.shift(1))
            stack.feature_order = (stack.feature_order * stack.index.to_series()).cummax()
            stack.feature_order = stack.index.to_series() - stack.feature_order
            stack = stack.merge(fomin, on=['assembly','nucleotide','type'], how='left')
            stack.fomin.fillna(0, inplace=True)
            stack.eval('feature_order = feature_order + fomin', inplace=True)
            stack.drop('fomin', axis=1, inplace=True)
            stack.feature_order = stack.feature_order.astype(int)
            stack.sort_values(stack.NDFProperties["order"], inplace=True)

            # Return
            if inplace:
                self.drop(self.index, axis=0, inplace=True)
                self[stack.columns] = stack.values
            else:
                return stack

    def writer(self, output_format='table', **kwargs):
        '''
        kwargs:
          of | output_format: Output format (default: table)
           Supported output formats:
            - table   : pandas supported tabular output
            - compact : single row description of gene neighborhoods
            - gi2operons | gi2operon | gi2op : old TASS text view format
        '''

        if 'of' in kwargs.keys() or 'output_format' in kwargs.keys():
            try:
                self._of = kwargs['of'].lower()
            except:
                self._of = kwargs['output_format'].lower()

        if self._of.lower() in ['gi2operon', 'gi2operons', 'gi2op']:
            self._gi2operon(**kwargs)

    def _gi2operon(self, **kwargs):
        self._exclude_col = []
        self._include_col = []
        self._new_info_ls = []

        # Put df
        df = self.copy()
        df = df.reset_index(drop = True)
        df = df.astype(str)

        # Convert block id to numeric and other columns to string
        to_str = ['plen', 'start', 'end', 'strand', 'query']
        for col in to_str:
            try:
                df = df.astype({col: str})
            except: pass
        try:
            df['plen'] = df['plen'].map(lambda x: '.' if x == '' else x)
        except: pass

        if self._exclude_col:
            if 'gi' in self._exclude_col:
                exclude_col.append('modified')
            self._exclude_col = [x for x in self._exclude_col if x in df.columns and x not in self._include_col]
            df = df.drop(columns = self._exclude_col)

        if self._include_col:
            include = []
            for x in self._include_col:
                if x in df.columns:
                    include.append(x)
                else:
                    logger.warning(f"Column named [{columns}] was not included the DataFrame")
            df = df[include]

        try:
            df['pid'] = df['pid'].map(lambda x: '.' if x == '' else x)

        except:
            pass

        col_len = self._collen(df)

    def _collen(self, df):
        collen_len = []
        try: self.df['start'].astype(str)
        except: pass
        try: df['end'].astype(str)
        except: pass
        try: df['plen'].astype(str)
        except: pass

        q = '-->'
        collen_len.append(len(q)) # query mark

        try:
            df['cds_loc'] = sub_df['start'].astype(str)+'..'+sub_df['end'].astype(str)
        except: pass

        try:
            len_cds = df['cds_loc'].str.len().max()
            collen_len.append(len_cds)
        except: pass

        direction = len('dir')
        collen_len.append(direction)

        try:
            len_plen = df['plen'].astype(str).str.len().max() if df['plen'].astype(str).str.len().max() >= len('len') else len('len')
            collen_len.append(len_plen)
        except: pass

        try:
            len_pid = df['pid'].str.len().max() if df['pid'].str.len().max() >= len('pid') else len('pid')
            collen_len.append(len_pid)
        except: pass

        try:
            len_type = df['type'].str.len().max() if df['type'].str.len().max() >= len('type') else len('type')
            collen_len.append(len_type)

        except: pass

        try:
            len_gene = df['gene'].str.len().max() if df['gene'].str.len().max() >= len('gene') else len('gene')
            collen_len.append(len_gene)

        except:
            pass

        try:
            len_locus = df['locus'].str.len().max() if df['locus'].str.len().max() >= len('locus') else len('locus')
            collen_len.append(len_locus)
        except: pass

        try:
            try:
                len_modified = df['modified'].str.len().max() if df['modified'].str.len().max() >= len('gi') else len('gi')
                collen_len.append(len_modified)
            except:
                len_modified = 1
                collen_len.append(len_modified)
        except: pass

        try:
            len_product = df['product'].str.len().max() if df['product'].str.len().max() >= len('product') else len('product')
            collen_len.append(len_product)
        except: pass

        if self._new_info_ls:
            for col in self._new_info_ls:
                new_col = df[col].str.len().max() if df[col].str.len().max() >= len(col) else len(col)
                collen_len.append(new_col)

        return collen_len



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

    def series_to_compact_frequency(self, concat_string=": ", normalize=False):
        """
        Summarize the value counts of a series in a single string.

        Each distinct value is written followed by its count inside
        parentheses, e.g. ``AAA(3): ABC_tran(1)``.

        Parameters
        ----------
        concat_string : str, default ": "
            Separator placed between each ``value(count)`` item.
        normalize : bool, default False
            If True, report relative frequencies, formatted with two
            decimal places, instead of absolute counts.

        Returns
        -------
        str
            The ``value(count)`` items, sorted by decreasing count and
            joined by ``concat_string``.

        Notes
        -----
        This is meant to be called on a :class:`pandas.Series`, since it
        relies on :meth:`pandas.Series.value_counts`.
        """
        if normalize:
            count_dict = dict(
                self.value_counts(normalize=normalize).map('{:,.2f}'.format)
            )
        else:
            count_dict = dict(self.value_counts(normalize=normalize))
        count_list = []
        for each_series in count_dict:
            count_list.append(f'{each_series}({count_dict[each_series]})')
        return concat_string.join(count_list)

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

        def compacted(dataframe, annotation=annotation):
            """
            Take a groupby slice of the dataframe and make a new collum
            containing  a string of  annotation  and strand direction.
            """

            #dataframe[annotation].fillna('?', inplace=True) changed the line to be compatible with pandas 3.0
            dataframe.fillna({annotation:'?'}, inplace=True)
            dataframe[annotation] = np.where(
                dataframe['type'] != 'CDS',
                dataframe['type'],
                dataframe[annotation]
            )
            if fill_unk:
                dataframe[annotation] = np.where(
                    dataframe[annotation] == '?',
                    dataframe[fill_unk],
                    dataframe[annotation])

            dataframe[annotation] = dataframe[annotation].astype(str)
            return np.where(
                dataframe['strand_group'] == dataframe['strand_group'].shift(-1),
                np.where(
                    dataframe[query] == 0,
                    np.where(
                        dataframe.strand == 1,
                        dataframe[annotation]+'->',
                        '<-'+dataframe[annotation]),
                    np.where(
                        dataframe.strand == 1,
                        '*' + dataframe[annotation]+'->',
                        '<-'+dataframe[annotation]+'*')
                ),
                np.where(
                    dataframe[query] == 0,
                    np.where(
                        dataframe.strand == 1,
                        dataframe[annotation]+'->|',
                        '<-'+dataframe[annotation] + '|'),
                    np.where(
                        dataframe.strand == 1,
                        '*' + dataframe[annotation]+'->|',
                        '<-'+dataframe[annotation]+ '*|'),
                )
            )

        def main_cnei(df_to_compact):
            """
            Main function to make the compact_neighborhood.
            It was designed to be used in a groupby.
            """
            #Checking the orientation of the first query, and if it is
            # in the reverse strand, change the orientation. This way all
            # the blocks will be displayed at the same orientation.
            #Somehow if you are trying to groupby one block without a query
            # an error will be raised.
            try:
                if df_to_compact[df_to_compact[query] == 1].loc[:, 'strand'].iloc[0] == 1:
                    _ = df_to_compact.copy()
                else:
                    _ = df_to_compact.copy().iloc[::-1]
                    _.strand = _.strand *- 1
            except IndexError:
                raise NameError('Query_Missing')

            #Fill the unannotated with an ? marker or a column in the Dataframe
            # that the users can select by the option fill_unk=column name.
            _['compact'] = compacted(_)
            _['pid_compact'] = compacted(_, annotation='pid')

            return _.groupby(
                group
            ).agg({'compact': lambda x: "".join(list(x))[:-1],
                   'block_id':'first',
                   'lineage':'first',
                   'organism':'first',
                   'seq_type':'first',
                   'classification':'first',
                   'pid_compact': lambda x: "".join(list(x))[:-1]
                  })

        # Main function, select if it will be grouped by block_id or strand,
        # and then make the groupby.
        group = 'block_id'
        if strand:
            group = 'strand_group'
            _1 = self.same_strand(query=query).copy()
        else:
            _1 = self.group_strand().copy()
        return _1.groupby(group).apply(main_cnei).reset_index(drop = True).set_index('block_id')

    #Alias for compact neighborhood
    compact = compact_neighborhood
    def full_neighborhood(self, column='cluster', value='value', stats=False, new_query=True):
        """
        Function that filter a given value in a collumn and fetch
        all the loci within the block_id
        """
        selected_df = self[self.block_id.isin(
            self[self[column] == value].block_id
        )].copy()
        if stats:
            return  selected_df.value_counts_full(stats)

        selected_df.loc[:, 'new_query'] = 0
        selected_df.loc[selected_df[column] == value, 'new_query'] = 1
        selected_df = selected_df[selected_df.block_id.isin(
            selected_df[selected_df[column] == value].block_id)]
        if new_query:
            return selected_df
        return selected_df.drop('new_query', axis=1)

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
        targets : str, pandas.Series or list of these, default ['query == 1']
            Rules used to identify target rows.

            Strings are evaluated with :meth:`pandas.DataFrame.eval` and
            boolean :class:`pandas.Series` are used as given.

            If all rules return True for a given row, the row is a
            target and its neighbors will be identified.
        before : int, default 3
            Keep at most this number of features *before* each target.
        after : int, default 3
            Keep at most this number of features *after* each target.
        min_block_distance : int, default 0
            Minimum distance between two blocks.

            Blocks that are separated by a smaller number of features
            will be merged into a single block.
        fttype : {'same', 'any'}, default 'same'
            How to count neighbors by feature type.

            - ``'same'``: consider only features of the same type as the
              target.
            - ``'any'``: ignore feature type and count all neighboring
              features (not yet supported).

        Returns
        -------
        pandas.DataFrame
            One row per region, with its location (``assembly``,
            ``nucleotide``, ``type``, ``block_id``), the targets it
            contains (``targets``, ``tiids``), its feature order limits
            (``foup``, ``fodown``) and its first and last internal IDs
            (``up``, ``down``). The ``origin`` column is non-zero for
            regions that cross the origin of circular replicons.

            An empty DataFrame is returned when no target is found or
            ``fttype`` is not supported.
        """
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
                except:
                    logger.error(f'Rule {code} failed', file=sys.stderr)
            elif isinstance(code,pd.Series):
                select &= code

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
            logger.error(f'Unsuported fttype {fttype}', file=sys.stderr)
            return pd.DataFrame()

        # Unkown fttype
        else:
            logger.error(f'Unknown fttype {fttype}', file=sys.stderr)
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
                     - None : ignore strand
                     - same : same strand as the targets
                     -    + : positive strand features and targets only
                     -    - : negative strand features and targets only

         - fttype : how to process feature types of neighbors
                    Supported values:
                     - same : consider only features of the same type as the target
                     - any  : ignore feature type and count all features when
                              setting neighborhood boundaries

        """
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
                except:
                    logger.error(f'Rule {code} failed', file=sys.stderr)
            elif isinstance(code,pd.Series):
                select &= code

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

###### New functions added 03/12/2021 Gian
    def jaccard(self, min_c80e3=3):
        import community
        import networkx as nx
        df = self.copy()
        viz = df.groupby('block_id').agg(vizi = ('c80e3', lambda x: set(pd.Series(x).dropna().drop_duplicates().sort_values().to_list())), vc = ('c80e3', 'nunique')).reset_index()
        viz2 = viz.query('vc >= @min_c80e3')
        viz2 = viz2[~viz2.vizi.astype(str).duplicated(keep='first')]
        viz2['t'] = 1
        t2 = viz2.merge(viz2, on='t')
        t2.block_id_x, t2.block_id_y = np.where(t2.block_id_x > t2.block_id_y , [t2.block_id_x, t2.block_id_y], [t2.block_id_y, t2.block_id_x])
        t2 = t2.drop_duplicates(['block_id_x', 'block_id_y'])
        logger.info(f'tamanho t2 depois da ordenacao:  {len(t2)}')
        t2 = t2[t2.block_id_x != t2.block_id_y]
        t2['uni'] = t2.apply(lambda x: x.vizi_x.union(x.vizi_y), axis=1)
        t2['inter'] = t2.apply(lambda x: x.vizi_x.intersection(x.vizi_y), axis=1)
        t2['Jaccard_index'] = t2.apply(lambda x: len(x.inter)/len(x.uni), axis=1)
        t2 = t2.sort_values('Jaccard_index')
        t3 = t2[['block_id_x', 'block_id_y', 'Jaccard_index']].rename({'block_id_x': 'source', 'block_id_y':'target', 'Jaccard_index':'ji'}, axis=1)
        G2 = nx.from_pandas_edgelist(t3[['source', 'target', 'ji']].query('ji >= 0.3'), edge_attr='ji')
        partition = community.best_partition(G2,weight='ji')
        c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
            {'index': 'block_id', 0: 'nei_c'}, axis=1)
        viz.vizi = viz.vizi.astype('str')
        v3 = viz.merge(c)
        v3.vizi = v3.vizi.astype('str')
        cc = [c  for c in sorted(nx.connected_components(G2), key=len, reverse=True)]
        ncc = pd.Series(cc).reset_index().explode([0]).rename({'index':'ncc', 0: 'block_id'}, axis=1)
        v4 = v3.merge(ncc, how='left')
        v5 = viz.merge(v4[['vizi', 'nei_c', 'ncc']], how='left')
        return v5

