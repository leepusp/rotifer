#!/usr/bin/env python3

# Data here is a child of Pandas DataFrame

from collections import OrderedDict
import warnings
import numpy as np
import pandas as pd
import rotifer
import rotifer.pipeline
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
        dataframe['strand_group'], _ = self._strand_runs(query=None)
        return dataframe

    def same_strand(self,
                    query='query'):
        """
        Using the strand group and a query it will filter the self with
        strands that present one  query at least.
        """
        dataframe = self.group_strand()
        _, keep = dataframe._strand_runs(query=query)
        return dataframe[keep]

    def _strand_runs(self, group='block_id', query='query'):
        """
        Number the runs of same-strand rows and flag the runs with a query.

        A run is a stretch of consecutive rows (in table order) of the
        same ``group`` value and the same ``strand``. Shared by
        :meth:`group_strand`, :meth:`same_strand`,
        :meth:`trim_blocks` and :meth:`neighbors` (``strand='same'``).

        Parameters
        ----------
        group : str, default 'block_id'
            Column whose value change also starts a new run.
        query : str or None, default 'query'
            Column whose positive values mark the rows a run must
            contain to be kept. None keeps every row.

        Returns
        -------
        tuple of numpy.ndarray
            ``(run, keep)``: ``run`` is the int64 run number of each row,
            starting at 1 and increasing in table order; ``keep`` is
            True for the rows of runs whose ``query`` sum is above 0.
            Both are empty for an empty table.
        """
        if len(self) == 0:
            return np.empty(0, dtype=np.int64), np.empty(0, dtype=bool)
        block = self[group]
        same = (block == block.shift(1)) & (self.strand == self.strand.shift(1))
        run = np.cumsum(~same.to_numpy(), dtype=np.int64)
        if query is None:
            return run, np.ones(len(run), dtype=bool)
        total = self[query].groupby(run).transform('sum')
        return run, (total > 0).to_numpy()


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
        g_df = g_df.rename({'block_id':'block_id_count', 'pid':'pid_count'}, axis=1)
        logger.info(f'Total of block_id in selection:{lbid}')
        logger.info(f'Total of pid in selection :{lpid}')

        return g_df

    def trim_blocks(self, targets=['query == 1'], before=1, after=1, *, strand='any'):
        """
        Cut every block down to the rows around its targets.

        Works inside the blocks this table already has: for each target
        row, keep up to ``before`` rows above it and ``after`` rows below
        it in the same ``block_id``, counting every row whatever its
        feature type. Windows never cross into another block, and
        overlapping windows give each row once. Blocks without a target
        are dropped.

        Use it to shrink neighbourhoods you already have, for example
        to look at the three genes on each side of a protein cluster.
        To search a whole genome and build new blocks, counting only
        features of the target's type and wrapping around circular
        replicons, use :meth:`neighbors` instead.

        Parameters
        ----------
        targets : str, pandas.Series or list of them, default ['query == 1']
            Rules that select target rows; see :meth:`neighbors`. A
            single string is one rule.
        before : int, default 1
            Rows kept above each target (towards the block's first row).
            0 keeps none.
        after : int, default 1
            Rows kept below each target (towards the block's last row).
            0 keeps none.
        strand : {'any', 'same'}, default 'any'
            ``'same'`` then keeps only the runs of same-strand rows that
            contain a target (see :meth:`same_strand`); ``'any'``
            ignores strand.

        Returns
        -------
        NeighborhoodDF
            The kept rows in their original order, with their original
            index labels and ``block_id``. ``query`` is 1 for target rows
            and 0 for the others (int64). Empty, with the same columns,
            when no row is a target.

        Raises
        ------
        ValueError
            If ``before`` or ``after`` is negative, ``strand`` is not
            ``'any'`` or ``'same'``, or a rule cannot be evaluated.
        KeyError
            If the table has no ``block_id`` column.

        See Also
        --------
        neighbors : Find neighbourhoods along whole replicons.
        """
        _check_choice('strand', strand, ['any', 'same'])
        if before < 0 or after < 0:
            raise ValueError(f'before and after must not be negative; got before={before!r}, after={after!r}')
        if 'block_id' not in self.columns:
            raise KeyError('block_id')
        hit = self._select_targets(targets).to_numpy()
        if not hit.any():
            out = self.iloc[:0].copy()
            out['query'] = np.empty(0, dtype=np.int64)
            return out

        # Each target is a window of row positions inside its own block
        block = self['block_id']
        position = block.groupby(block, dropna=False).cumcount().to_numpy()
        windows = pd.DataFrame({
            'block_id': block.to_numpy()[hit],
            'up': position[hit] - before,
            'down': position[hit] + after,
        })
        windows['block'] = windows['block_id']
        work = self.assign(_trim_row=np.arange(len(self)), _trim_position=position)
        out = work._expand_windows(windows, '_trim_position', ['block_id'])

        # Back to table order, with query marking the targets
        out = out.iloc[np.argsort(out['_trim_row'].to_numpy(), kind='stable')]
        out['query'] = hit[out['_trim_row'].to_numpy()].astype(np.int64)
        out = out.drop(columns=['_trim_row', '_trim_position', 'block'])
        if strand == 'same':
            run, keep = out._strand_runs('block_id', 'query')
            out = out[keep]
        return out


    def select_neighbors(self, column='query', value=1, new_query=True, **kwargs):
        """
        Select the rows around ``column == value`` (deprecated).

        .. deprecated::
           Use :meth:`trim_blocks`. ``select_neighbors`` will be removed
           on 2026-11-30 (tracking issue #175).

        A shim over :meth:`trim_blocks` that keeps the old interface and
        output: ``query`` is left unchanged and the targets are marked in
        a ``new_query`` column instead, and ``strand=True`` adds the
        ``strand_group`` column that :meth:`same_strand` writes.

        Parameters
        ----------
        column : str, default 'query'
            Column compared with ``value`` to find the targets.
        value : scalar, default 1
            Value that marks a target row.
        new_query : bool, default True
            Keep the ``new_query`` column (1 for targets, 0 otherwise).
        **kwargs
            ``before`` and ``after`` (int, default 1), ``strand`` (bool,
            default False: True keeps only same-strand runs with a
            target) and ``stats`` (column name, default False: return
            :meth:`value_counts_full` of that column instead of rows).
            Other keywords raise TypeError.

        Returns
        -------
        NeighborhoodDF or pandas.DataFrame
            The selected rows, or the frequency table when ``stats`` is
            set.

        Warns
        -----
        FutureWarning
            On every call, pointing at the caller's line.
        """
        after = kwargs.pop('after', 1)
        before = kwargs.pop('before', 1)
        strand = kwargs.pop('strand', False)
        stats = kwargs.pop('stats', False)
        if kwargs:
            raise TypeError(f"select_neighbors() got an unexpected keyword argument {next(iter(kwargs))!r}")
        warnings.warn(
            "NeighborhoodDF.select_neighbors is deprecated and will be removed on 2026-11-30 "
            "(see #175); use NeighborhoodDF.trim_blocks",
            FutureWarning, stacklevel=2)

        # The old window ran `after` rows above and `before` rows below a target
        work = self.assign(_select_row=np.arange(len(self)))
        trimmed = work.trim_blocks([self[column] == value], before=after, after=before)
        selected_df = self.iloc[trimmed['_select_row'].to_numpy()].copy()
        selected_df['new_query'] = trimmed['query'].to_numpy()
        if strand:
            selected_df = selected_df.same_strand(query='new_query')
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

    def _select_targets(self, targets):
        """
        Mark the rows that anchor a neighbourhood.

        Shared by :meth:`neighbors` and :meth:`trim_blocks`. Each rule
        is evaluated on this table and the results are combined with a
        logical AND, so a row is a target only when every rule holds.

        Parameters
        ----------
        targets : str, pandas.Series or list of them
            Rules that select target rows. A string is evaluated with
            :meth:`pandas.DataFrame.eval` (for example ``"query == 1"``
            or ``"pid in ['WP_1.1', 'WP_2.1']"``); a single string is
            one rule, not a list of characters. A boolean Series is
            aligned on this table's index; labels it lacks count as
            False. A list combines its items.

        Returns
        -------
        pandas.Series
            Boolean, indexed like this table, named ``query``. All False
            when no row matches.

        Raises
        ------
        ValueError
            If a string rule cannot be evaluated (unknown column, syntax
            error) or does not give one boolean per row. The message
            names the rule.
        TypeError
            If an item is neither a string nor a boolean Series.
        """
        if not isinstance(targets, list):
            targets = [targets]
        select = pd.Series(True, index=self.index, name='query')
        for code in targets:
            if isinstance(code, str):
                try:
                    result = self.eval(code)
                except Exception as exc:
                    raise ValueError(f'Rule {code!r} failed: {exc}') from exc
                if not (isinstance(result, pd.Series) and len(result) == len(self)
                        and pd.api.types.is_bool_dtype(result.dtype)):
                    raise ValueError(f'Rule {code!r} does not give one boolean per row')
            elif isinstance(code, pd.Series):
                if not pd.api.types.is_bool_dtype(code.dtype):
                    raise TypeError(f'targets must be strings or boolean Series, not a Series of {code.dtype}')
                result = code
            else:
                raise TypeError(f'targets must be strings or boolean Series, not {type(code).__name__}')
            if not result.index.equals(self.index):
                result = result.reindex(self.index, fill_value=False)
            select &= result.fillna(False).astype(bool).to_numpy()
        return select


    def _expand_windows(self, windows, position, keys):
        """
        Return the rows of this table that fall inside closed windows.

        Shared by :meth:`neighbors` and :meth:`trim_blocks`. Each
        window covers the rows whose ``position`` value lies between its
        ``up`` and ``down`` values, both inclusive, among rows with the
        same ``keys`` values. Overlapping windows that share a
        ``block`` value give each row once.

        Parameters
        ----------
        windows : pandas.DataFrame
            One row per window, with the ``keys`` columns, ``up`` and
            ``down`` (integers, ``up <= down``) and ``block`` (the label
            given to the rows of that window). Other columns are copied
            to the output rows. It is not modified.
        position : str
            Column of this table holding each row's integer coordinate
            (``internal_id`` for :meth:`neighbors`, a row position for
            :meth:`trim_blocks`).
        keys : list of str
            Columns that must match between a window and its rows.

        Returns
        -------
        pandas.DataFrame
            The matching rows of this table, every column kept, plus the
            columns of ``windows`` that this table lacks. Rows are
            ordered by ``block`` and then by ``position``; the index is
            the original row label. Empty, with those columns, when no
            row matches.
        """
        w = windows.reset_index(drop=True)
        columns = [c for c in w.columns if c not in self.columns]
        pos = self[position].to_numpy().astype(np.int64)
        up = w['up'].to_numpy().astype(np.int64)
        down = w['down'].to_numpy().astype(np.int64)

        # One integer per key combination, so windows and rows can be sorted together
        both = pd.concat([self[keys], w[keys]], ignore_index=True)
        code = both.groupby(keys, sort=False, dropna=False).ngroup().to_numpy().astype(np.int64)
        row_code, win_code = code[:len(self)], code[len(self):]

        low = pos.min() if len(pos) else 0
        span = (pos.max() - low + 1) if len(pos) else 1
        row_key = row_code * span + (pos - low)
        order = np.argsort(row_key, kind='stable')
        sorted_key = row_key[order]
        lo = win_code * span + (np.clip(up, low, low + span - 1) - low)
        hi = win_code * span + (np.clip(down, low, low + span - 1) - low)
        start = np.searchsorted(sorted_key, lo, side='left')
        count = np.maximum(np.searchsorted(sorted_key, hi, side='right') - start, 0)
        count[(down < low) | (up > low + span - 1)] = 0

        # Rows of every window, then each row once per block
        which = np.repeat(np.arange(len(w)), count)
        offset = np.arange(count.sum()) - np.repeat(np.cumsum(count) - count, count)
        rows = order[start[which] + offset]
        label = pd.factorize(w['block'], sort=True)[0][which]
        first = ~pd.DataFrame({'block': label, 'row': rows}).duplicated().to_numpy()
        which, rows, label = which[first], rows[first], label[first]
        final = np.lexsort((np.arange(len(rows)), pos[rows], label))
        which, rows = which[final], rows[final]

        out = self.iloc[rows].copy()
        for name in columns:
            out[name] = w[name].to_numpy()[which]
        return out


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
        dflim = dflim.reset_index()
        dflim.columns = ['assembly','nucleotide','block','type','fomin','fomax']
        iidlim = df.groupby(['assembly','nucleotide','block_id']).agg({'internal_id':['min','max']})
        iidlim = iidlim.reset_index()
        iidlim.columns = ['assembly','nucleotide','block','iidmin','iidmax']
        bidlim = df.groupby(['assembly']).agg({'block_id':[lambda x: 1,'nunique']})
        bidlim = bidlim.reset_index()
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

        # Make sure data is sorted and has compatible internal ids
        blks = self.sort_values(['assembly','nucleotide','block_id','internal_id'], ascending=[True,True,True,True])

        # Boolean pandas.Series that marks targets
        select = blks._select_targets(targets)

        # Find targets and their neighborhood
        if not select.any():
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
            blks = blks.sort_values([*cols,'feature_order'])
            bid = (blks.nucleotide != blks.nucleotide.shift(1))
            bid = bid | (blks.is_fragment & (blks.block_id != blks.block_id.shift(1)))
            bid = bid | (blks.type != blks.type.shift(1))
            bid = bid | ((blks.foup - blks.fodown.shift(1) - 1) > min_block_distance)
            blks = blks.rename({'block_id':'block'}, axis=1)
            blks['block_id'] = bid.cumsum()
            blks = blks.groupby([*cols,'block','block_id']).agg({
                'feature_order': lambda x: ", ".join(sorted(set([ str(y) for y in x]))),
                'internal_id': lambda x: ", ".join(sorted(set([ str(y) for y in x]))),
                'foup': 'min',
                'fodown': 'max',
                'is_fragment': 'all'
            }).reset_index()
            blks = blks.merge(dflim, on=['assembly','nucleotide','block','type'], how='left')
            blks = blks.rename({'feature_order':'targets','internal_id':'tiids'}, axis=1)
            blks['origin'] = 0

            # Initiate analysis of circular replicons
            circular = blks.query('topology == "circular" and not is_fragment')
            if len(circular) > 0:
                circular = circular.groupby(['assembly','nucleotide','type'])
                circular = circular.agg({'block_id':['min','max'],'foup':'min','fodown':'max','fomin':'min','fomax':'max'})
                circular = circular.reset_index()
                circular.columns = ['assembly','nucleotide','type','bidmin','bidmax','foup','fodown','fomin','fomax']

                # Merge first and last blocks of circular replicons, if too close
                tooclose = circular.eval(f'(fomax - fodown) + (foup - fomin) <= {min_block_distance}')
                if tooclose.any():
                    bidmaxIndex = blks[blks.block_id.isin(circular[tooclose].bidmax)].index
                    bidminIndex = blks[blks.block_id.isin(circular[tooclose].bidmin)].index
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
                    afterOrigin.foup   = afterOrigin.fomin
                    afterOrigin.fodown = afterOrigin.fomin + afterOrigin.fodown - afterOrigin.fomax - 1
                    blks.loc[afterOriginIndex,'fodown'] = blks.loc[afterOriginIndex,'fomax'] # Set fodown = fomax
                    blks = pd.concat([blks, afterOrigin]).reset_index(drop=True)

            # Truncate feature_id boundaries
            blks.foup   = np.where(blks.foup < blks.fomin, blks.fomin, blks.foup)
            blks.fodown = np.where(blks.fodown > blks.fomax, blks.fomax, blks.fodown)

            # Sort again and add row ID, i.e. within block row number
            blks = blks.sort_values([*cols,'block_id','foup'], ascending=[True,True,True,True,True,False])
            blks = blks.reset_index(drop=True)
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

        # Return a summary of all regions
        return blks.drop_duplicates()

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
        # Make sure data is sorted and has compatible internal ids
        df = self.sort_values(['assembly','nucleotide','block_id','internal_id'], ascending=[True,True,True,True])

        # Boolean pandas.Series that marks targets
        select = df._select_targets(targets)

        # Find targets and their neighborhood
        if not select.any():
            logger.error('No anchors to search for neighbors were found! Revise your list of targets!')
            return seqrecords_to_dataframe([])

        # Initialize dataframe for each region (blocks)
        blks = df.vicinity(select, before, after, min_block_distance, fttype)
        if blks.empty:
            return seqrecords_to_dataframe([])

        # Expand list of features within each neighborhood, tag queries
        # IMPORTANT: row number may increase if min_block_distance < 0!
        windows = blks.filter(['assembly','nucleotide','block_id','rid','origin','is_fragment','up','down'])
        windows['block'] = windows['block_id']
        dropList = [ x for x in ['query','block_id','origin','rid','is_fragment'] if x in df.columns ]
        # Rows that share an internal_id (a CDS split at the origin) keep their start, end order
        copy = df.drop(dropList, axis=1).sort_values(['assembly','nucleotide','start','end'])
        copy = copy._expand_windows(windows, 'internal_id', ['assembly','nucleotide'])
        # Tag queries. A CDS split at the origin has two rows with one internal_id; the
        # merge keeps a row for each distinct query value, as it always has
        marks = df.filter(['assembly','nucleotide','internal_id']).assign(query=select).drop_duplicates()
        copy = copy.merge(marks, on=['assembly','nucleotide','internal_id'], how='left')
        copy['query'] = np.where((~copy['query'].isna()) & copy['query'], 1, 0).astype(int)
        copy = copy.sort_values(['assembly','nucleotide','block_id','rid','internal_id'], ascending=[True,True,True,True,True])
        copy = copy[df.columns].copy()

        # Evaluate other criteria for defining blocks
        if strand == 'same':
            # Keep runs of same-strand rows that contain a target; they are the new blocks
            run, keep = copy._strand_runs('block_id', 'query')
            run = run[keep]
            copy = copy[keep].copy()
            copy['block_id'] = np.cumsum(np.r_[True, run[1:] != run[:-1]]) if len(run) else run

        # Set within block row number
        copy = copy.reset_index(drop=True)
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
