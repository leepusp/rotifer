#!/usr/bin/env python3
"""Frozen copy of NeighborhoodDF.compact_neighborhood as of #145.

test_compact_equivalence.py compares the current implementation with this
one. Only the Jupyter branch at the top was left out; the rest is verbatim.
Do not change it; remove it once the whole-table implementation (#144) has
been in use for a while.
"""

import numpy as np
import pandas as pd


def compact_reference(self,
                      annotation='pfam',
                      query='query',
                      strand=False,
                      fill_unk=False):
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

        return _.groupby(group).agg(agg_spec)

    # One row per group; also the schema of an empty result.
    agg_spec = {'compact': lambda x: "".join(list(x))[:-1],
                'block_id':'first',
                'lineage':'first',
                'organism':'first',
                'seq_type':'first',
                'classification':'first',
                'pid_compact': lambda x: "".join(list(x))[:-1]
               }

    # Main function, select if it will be grouped by block_id or strand,
    # and then make the groupby.
    group = 'block_id'
    if strand:
        group = 'strand_group'
        _1 = self.same_strand(query=query).copy()
    else:
        _1 = self.group_strand().copy()
    if _1.empty:
        return pd.DataFrame(columns=list(agg_spec)).set_index('block_id')

    # Iterate the groups instead of groupby().apply(): pandas 3 no longer
    # passes the grouping column to apply(), and main_cnei needs it.
    compacted_groups = [main_cnei(df_to_compact) for _, df_to_compact in _1.groupby(group)]
    return pd.concat(compacted_groups).reset_index(drop = True).set_index('block_id')
