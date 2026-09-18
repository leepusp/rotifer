#!/usr/bin/env python3

# Core functions for conservation calculation
import numpy as np
import pandas as pd
class conservation:
    # SE col, sim_matrix, bg_distr, seq_weights
    def __init__(self, aln, matrix, pseudocount = 0.0000001):
        self.aln = aln
        self.matrix = matrix
        self._seq_weights = self._calculate_sequence_weights()
        self.amino_acids = ['A', 'R', 'N', 'D',
                            'C', 'Q', 'E', 'G',
                            'H', 'I', 'L', 'K',
                            'M', 'F', 'P', 'S',
                            'T', 'W', 'Y', 'V',
                            '-']

        self.pseudocount = pseudocount

    def _weighted_freq_count_pseudocount_run(self, col):
        '''
        OK
        Intermediate function used in weighted_freq_count_pseudocount
        '''

    #    if the weights do not match, use equal weight
        if len(self._seq_weights) != len(col):
            self._seq_weights = [1.] * len(col)

        aa_num = 0
        freq_counts = len(self.amino_acids)*[self.pseudocount] # in order defined by amino_acids

        for aa in self.amino_acids:
            for j in range(len(col)):
                if col[j] == aa:
                    freq_counts[aa_num] += 1 * self._seq_weights[j]

            aa_num += 1

        for j in range(len(freq_counts)):
            freq_counts[j] = freq_counts[j] / (sum(self._seq_weights) + len(self.amino_acids) * self.pseudocount)

        return freq_counts

    def weighted_freq_count_pseudocount(self):
        """
        OK
        Return the weighted frequency count for a column--with pseudocount.
        OUTPUT:
        df: index is amino acid each column corresponds to the
        alignment position (0 index)
        """

        df = pd.DataFrame(index = self.amino_acids)

        for r in range(self.matrix.shape[1]):
            sub = pd.DataFrame(self._weighted_freq_count_pseudocount_run(self.matrix[r].values),
            index = self.amino_acids, columns = [r])

            df = df.join(sub)

        return df

    def _weighted_gap_penalty(self,col):
        """
        Using original
        Calculate the simple gap penalty multiplier for the column. If the
        sequences are weighted, the gaps, when penalized, are weighted
        accordingly. """

        # if the weights do not match, use equal weight
        if len(self._seq_weights) != self.aln.shape[0]:
            self._seq_weights = [1.] * self.aln.shape[0]

        gap_sum = 0.
        for i in range(len(self.matrix[col].values)):
            if self.matrix.iloc[i, col] == '-':
                gap_sum += self._seq_weights[i]

        return 1 - (gap_sum / sum(self._seq_weights))


    def _calculate_sequence_weights(self):
        """Ok
        Calculate the sequence weights using the Henikoff '94 method
        for the given msa.
        """

        counts = self.matrix.apply(pd.Series.value_counts).drop('-') # drop gaps

        observed = counts.apply(lambda x: sum([1 for y in x if pd.notnull(y)]), 0)

        seq_weights = []
        counts = self.matrix.apply(pd.Series.value_counts).drop('-')
        counts = counts.sort_index()

        observed = counts.apply(lambda x: sum([1 for y in x if pd.notnull(y)]), 0)
        observed = observed.sort_index()

        counts = counts.mul(observed.values, axis=1)
        counts = counts.T

        counts['-'] = np.nan

        for row in self.matrix.index:
            tmp = self.matrix.loc[row].to_frame().join(counts)
            tmp['look'] = 1/tmp.lookup(tmp.index, tmp[row])

            seq_weights.append(tmp['look'].sum()/self.matrix.shape[1])

        # seq_weights each element is the sequence weight same order in the rows
        return seq_weights


    # working here
    # Maybe it is possible to get better results
    # def _weighted_freq_count_pseudocount(aln, matrix, seq_weights = '', pc_amount = .0000001):
    #     ''
    #     aln alignment
    #     matrix: matrix
    #     col: each col of the matrix
    #     seq_weights = _calculate_sequence_weights otherwise = 1
    #     pc_amount = pseudocount
    #     ''
    #     # AA list including gaps
    #     amino_acids = ['A', 'R', 'N', 'D',
    #     'C', 'Q', 'E', 'G',
    #     'H', 'I', 'L', 'K',
    #     'M', 'F', 'P', 'S',
    #     'T', 'W', 'Y', 'V', '-']
    #
    #     pc_amount = .0000001
    #     if not seq_weights:
    #         seq_weights = _calculate_sequence_weights(aln, matrix)
    #
    #         # if the weights do not match, use equal weight
    #
    #     if len(seq_weights) != (aln.shape[0]):
    #         seq_weights = [1.] * (aln.shape[0])
    #
    #     d_final = pd.DataFrame()
    #
    #     for row in range(aln.shape[0]):
    #         #tmp = matrix.loc[row].to_frame()
    #         tmp = aln.matrix.loc[0].to_frame(name = 'aa')
    #         tmp['look'] = 1*seq_weights[row]
    #         #tmp['look'] = 1
    #         #tmp['look'] = tmp['look']*seq_weights[row]
    #
    #         #tmp = tmp[[row, 'look']]
    #         #tmp.columns = ['aa', 'look']
    #         tmp = tmp.reset_index()
    #
    #         d_final = pd.concat([d_final, tmp])
    #
    #     g = d_final.groupby(['index', 'aa']).agg({'look': 'sum'})
    #     g2 = g.reset_index().groupby('index')
    #     mock = pd.DataFrame(index = amino_acids)
    #     x = 0
    #     mock = g.reset_index().pivot(columns='aa', index = 'index', values = 'look').T
    #     return mock.fillna(pc_amount)
    #

    def _zscore(self, res):
        from scipy.stats import zscore
        return zscore(res)
