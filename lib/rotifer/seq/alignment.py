#!/usr/bin/env python3

import os
import sys
# TESTING LINES

from subprocess import Popen, PIPE, STDOUT
from tempfile import mkstemp
import rotifer.core.functions as rcf
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import matplotlib.pyplot as plt
from rotifer.core.functions import loadClasses

from Bio import SeqIO
from io import StringIO as iStringIO

__version__ = 0.10
__authors__ = 'Gilberto Kaihami; Gianlucca Nicastro; Robson Souza'


class IO:
    def __init__(self, input_alignment, input_format = '', **kwargs):
        import warnings
        warnings.filterwarnings("ignore")
        self._pos = 0
        self._load_io()
        self.input_format = input_format
        self._kwargs = kwargs

        if isinstance(input_alignment, str):
            self.input_alignment = [iStringIO(input_alignment)]

        else:
            self.input_alignment = []
            for aln in input_alignment:
                if isinstance(aln, str) and not os.path.isfile(aln):
                    self.input_alignment.append(iStringIO(aln))
                else:
                    self.input_alignment.append(aln)

    def _load_io(self,sources = [],
                     additional_sources = [],
                     verbose = 0,
                     log_file = '',
                     **kwargs):

        try:
            base = os.path.join(os.path.realpath(os.path.join(os.path.abspath(__file__), '..', '..') ), 'seq')

        except:
            base = os.path.dirname(os.path.realpath(__name__))

        self._verbose = verbose
        self._log_file = log_file
        self.sources =[]
        self.sources.extend(sources)
        self.config = {}
        self._switch_source = loadClasses( base+'.alignment.io')
        if additional_sources:
            for additional_source in additional_sources:
                self._switch_source.update(loadClasses(additional_source))

    def __iter__(self):
        return self

    def __next__(self):
        try:
            alignment = self.input_alignment[self._pos]

        except IndexError:
            raise StopIteration
        self._pos += 1

        df, metadata = self._switch_source[self.input_format](alignment,
                                                 input_format = self.input_format,
                                                 **self._kwargs
                                                ).run()

        return MSA(df, metadata = metadata)

    def show_options(self):
        print(self._switch_source)

class MSA(pd.DataFrame):
    '''
    This class represents alignments as a table (pandas DataFrame)
    It is possible to use all power of numpy/pandas, like loc and iloc
    Also added some custom functions which includes:
        - slice
        - hydrophobic
        - plot_hydrophobic
        - consensus
        - conservation (including:
                            - Shannon entropy (classic)
                            - Jensen Shannon Divergence
                            - Property Entropy
                              )

    A sequence alignment need at least 2 features:
        - sequence_id
        - sequence

    These features are represented as metadata. If the user give this information
    it is possible to work with an n x m dimensional array, keeping all additional information (e.g. taxonomy).

    To overcome some issues here we decided to use MultiIndex.
    For example:

    Given a sequence fasta:
    >A
    PCK-ACG-
    >B
    PCK-ACG-
    >C
    ACT-RCG-
    >D
    ACTDRCG-

    It will be transformed as a MSA object.

    +-------------+-------------------------------+
    | sequence_id | sequence                      |
    | name        | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
    | A           | P | C | K | - | A | C | G | - |
    | B           | P | C | K | - | A | C | G | - |
    | C           | A | C | T | - | R | C | G | - |
    | D           | A | C | T | D | R | C | G | - |
    +-------------+-------------------------------+

    metadata = {'sequence_id': 'name',
                'sequence':'sequence'}

    Here the sequence_id and sequence (First row of the table) are the the first level MultiIndex

    This can be a major advantage since later it is possible to get some useful information from fragments of the df
    without compromissing the overall performance.

    For example given the following df:

    In the former example it is clear the advantages of using a dataframe.

    We can fast and simple calculate the conservation (entropy) for each clade:

    se_clade_1 = df[df['clade'] == 1].conservation()
    se_clade_2 = df[df['clade'] == 2].conservation()

    Since conservation method returns a numpy array it is possible to:
    delta_se = se_clade_1 - se_clade_2
    returning the difference between the entropy for each amino acid in each clade.

    TODO: Still need to add new methods for consevation, improve documentation and check if all math is right.
    TODO: The conservation algo can be improved.
    TODO: Insert metadata access when needed

    MAYBE: Add clustering?

    # Notes:
    # The hydrophobic graph and calculation ideia were derived from scikit-bio (https://github.com/biocore/scikit-bio)
    # The conservation measurements code were adapted from:
    #   Capra JA and Singh M.
    #   Predicting functionally important residues from sequence
    #   conservation. Bioinformatics. 23(15): 1875-1882, 2007

    '''
    _metadata = ['alignment',
                 'verbose',
                 'log_file',
                 'format',
                 'metadata',
                 'matrix',
                 'a',
                 'hydrophobic_matrix']

    def __init__(self, *args, **kwargs):
        import warnings
        warnings.filterwarnings("ignore")
        self.alignment = kwargs.pop('alignment', None)
        self.verbose = kwargs.pop('verbose', None)
        self.log_file = kwargs.pop('log_file', None)
        self.format = kwargs.pop('format', None)
        self.metadata = kwargs.pop('metadata', None)
        self.matrix = ''
        self.a = ''
        self.hydrophobic_matrix = ''

        super(MSA, self).__init__(*args, **kwargs)

    # Pandas contructor
    @property
    def _constructor(self):
        return MSA

    @classmethod
    def read(cls, input_alignment, input_format = 'fasta', **kwargs):
        '''
        Another constructor
        We need to decide if we are going to remove this method.
        Since we have aligment.IO
        '''

        try:
            input_file_checker = os.path.isfile(input_alignment)

        except:
            input_file_checker = False
        alignment = input_alignment

        if input_format in ['table', 'dataframe', 'df']:
            try:
                metadata = kwargs.pop('metadata', None)

            except:
                pass

            if input_file_checker:
                df = pd.read_csv(alignment, **kwargs)
                if 'sequence' in metadata.keys():
                    df[metadata['sequence']] = df[metadata['sequence']].str.ljust(df[metadata['sequence']].str.len().max(), '-')

                else:
                    df[df.columns[1]] = df[df.columns[1]].str.ljust(df[df.columns[1]].str.len().max(), '-')
                    metadata = {'sequence_id': df.columns[0],
                                'sequence': df.columns[1]
                                }

            else:
                pass

            ma = _seq2matrix(df[['sequence']], metadata = metadata)
            df = df.drop(columns = ['sequence'])

            tuple_columns = ()

            for col in df.columns:
                if col == 'sequence_id':
                    tuple_columns += ('_sequence_id', col)
                else:
                    tuple_columns += ('other_information', col)

            df.columns = pd.MultiIndex.from_tuples(
                                            [tuple_columns]
                                            )

            df = pd.concat([df,ma], axis = 1)

            return cls(df,columns = df.columns, metadata = metadata)

        else:

            idd=[]
            seq=[]

            if (input_file_checker):
                fasta_sequences = SeqIO.parse(alignment, input_format)

            else:
                fasta_sequences = SeqIO.parse(iStringIO(alignment), input_format)

            for fasta in fasta_sequences:
                idd.append(fasta.id)
                seq.append(str(fasta.seq))

            data={'sequence_id':idd, 'sequence':seq}

            df = pd.DataFrame(data=data)
            df['sequence'] = df['sequence'].str.ljust(df['sequence'].str.len().max(), '-')
            metadata = {'sequence_id': 'sequence_id',
                        'sequence': 'sequence'}

            ma = _seq2matrix(df[['sequence']], metadata = metadata)
            df = df.drop(columns = ['sequence'])

            tuple_columns = ()

            for col in df.columns:
                if col == 'sequence_id':
                    tuple_columns += ('_sequence_id', col)
                else:
                    tuple_columns += ('other_information', col)

            df.columns = pd.MultiIndex.from_tuples(
                                            [tuple_columns]
                                            )

            df = pd.concat([df,ma], axis = 1)

            return cls(df,columns = df.columns, metadata = metadata)


    def slice(self, seq_col = '', start = 0, end = 0, **kwargs):
        '''
        Slice the object, could be axis = 0 or axis = 1,
        TODO
        '''
        df = self
        if 'sequence' in self.metadata.keys():
            seq_col = self.metadata['sequence']

        else:
            seq_col = seq_col

        return df[seq_col].str[start:end]


    def hydrophobic(self, **kwargs):
        '''
        # Data derived from AAIndex:
        # http://www.genome.jp/dbget-bin/www_bget?aaindex:ARGP820101
        '''
        import numpy as np
        from collections import defaultdict

        df = self

        self._hydrophobicity_idx = defaultdict(lambda: np.nan)
        self._hydrophobicity_idx.update({'A': 0.61, 'L': 1.53, 'R': 0.60, 'K': 1.15, 'N': 0.06, 'M': 1.18,
                                  'D': 0.46, 'F': 2.02, 'C': 1.07, 'P': 1.95, 'Q': 0., 'S': 0.05,
                                  'E': 0.47, 'T': 0.05, 'G': 0.07, 'W': 2.65, 'H': 0.61, 'Y': 1.88,
                                  'I': 2.22, 'V': 1.32, '-': 1.15 })
        self._hydrophobicity_labels=['Hydrophilic', 'Medium', 'Hydrophobic']


        self.hydrophobic_matrix =  df['_rotifer.sequence'].apply(lambda x: [self._hydrophobicity_idx[y] for y in x])

        self.hydrophobic_matrix = self.hydrophobic_matrix.fillna(1.15)

    def plot_hydrophobic(self, figsize = (15,10),
                         cmap = sns.diverging_palette(240, 10, as_cmap=True),
                         dpi = 300,
                         cbar_orientation = 'horizontal',
                         cbar_fraction = 0.0046,
                         cbar_pad = 0.04,
                         consensus_level = [0.9],
                         slice = ()
                         ):
        from matplotlib import colors, cm

        self.hydrophobic()

        df = self

        fig, ax = plt.subplots(figsize = figsize, dpi = dpi)

        if slice:
            matrix = self.hydrophobic_matrix.loc[:, slice[0]:slice[1]]
        else:
            slice = (0, self.hydrophobic_matrix.shape[1])
            matrix = self.hydrophobic_matrix

        cax = ax.imshow(matrix,  cmap=cmap)

        values = list(self._hydrophobicity_idx.values())

        ax.set_xticks(range(df['_rotifer.sequence'].loc[:, slice[0]:slice[1]].shape[1]))

        # Need to fix this

        ax.set_yticklabels(df.ID.values, size=7)

        ax.set_xticklabels(self.consensus([0.9])[0.9][slice[0]:slice[1]], size=7)

        cbar = fig.colorbar(cax,
                            ticks=[min(values),
                                  np.nanmedian(values),
                                  max(values)],
                            orientation= cbar_orientation,
                            fraction=cbar_fraction, pad = cbar_pad
                            )

        if cbar_orientation == 'horizontal':
            cbar.ax.set_xticklabels(self._hydrophobicity_labels) # horizontal colorbar
        if cbar_orientation == 'vertical':
            cbar.ax.set_yticklabels(self._hydrophobicity_labels) # vertical colorbar

        fig.tight_layout()

        return ax

    def aa_frequency(self, gap = True, **kwargs):
        df = self

        dc_groups = {
            'a' : ['F','Y', 'W', 'H'],
            'l' : ['I','V','L'],
            'h' : ['F', 'Y', 'W', 'H' 'A', 'C', 'F', 'M', 'W', 'Y'],
            '+' : ['H', 'K', 'R'],
            '-' : [ 'D', 'E'],
            'c' : ['H', 'K', 'R', 'D', 'E'],
            'p' : ['H', 'K', 'R', 'D', 'E','Q', 'N', 'S', 'T','C'],
            'o' : ['S','T'],
            'u' : ['G', 'A', 'S'],
            's' : ['G', 'A', 'S','V', 'T', 'D', 'N', 'P', 'C'],
            'b' : ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        }

        all_aa = ['G','A','V','I','L','M','F','Y',
                  'W','H','C','P','K','R','D','E',
                  'Q','N','S','T', '-']


        residue_matrix = pd.DataFrame(index=all_aa)

        self.group_size = pd.DataFrame(columns = ['a','l','h',
                                                '+', '-', 'c',
                                                'p', 'o',
                                                'u', 's', 'b'])

        self.group_ranking = {k:v for k,v in  enumerate(['a', 'l', 'h', '+',
                                          '-', 'c', 'p', 'o',
                                          'u', 's', 'b'])}

        a = df['_rotifer.sequence'].apply(pd.Series.value_counts)

        if gap:
            residue_matrix = residue_matrix.join(a)
            residue_matrix = residue_matrix.rename(index={'-': '_'})

        else:
            residue_matrix = residue_matrix.join(a)
            residue_matrix = residue_matrix.rename(index={'-': '.'})

        for group, values in dc_groups.items():
            #group is g
            self.group_size[group] =df['_rotifer.sequence'].isin(values).sum()

        # Improve names
        self.a = residue_matrix/len(df['_rotifer.sequence'])
        self.c = self.group_size.T/len(df['_rotifer.sequence'])
        self.c = self.c.reset_index(drop = True)

    # def aa_frequency2(self, gap = True, **kwargs):
    #
    #     df = self
    #
    #     aromatic = ['F','Y', 'W', 'H']
    #
    #     alifatic = ['I','V','L']
    #     hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
    #
    #     positive = ['H', 'K', 'R']
    #     negative = [ 'D', 'E']
    #     charged = positive + negative
    #     polar = charged + ['Q', 'N', 'S', 'T','C']
    #
    #     alcohol = ['S','T']
    #
    #     tiny = ['G', 'A', 'S']
    #     small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
    #     big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
    #
    #     all_aa = ['G','A','V','I','L','M','F','Y',
    #               'W','H','C','P','K','R','D','E',
    #               'Q','N','S','T', '-']
    #
    #     self.group_size = pd.DataFrame.from_dict({'a': len(aromatic),
    #                                          'l': len(alifatic),
    #                                          'h': len(hydrophobic),
    #                                          '+': len(positive),
    #                                          '-': len(negative),
    #                                          'c': len(charged),
    #                                          'p': len(polar),
    #                                          'o': len(alcohol),
    #                                          'u': len(tiny),
    #                                          's': len(small),
    #                                          'b': len(big),
    #                                         }, orient='index').rename(
    #                                             {0:'size'}, axis=1)
    #
    #     self.group_ranking = {k:v for k,v in  enumerate(['a', 'l', 'h', '+',
    #                                       '-', 'c', 'p', 'o',
    #                                       'u', 's', 'b'])}
    #
    #
    #     residue_matrix = pd.DataFrame(index=all_aa)
    #
    #     group_matrix=pd.DataFrame(index=['a', 'l', 'h', '+',
    #                                  '-', 'c', 'p', 'o',
    #                                  'u', 's', 'b'])
    #
    #     for x in range(len(df['_rotifer.sequence'].T)):
    #         a = df['_rotifer.sequence'][x].value_counts().to_frame()
    #
    #         if gap:
    #             residue_matrix = residue_matrix.join(a.rename(index={'-': '_'}))
    #             residue_matrix = residue_matrix.rename(index={'-': '_'})
    #         else:
    #             residue_matrix = residue_matrix.join(a.rename(index={'-': '.'}))
    #             residue_matrix = residue_matrix.rename(index={'-': '.'})
    #
    #         d = {
    #             'a':np.where(df['_rotifer.sequence'][x].isin(aromatic),1,0).sum(),
    #             'l':np.where(df['_rotifer.sequence'][x].isin(alifatic),1,0).sum(),
    #             'h':np.where(df['_rotifer.sequence'][x].isin(hydrophobic),1,0).sum(),
    #             '+':np.where(df['_rotifer.sequence'][x].isin(positive),1,0).sum(),
    #             '-':np.where(df['_rotifer.sequence'][x].isin(negative),1,0).sum(),
    #             'c':np.where(df['_rotifer.sequence'][x].isin(charged),1,0).sum(),
    #             'p':np.where(df['_rotifer.sequence'][x].isin(polar),1,0).sum(),
    #             'o':np.where(df['_rotifer.sequence'][x].isin(alcohol),1,0).sum(),
    #             'u':np.where(df['_rotifer.sequence'][x].isin(tiny),1,0).sum(),
    #             's':np.where(df['_rotifer.sequence'][x].isin(small),1,0).sum(),
    #             'b':np.where(df['_rotifer.sequence'][x].isin(big),1,0).sum(),
    #             }
    #
    #         group_matrix = group_matrix.join(
    #             pd.DataFrame.from_dict(d,orient='index').rename({0:x}, axis=1)
    #                 )
    #
    #     # Improve names
    #     self.a = residue_matrix/len(df['_rotifer.sequence'])
    #     self.c = group_matrix/len(df['_rotifer.sequence'])
    #     self.c = self.c.reset_index(drop = True)
    #
    def consensus(self, thresholds,
                  output = 'dict',
                  **kwargs):
        '''
        output: dict/dataframe/matrix
        '''

        df = self

        if isinstance(self.a, pd.DataFrame):
            pass
        else:
            self.aa_frequency(**kwargs)

        dc = {}
        if not isinstance(thresholds, (list, tuple)):
            thresholds = [thresholds]
        for threshold in thresholds:
            residue_matrix = self.a.fillna(0)

            residue_res=  residue_matrix[residue_matrix >= threshold].idxmax().dropna()

            group_matrix = self.c.fillna(0)

            group_res = (group_matrix[group_matrix >= threshold])
            group_res = group_res.apply(lambda x: [idx for idx, v in enumerate(x) if pd.notnull(v) ]
                            ).apply(lambda x: min(x) if x else np.nan).dropna()

            # self.debug = group_res

            group_res = group_res.map(self.group_ranking)

            e = pd.DataFrame(index = [x for x in range(0, df['_rotifer.sequence'].shape[1])])

            e = e.join(pd.concat([residue_res, group_res], axis = 1 ))

            e = e.fillna(-999)

            e['gap'] = '.'

            try:
                res = np.where(e[0] != -999,
                            np.where(e[0] != '_',
                                       e[0],
                                       np.where((e[0] == '_') & (e[1] != -999),
                                                        e[1],
                                                         e[0]
                                                )
                                     ),
                           np.where(e[1] != -999, e[1], '.')
                           )

            except:
                res = '.' * df['_rotifer.sequence'].shape[1]

            dc[threshold] = ''.join(res)


        if output.lower() == 'dict':
            return dc

        elif output.lower() in ['df', 'dataframe', 'data frame']:
            return MSA(pd.DataFrame({'consensus_level': [x for x in dc.keys()],
                                          'sequence':[x for x in dc.values()]}),
                            metadata = self.metadata)

        elif output.lower() == 'matrix':
            return self.seq2matrix(df = pd.DataFrame({'consensus_level': [x for x in dc.keys()],
                                                 'sequence':[x for x in dc.values()]},
                                                ), no_self = True)


    def _load_stats(self,sources = [],
                     additional_sources = [],
                     verbose = 0,
                     log_file = '',
                     **kwargs):

        try:
            base = os.path.join(os.path.realpath(os.path.join(os.path.abspath(__file__), '..', '..') ), 'seq')
        except:
            base = os.path.dirname(os.path.realpath(__name__))
        self._verbose = verbose
        self._log_file = log_file
        self.sources =[]
        self.sources.extend(sources)
        self.config = {}
        self._switch_source = loadClasses( base+'.conservation')
        if additional_sources:
            for additional_source in additional_sources:
                self._switch_source.update(loadClasses(additional_source))

    def stats(self, metric = 'shannon_entropy',
                     gap_penalty = 1,
                     pseudocount = 0.0000001,
                     **kwargs
                     ):

        df = self
        self._load_stats()

        # Some shortcuts

        if metric == 'se':
            metric = 'shannon_entropy'
        elif metric == 'pe':
            metric = 'property_entropy'
        elif metric == 'jsd':
            metric = 'jensen_shannon_divergence'
        elif metric == 'pre':
            metric = 'property_relative_entropy'
        elif metric == 're':
            metric = 'relative_entropy'

        return self._switch_source[metric](df['_sequence_id'], df['_rotifer.sequence'],
                                           gap_penalty = gap_penalty,
                                           pseudocount = pseudocount,
                                           **kwargs).run()

    def conservation_types(self):
        self._load_stats()
        return list(self._switch_source.keys())

    def seq2colors(self):
        df = self
        colors, cmap, bounds, norm = _colors_dict_and_cmap()
        _ = list(df['_rotifer.sequence'].map(lambda x: [int(colors[y])for y in x]))
        return (_, cmap, bounds, norm)

    def plot_logo(self, font_family = 'Arial',
                        data_type='bits',
                        seq_type='dna',
                        yaxis='bits',
                        colorscheme='classic',
                        nrows=1,
                        # ncols=1,
                        padding=0,
                        draw_range=None,
                        coordinate_type='data',
                        draw_axis=False,
                        fontfamily='Arial',
                        debug=False,
                        ax=None,
                        dpi = 96):


        import math
        import seaborn
        import matplotlib.pyplot as plt
        import matplotlib.patheffects
        import matplotlib as mpl
        from matplotlib.font_manager import FontProperties
        from matplotlib import transforms

        df = self
        bits = self.stats('bit_score')

        # Fix this

def _seq2matrix(df, seq_col = '', **kwargs):

    try:
        metadata = kwargs['metadata']
    except:
        metadata = {}

    if 'sequence' in metadata.keys():
        seq_col = metadata['sequence']
    else:
        seq_col = seq_col
    df = df[seq_col].map(lambda x: list(x)).apply(pd.Series)
    df.columns = pd.MultiIndex.from_tuples([
                                ('_rotifer.sequence', col) for col in df.columns
                                        ]
                                    )
    return df

def _colors_dict_and_cmap():
    from matplotlib import colors as ccolors
    colors = {'A':1,
                    'I':1,
                    'L':1,
                    'M':1,
                    'F':1,
                    'W':1,
                    'V':1,
                    'K':2,
                    'R':2,
                    'E':3,
                    'D':3,
                    'N':4,
                    'Q':4,
                    'S':4,
                    'T':4,
                    'C':5,
                    'G':6,
                    'P':7,
                    'H':8,
                    'Y':8,
                    'X':0,
                    'a':8, #aromatic
                    'l':1, #alifatic
                    'h':1, #hydrophobic
                    '+':2, #positive
                    'c':8, #charged
                    'p':4 , #polar
                    'o':4, #alcohol
                    'u':1, #tiny
                    's':1, #small
                    'b':1, #big
                    '.':0, #gap
                    '-':3 #negative
                   }

    cmap = ccolors.ListedColormap(['white', (134/256,159/256,234/256), # blue 1
                                                              (220/256,54/256,35/256), # red 2
                                                              (178/256, 80/256,186/256),# magenta 3
                                                              (89/256, 188/256, 60/256), # green 4
                                                              (225/256, 134/256, 131/256), # pink 5
                                                              (227/256, 149/256, 87/256), # orange 6
                                                              (192/256, 192/256, 61/256), # yellow 7
                                                              (76/256, 161/256, 163/256) #cyan 8
                                                              #9 charged
                                                                                        #10 tiny
                                                                                        # 11 small
                               #                                                          # 12 big
                                                                                                                      ])
    bounds=[0,1,2,3,4,5,6,7,8]
    norm = ccolors.BoundaryNorm(bounds, cmap.N)
    return (colors, cmap, bounds, norm)

if __name__ == '__main__':
    pass
