#!/usr/bin/env python3

import os
import sys
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkstemp
import rotifer.core.functions as rcf
import pandas as pd


class protein:
    '''
    Common protein sequence operations
    '''
    def __init__(self, *args, verbose = 0):
        '''
        Accepts one or more fasta sequences.
        Transform the input in a list, one element is one line in the fasta.
        --------
        EXAMPLE:
        protein1 = """>Protein1
        MATV"""

        protein2 = ">Protein2\nMATIU"

        file_protein = 'PATH_TO_FILE'

        protein(protein1) # add one fasta sequence

        protein(protein1, protein2) # Add two fasta sequence

        protein([protein1, protein2]) # Add a list containing two fasta sequence

        protein(protein1, file_protein) # Add one fasta string and one file containing fasta sequence

        protein(file_protein) # Add a file containing fasta sequence
        '''
        # check ipt
        self.fasta = rcf._fasta_ls(rcf._flatten(args))
        self.df = tb.fasta2df(self.fasta)

    def slice(self, start = 0, end = 1, inplace = False, col_name = ''):
        '''
        Slice sequence
        '''
        if inplace:
            if col_name:
               self.df[col_name]= self.df.Seq.str[start:end]
            else:
                self.df[str(len(self.df))] = self.df.Seq.str[start:end]

        else:
            return self.df.Seq.str[start:end]

    def seq_len(self, col = 'Seq'):
        '''
        Get the seq length
        -----------
        PARAMETERS:
        col: Column name to get the string length
        ----------
        RETURNS:
        The length of a sequence
        '''
        return self.df[col].str.len()

    def _tmp(self):
        fd, path = mkstemp()
        return path

