#!/usr/bin/env python3
import os
import sys
import pandas as pd
from rotifer.core.functions import loadClasses
from rotifer.genome.data import NeighborhoodDF

# Core class
class database:
    '''
    What defines a genome block?
        - nucleotide
        - type
        - start
        - end
        - strand
        - feature_id
        - block_id
    '''

    def __init__(self, sources = [],
                 additional_sources = [],
                 verbose = 0, log_file = ''):
        '''
        '''

        # from rotifer.genome.ch import ClickHouse
        try:
            base = os.path.join(os.path.realpath(os.path.join(os.path.abspath(__file__), '..', '..') ), 'genome')
        except:
            base = os.path.dirname(os.path.realpath(__name__))

        self._verbose = verbose
        self._log_file = log_file
        self.sources =[]
        self.sources.extend(sources)
        self.config = {}


        self._switch_source = loadClasses( base+'.db' )

        if additional_sources:
            for additional_source in additional_sources:
                self._switch_source.update(loadClasses(additional_source))


    def open(self, source, config = {}):
        '''
        source: source type
        config: configuration
        '''
        self.target = self._switch_source[source](config,
                                             verbose = self._verbose,
                                             log_file = self._log_file)
        if source not in self.sources:
            self.sources.append(source)

        self.config.update({source: config})

    def submit(self, accs = '',
               input_type =['protein_acc'],
               **parameters):
        '''
        accs:       A list of Accession (check valid input types)
        input_type: Accession input type
                    - Accepted types:
                        - protein_acc
                        - locus
                        - assembly
                        - nucleotide
        block_id: Initial block id (default: 0)
        filterby: Several filters
        above:    How many CDS above (default: maximum)
        below:    How many CDS below (default: maximum)
        '''
        self._accs = accs
        self._input_type = input_type
        self._parameters = parameters
        self.target.submit(accs = self._accs,
                           input_type = input_type,
                           **parameters)
        pass

    def fetch_next(self, how = 'block'):
        # genome, block, nucleotide, assembly
        return self.target.fetch_next(how = how)

    def missing(self):
        return self.target.missing()

    def configuration(self):
        print(self.target.config)


if __name__ == '__main__':
    s = database(sources = ['clickhouse'], verbose = 3)

    config = {'clickhouse': {'uri': 'clickhouse://default:@localhost/rotifer',
                             'table_name': 'genomes'},
              'gff': {'b':2}
             }

    # Need to accept a list of modified
    # Need to filter by best

    for source in s.sources:
        print(config[source])
        s.open(source, config = config[source])
        s.submit(['WP_077633710.1'], input_type = ['protein_acc'], above = 9, below = 0,
                 verbose = 10, filterby = {'nucleotide': ['NZ_BHFM01000049.1'],
                                           'nuc_asm': ['GCF_900015815.1']})
        # s.submit(['WP_12749045.1','GCF_000067045.1','GCF_002094795.1', 'a'], block_id = 100, input_type = ['assembly', 'nucleotide', 'locus'], verbose = 3)
        # s.submit(['WP_127490459.1', 'AAP99025.1', 'AAP98158.1'], input_type = ['protein_accession'])
    # a = next(s.fetch_next('block'))
    # b = NeighborhoodDF(a)
    #
    # d = database('clickhouse', {1:2}, additional_source = '~/mymodules/test/dyn/dyn')
