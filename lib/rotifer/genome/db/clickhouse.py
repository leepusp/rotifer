#!/usr/bin/env python3
import pandas as pd
import rotifer
from rotifer.core.functions import loadConfig
logger = rotifer.logging.getLogger(__name__)

class clickhouse:
    def __init__(
            self,
            config= {},
            table_name = '',
            uri = '',
            ):

        '''
        Connect to clickhouse
        '''

        from rotifer.alchemy.connect import clickhouse

        logger.info(__name__)
        try:
            if config.startswith(':'):
                config = (loadConfig(config))
        except:
            pass

        self.config = config
        self.uri = 'clickhouse://default:@localhost/rotifer'
        self.table_name = 'genomes'

        for c in config.keys():
            if c in ['uri', 'table_name']:
                try:
                    setattr(self,c, config[c])
                except:
                    if self.uri == '' or self.table_name == '':
                        logger.warning(f'Argument not {c} is empty')
            else:
                logger.warning( f'Argument not used in the config[{c}]')
                logger.error(f'Loading CH configuration argument not used [{c}]')

        try:
            logger.error(f'Connecting to ClickHouse DB')
            db =  clickhouse(uri = self.uri,
                                    table_name = self.table_name)
            self._conn = db.conn
            self._db_columns = db.table_columns()
            logger.error(f'Connected to ClickHouse DB')
        except:
            logger.error(f'Unable to connect to ClickHouse Database')

    def submit(self, accs,
               input_type = ['protein_acc'],
               **parameters):

        '''
        Submit/commit a query
        ----------
        accs:       A list of accessions (check valid input types)
        input_type: Accession input type
                        - Accepted types:
                            - protein_acc
                            - locus
                            - assembly
                            - nucleotide
                            - organism
                            - gene

        parameters:
            block_id:
            filterby:
            above:
            below:

        '''

        _switch_ipttype = {
                           'protein_acc': 'pid',
                           'locus': 'locus',
                           'assembly': 'assembly',
                           'nucleotide': 'nucleotide',
                           'organism': 'organism',
                           'gene': 'gene'
                          }

        # Key is the column value is a list
        try:
            self._filterby = parameters['filterby']
        except:
            self._filterby = {}

        try:
            self._above = parameters['above']
        except:
            self._above = 'max'

        try:
            self._below = parameters['below']
        except:
            self._below = 'max'

        try:
            self._distance = parameters['distance']
        except:
            self._distance = 0

        try:
            self._block_id = parameters['block_id']
        except:
            self._block_id = 0

        # Abstraction of the input type
        # Convert an arbitrary input type to a valid ipt type
        self._input_type = []
        for ipt in input_type:
            try:
                self._input_type.append(_switch_ipttype[ipt])
            except:
                logger.warning( f'ClickHouse invalid input type [{ipt}]')
                logger.error(f'Invalid input type for ClickHouse [{ipt}], valid options [{", ".join(list(_switch_ipttype.keys()) )}]')
        self.accs = accs

        # Testing
        if accs == '':
            self.accs = ['XP_753543.1', 'OXN08776.1']
        else:
            self.accs = accs

        # self._input_type = inputtype

        self._nucleotide = self._find_nucletides()

    def fetch_all(self):
        # Get all hits
        pass

    def fetch_next(self, how, **kwargs):

        # genome, assembly, block
        # Generator
        try:
            _switch_fetch = {
                            'block': self._next_block,
                            'nucleotide': self._next_nucleotide,
                            'assembly': self._next_assembly,
                            'genome': self._next_assembly
                             }

            return _switch_fetch[how]()
        except:
            pass


    def configuration(self):
        # Show database configuration
        print(self.config)

    def missing(self):
        res = []
        for ipt in self._input_type:
            for x in range(0, len(self._input_type)):
                q_ipt_type = 'WHERE '
                s = self._conn.execute(f'''SELECT DISTINCT {ipt}
                                            FROM {self.table_name}
                                            WHERE {ipt} IN
                                             ({', '.join(self.accs_formated)})
                                             ''').fetchall()
                if s:
                    _ = list(zip(*s))[0]
                    res.extend(_)
                    break

        # List comprehension find missing values
        miss = [x for x in self.accs if x not in res]

        return miss

    def _find_nucletides(self):
        # Get assemblies here

        ### BLOCK FOR QUERY
        # Build query for type
        # Check if valid


        ### END BLOCK FOR QUERY
        # Best nucleotides
        # nuc_asm in nucs_found
        # self._filterby key is the column

        self.nucleotide_in_db = []
        n_partition = 2000
        accs_divided = [self.accs[x:x+n_partition] for x in range(0, len(self.accs), n_partition)]

        # Divide queris (2000)
        q_ipt_type = ''
        for x in range(0, len(self.accs), n_partition):
            self.accs_formated = ["'" + x + "'" for x in self.accs[x:x+n_partition]]

            for ipt_type in self._input_type:
                if ipt_type in self._db_columns:
                    q_ipt_type += f"{ipt_type} IN ({','.join(self.accs_formated)}) OR "
                else:
                    logger.warning(f'No column named [{ipt_type}]')

            self.q_ipt_type = q_ipt_type.rstrip(' OR ')
            print(self._filterby)
            try:
                if self._filterby:
                    filter_str = ''
                    for k in self._filterby.keys():
                        print(k)
                        if isinstance(self._filterby[k], (list,tuple)):
                            filter_ls = self._filterby[k]
                        else:
                            filter_ls = [self._filterby[k]]
                        print(filter_ls)
                        filter_str += f'''{k} IN ({','.join(["'"+str(x)+"'" for x in filter_ls])}) OR '''
                        print(filter_str)
                    print(filter_str)
                    filter_str = filter_str.rstrip(' OR ')
                    print(filter_str)
                    print('Query2')
                    print(f'''SELECT DISTINCT nucleotide, nuc_asm
                                                        FROM {self.table_name}
                                                        WHERE ({filter_str}) AND
                                                        ({self.q_ipt_type})
                                            ''')
                    # Build filter
                    _ = self._conn.execute(f'''SELECT DISTINCT nucleotide, nuc_asm
                                                        FROM {self.table_name}
                                                        WHERE ({self.q_ipt_type}) AND
                                                        ({filter_str})
                                            ''')

                else:
                    print(self.q_ipt_type)
                    print('Debug query')
                    print(f'''SELECT DISTINCT nucleotide, nuc_asm
                              FROM {self.table_name}
                              {self.q_ipt_type}''')
                    _ = self._conn.execute(f'''SELECT DISTINCT nucleotide, nuc_asm
                                                        FROM {self.table_name}
                                                        {self.q_ipt_type}
                                             ''')
                self.nucleotide_in_db.extend([(nuc, nuc_asm) for nuc, nuc_asm in _])
            except:
                logger.error(f'Issue with partition {x}')
        return self.nucleotide_in_db

    def _next_block(self):
        # Add filter by query
        # self.accs
        # self.q_ipt_type
        # self.q_ipt_type = q_ipt_type.rstrip(' OR ')


        # This is a list of tuples
        # original_acc, how to use
        for nucleotide,nuc_asm in self._nucleotide:

            # Filter here?

            try:
                accs_indexes = self._conn.execute(f"""SELECT feature_order from genomes
                                            where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                            ({self.q_ipt_type.replace('WHERE', '')})
                                            ORDER BY feature_order""").fetchall()
                accs_indexes= sorted(list(zip(*accs_indexes))[0])
                # accs_indexes= sorted(list(zip(*accs_indexes))[0])


                # accs_indexes = self._conn.execute(f"""SELECT feature_order FROM {self.table_name}
                #                              WHERE nuc_asm = '{nuc_asm}'
                #                                   AND nucleotide = '{nucleotide}'
                #                             {self.q_ipt_type.replace('WHERE','AND')}""")
                #

                max_index, min_index = self._conn.execute(f"""SELECT max(feature_order), min(feature_order) from {self.table_name}
                                                where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                                type = 'CDS'""").fetchone()

                max_internal_id, min_internal_id =  self._conn.execute(f"""SELECT max(internal_id), min(internal_id)
                                             FROM {self.table_name}
                                    where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}'""").fetchone()

                topology = self._conn.execute(f"""SELECT distinct(topology) from {self.table_name}
                                            where nuc_asm = '{nuc_asm}' and nucleotide = '{nucleotide}' and
                                            type = 'CDS'""").fetchone()[0]


                if self._below == 'max':
                    self._below = min_index
                if self._above == 'max':
                    self._above = max_index

                columns = ['internal_id','nucleotide',
                           'assembly', 'topology', 'start', 'end',
                           'strand','pid', 'type', 'plen', 'locus',
                           'seq_type', 'gene', 'product','organism',
                           'classification']


                intervals = []

                for value in accs_indexes:
                    find_up = value - self._above
                    find_down = value + self._below

                    if intervals:
                        if intervals[-1][1] - (find_up + 1) >= self._distance:
                            if (self._above + self._below + 1) >= max_index:
                                intervals[len(intervals)-1] = [0, max_index]
                            else:
                                intervals[len(intervals) -1 ] = [intervals[len(intervals)-1][0], find_down]
                        else:
                            if (self._above + self._below + 1) >= max_index:
                                intervals.append([0, max_index])
                            else:
                                intervals.append([find_up, find_down])

                    else:
                        if (self._above + self._below + 1) >= max_index:
                            intervals.append([0, max_index])
                        else:
                            intervals.append([find_up, find_down])

                for up_index, down_index in intervals:
                    sub_df = pd.DataFrame()
                    get_more_up = ''
                    get_more_down = ''
                    if up_index < min_index:
                        if topology != 'linear':
                            get_more_up = max_index+ up_index
                        up_index = min_index

                    if down_index > max_index:
                        if topology != 'linear':
                            get_more_down = down_index - max_index
                        down_index = max_index

                    if get_more_up:
                        _ = self._conn.execute(f"""
                                            SELECT internal_id from {self.table_name}
                                            WHERE nuc_asm = '{nuc_asm}' and
                                            nucleotide = '{nucleotide}' and
                                            type = 'CDS'
                                            and feature_order in ({max_index}, {get_more_up})
                                              """).fetchall()
                        _ = [x[0] for x in _]

                        if len (_) == 1:
                            max_order = min_order = _[0]
                            min_order = max_internal_id # check if exist

                        else:
                            max_order, min_order = sorted(_)
                            min_order = max_internal_id

                        # check this query
                        q = self._query_region(nuc_asm = nuc_asm,
                                               nucleotide = nucleotide,
                                               max_order = max_order,
                                               min_order = min_order)

                        sub_df = pd.DataFrame(q,
                                              columns = columns)

                    if not sub_df.empty:
                        _ = self._conn.execute(f"""
                                            SELECT internal_id from {self.table_name}
                                            WHERE nuc_asm = '{nuc_asm}' and
                                            nucleotide = '{nucleotide}' and
                                            type = 'CDS'
                                            and feature_order in ({up_index}, {down_index})
                                              """).fetchall()

                        _ = [x[0] for x in _]

                        if len (_) == 1:
                            max_order = min_order = _[0]
                        else:
                            max_order, min_order = sorted(_)

                        if get_more_up:
                            max_order = min_internal_id
                        if get_more_down:
                            min_order = max_internal_id

                        q = self._query_region(nuc_asm = nuc_asm,
                                               nucleotide = nucleotide,
                                               max_order = max_order,
                                               min_order = min_order)

                        t = pd.DataFrame(q,
                                              columns = columns)

                        sub_df = pd.concat([sub_df, t])

                    else:
                        _ = self._conn.execute(f"""
                                            SELECT internal_id from {self.table_name}
                                            WHERE nuc_asm = '{nuc_asm}' and
                                            nucleotide = '{nucleotide}' and
                                            type = 'CDS'
                                            and feature_order in ({up_index}, {down_index})
                                              """).fetchall()
                        _ = [x[0] for x in _]
                        if len (_) == 1:
                            max_order = min_order = _[0]
                        else:
                            max_order, min_order = sorted(_)

                        if get_more_up:
                            max_order = min_internal_id

                        if get_more_down:
                            min_order = max_internal_id

                        q = self._query_region(nuc_asm = nuc_asm,
                                               nucleotide = nucleotide,
                                               max_order = max_order,
                                               min_order = min_order)

                        sub_df = pd.DataFrame(q,
                                              columns = columns)

                    if get_more_down:
                        _ = self._conn.execute(f"""
                                            SELECT internal_id from {self.table_name}
                                            WHERE nuc_asm = '{nuc_asm}' and
                                            nucleotide = '{nucleotide}' and
                                            type = 'CDS'
                                            and feature_order in ({get_more_down})
                                              """).fetchall()
                        _ = [x[0] for x in _]
                        if len (_) == 1:
                            max_order = min_order = _[0]
                            max_order = min_internal_id
                        else:
                            max_order, min_order = sorted(_)
                            max_order = min_internal_id

                        q = self._query_region(nuc_asm = nuc_asm,
                                               nucleotide = nucleotide,
                                               max_order = max_order,
                                               min_order = min_order)

                        if not sub_df.empty:
                            t = pd.DataFrame(q,
                                              columns = columns)

                            sub_df = pd.concat([sub_df, t])


                    if not sub_df.empty:
                        self._block_id +=1
                        sub_df['block_id'] = sub_df.shape[0]*[self._block_id]

                        sub_df['query'] = sub_df['pid'].map(lambda x: 1 if x in self.accs else 0)
                        if self._block_id == 1:
                            print_header = True
                        else:
                            print_header = False

                        self.header = 'nucleotide start end strand block_id query pid type plen locus seq_type assembly gene modified product organism classification'.split()

                        try:
                            dc = {k:v for k,v in original_acc}
                            sub_df['modified'] = sub_df['pid'].map(lambda x: dc[x] if x in dc.keys() else '.')

                        except:
                            sub_df['modified'] = '.'

                        sub_df = sub_df.drop_duplicates()

                        yield sub_df[self.header] # this is a object for data

            except:
                pass
        self._nucleotide

    def _query_region(self, nuc_asm,
                      nucleotide,
                      max_order,
                      min_order):
        q = (self._conn.execute(f"""SELECT
                                      internal_id,nucleotide,
                                               assembly, topology, start, end,
                                               strand,pid, type, plen, locus,
                                               seq_type, gene, product,organism,
                                               taxonomy
                              FROM {self.table_name} where nuc_asm = '{nuc_asm}'
                                                  and nucleotide = '{nucleotide}'
                                                  and internal_id >= {max_order}
                                                  and internal_id <= {min_order}
                              ORDER BY internal_id""").fetchall())
        return q

#gdb.open(source = source, config = entrezopts)
#s.submit(accs, input_type = ['protein accession'], filterby=args.filterby, after=7, before=7)
# rgd = s.fetch_all()
# for rgd in s.fetch_next("genome"): # fetch next rotifer.genome.data (i.e. collection of dataframes)
# for block in rgd.blocks():
