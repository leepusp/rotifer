#!/usr/bin/env python3
import os
import sys
from subprocess import Popen, PIPE, STDOUT
from multiprocessing import Process, Manager
from tempfile import mkstemp
import rotifer.core.functions as rcf
import rotifer.table.table as tb
import pandas as pd
import numpy as np
import copy
from pandas.compat import StringIO
import matplotlib.pyplot as plt
import seaborn as sns
from clickhouse_driver import Client

class search:
    def __init__(self, *args, methods = {}):
        '''
        The search class accepts as input a fasta file or fasta sequence.
        -----------
        PARAMETERS:
        args:     A list of fasta file or fasta sequence (string/file)
        methods:: Set models and parameters              (dict)
        -----------
        ATTRIBUTES:
        fasta:   A list with each line one element of the list
        df:      A dataframe containg the Header, ID and Sequence
        -----------
        DETAILS:
        methods:
          - hmmscan
              - threads/thread: number of threads (default: 8)
              - db/database:      model database (default: :db.models:pfam)
              - of/output_format:  output format ('stdout', 'string', 'list') (default: 'list')
        '''

        self.fasta = rcf._fasta_ls(rcf._flatten(args))
#        methods = ['hmmscan', 'phobius', 'tmhmm', 'rpsblast',
#                   'blast', 'mmseqs', hmmsearch]
#        PARAMETERS: DB
#                    threads
        _switch_case = {'hmmscan': self._hmmscan}
        # self.hmmscan = ''

        #                 'hmmsearch': self._hmmsearch,
        #                 'phobius': self._phobius,
        #                 'rpsblast': self._rpsblast,
        #                 'blast': self._blast,
        #                 'mmseqs': self._mmseqs
        #                 }
        for method in methods.keys():
            _switch_case[method](methods[method])

    def _hmmscan(self, kw):
        '''
        Run hmmscan
        -----------
        PARAMETERS:
        kw: a dictionary containing the parameters:
            - threads       (default: 8)
            - database      (default: :db.models:pfam)
            - output format (default: list) (options: string, list, or stdout)
        -----------
        ATTRIBUTES:
        hmmscan: a list or string containing hmmscan result
        '''
        threads = 8
        db = ':db.models.pfam'
        output = 'list'

        if 'db' in kw.keys() or 'database' in kw.keys():
            try:
                db = kw['db']
            except:
                db = kw['database']
        else:
            pass

        db = self._dbfinder(db)

        if 'threads' in kw.keys() or 'thread' in kw.keys():
            try:
                threads = kw['threads']
            except:
                threads = kw['thread']

        if 'of' in kw.keys() or 'output_format' in kw.keys():
            try:
                output = kw['of']
            except:
                output = kw['output_format']

        ####################
        ### Multiprocess ###
        ####################

        # Create a manager for multiprocessing output
        manager = Manager()
        q = manager.list()
        _hmmscan_list = manager.list()

        # Split sequences in sub sequences
        df = tb.fasta2df(self.fasta)
        df['2write'] = '>'+df['ID']+'\n'+df['Seq']

        n  = (df['2write'].shape[0])// threads if (df['2write'].shape[0])//threads > 0 else 1

        sub_fasta = [df['2write'].values[x:x+n] for x in range(0, df['2write'].shape[0], n)
                     ]
        jobs = []

        for x in range(0, len(sub_fasta)):
            p = Process(target = self._multi_hmmscan, args = (sub_fasta[x], db, output, _hmmscan_list))
            p.start()
            jobs.append(p)

        try:
            for p in jobs:
                p.join()
        except KeyboardInterrupt:
            for p in jobs:
                p.terminate()
            raise error.with_traceback(sys.exc_info()[2])

        ########################
        ### End multiprocess ###
        ########################

        # Decide output format [string or list]
        if output == 'string':
            self.hmmscan = ''.join([x for x in _hmmscan_list])

            return self.hmmscan
        elif output == 'list':
            self.hmmscan = [x for x in _hmmscan_list]
            return self.hmmscan
        else:
            pass

    def _multi_hmmscan(self, fasta, db,output, ls):
        '''
        A function for multiprocessing hmmscan
        ----------
        PARAMETERS:
        fasta: a sub sequence fasta (list)
        db:    the database to search
        ls:    append result to a list
        '''
        try:
            path_fasta = self._tmp()
            with open(path_fasta, 'a') as f:
                f.write('\n'.join(fasta))

            cmd = '''hmmscan {db} {path_fasta}'''.format(path_fasta = path_fasta, db = db)

            p = Popen(cmd, shell = True, stdout = PIPE, stderr = PIPE)

            arch = p.communicate()[0].decode('utf-8')
            if output == 'stdout':
                print(arch.rstrip('\n'))
            else:
                ls.append(arch)
        except KeyboardInterrupt:
            raise error.with_traceback(sys.exc_info()[2])
        finally:
            os.remove(path_fasta)

    def parser(self, *args, **kwargs):
        '''
        Parse search results

        ---------
        PARAMETERS:
        args : a List of what to parse (not working)
        kwargs: a kwarg containing a dictionary with parameters
        ---------
        ATTRIBUTES:
        *_parsed: the parsed result (where *: is the method)
        ---------
        EXAMPLE:
        parser(hmmscan = {Parameters})

        to_search = search(fasta_file, methods = {'hmmscan': {} })
        to_search.parser(hmmscan = {'threads': 8,
                                    'of': 'raw'})
        # Access result using:
        to_search.hmmscan_parsed
        '''

        _switch_parser = {'hmmscan': self._parse_hmmscan}
        _vars = copy.deepcopy(list(self._vars()))

        self.a = ''
        for var in _vars:
            if var in _switch_parser.keys():
                try:
                    _switch_parser[var](kwargs[var])

                except:
                    pass

    def _parse_hmmscan(self, kw):
        '''
        parse hmmscan result

        Ouput format:
            - raw/table
            - arch/architecture/architecture2table
            - dom/domain/domain2architecture
        '''
        threads = 8
        options = ('raw', 'table',
                   'architecture', 'arch',
                   'domain', 'dom')

        option = 'raw'

        if 'thread' in kw.keys() or 'threads' in kw.keys():
            try:
                threads = kw['thread']
            except:
                threads = kw['threads']

        if 'of' in kw.keys() or 'output_format' in kw.keys():
            try:
                option = kw['of']
            except:
                option = kw['output_format']

        if 'insert2db' in kw.keys():
            insertdb = kw['insert2db']
        else:
            insertdb = False

        if 'user' in kw.keys():
            user = kw['user']
        else:
            user = 'rotifer'

        if 'table' in kw.keys():
            table = kw['table']

        else:
            table = ''

        # User = rotifer
        # table = hmmscan

        ####################
        ### Multiprocess ###
        ####################

        # Create a manager for multiprocessing output

        # Split sequences in sub sequences
        if isinstance(self.hmmscan, list):
            manager = Manager()
            q = manager.list()
            _hmmscan_parsed_list = manager.list()
            n = len(self.hmmscan)//threads if len(self.hmmscan)//threads >0 else 1
            sub_list_hmmscan = [self.hmmscan[x:x+n] for x in range(0, len(self.hmmscan), n)]

            jobs = []
            for x in range(0, len(sub_list_hmmscan)):
                p = Process(target = self._multi_parse_hmmscan, args = (sub_list_hmmscan[x],
                                                                        option,
                                                                        _hmmscan_parsed_list,
                                                                        insertdb,
                                                                        user,
                                                                        table)
                            )

                p.start()
                jobs.append(p)
            try:
                for p in jobs:
                    p.join()
            except KeyboardInterrupt:
                for p in jobs:
                    p.terminate()
                raise error.with_traceback(sys.exc_info()[2])
            df = pd.DataFrame()
            for result in _hmmscan_parsed_list:
                if df.empty:
                    df = result
                else:
                    df = pd.concat([df, result])
            self.hmmscan_parsed = df
        ########################
        ### End multiprocess ###
        ########################

        else:
            # Single hmmscan parser
            self.hmmscan_parsed =  self._single_parse_hmmscan(option)
        return self.hmmscan_parsed

    def _single_parse_hmmscan(self, option):
        '''
        Parse hmmscan one thread
        '''
        try:
            path_hmmscan = self._tmp()
            with open(path_hmmscan, 'a') as f:
                f.write(self.hmmscan)

            if option in ['raw', 'table']:
                cmd = '''hmmer2table -s -y {path_hmmscan}'''.format(path_hmmscan = path_hmmscan)
            elif option in ['architecture', 'architecture2table', 'arch']:
                cmd = '''hmmer2table -s {path_hmmscan} | domain2architecture'''.format(path_hmmscan = path_hmmscan)
            elif option in ['domain', 'domain2architecture', 'dom']:
                cmd = '''hmmer2table -s {path_hmmscan} | domain2architecture | architecture2table '''.format(path_hmmscan = path_hmmscan)
            else:
                pass

            p = Popen(cmd, shell = True, stdout = PIPE, stderr = PIPE)
            df = pd.read_csv(StringIO(p.communicate()[0].decode('utf-8')), sep = '\t')
            return df

        except KeyboardInterrupt:
            raise error.with_traceback(sys.exc_info()[2])

        finally:
            try:
                os.remove(path_hmmscan)
            except: pass

    def _multi_parse_hmmscan(self,sub_hmmscan, option, ls, insert2db = False,
                             user = 'rotifer', table = 'hmmscan'):
        '''
        Multithread processing
        ----------
        PARAMETERS:
        sub_hmmscan: a fraction of hmmscan result
        option:      output option
        ls:          return result to ls
        '''
        try:
            path_hmmscan = self._tmp()
            with open(path_hmmscan, 'a') as f:
                f.write(''.join(rcf._flatten(sub_hmmscan)))

            if option in ['raw', 'table']:
                cmd = '''hmmer2table -s -y {path_hmmscan}'''.format(path_hmmscan = path_hmmscan)

            elif option in ['architecture', 'architecture2table', 'arch']:
                cmd = '''hmmer2table -s {path_hmmscan} | domain2architecture'''.format(path_hmmscan = path_hmmscan)

            elif option in ['domain', 'domain2architecture', 'dom']:
                cmd = '''hmmer2table -s {path_hmmscan} | domain2architecture | architecture2table '''.format(path_hmmscan = path_hmmscan)

            else:
                pass

            p = Popen(cmd, shell = True, stdout = PIPE, stderr = PIPE)
            out = pd.read_csv(StringIO(p.communicate()[0].decode('utf-8')), sep = '\t')
            # add datetime column here

            if insert2db:
                # Kill nice

                out_dict = out.apply(pd.Series.to_dict, axis=1)
                client = Client('localhost')
                client.execute(f'Insert into {user}.{table} values',
                               list(out_dict.values))

            ls.append(out)

        except KeyboardInterrupt:
            raise error.with_traceback(sys.exc_info()[2])

        finally:
            try:
                os.remove(path_hmmscan)

            except: pass

    def _vars(self):
        '''
        Get vars keys
        '''
        return vars(self).keys()

    def _dbfinder(self, db):
        '''
        '''
        if ':' not in db:
            db = db

        else:
            db = rcf.loadConfig(db)

        return db

    def _tmp(self):
        fd, path = mkstemp()
        return path

class domain:
    '''
    Work with domains
    '''
    def __init__(self, domain):
        '''
        Load a pandas dataframe
        ----------
        PARAMETERS:
        domain: a pandas dataframe
        ATTRIBUTES:
        domain: a pandas dataframe
        '''

        self.domain = domain.copy()

    def domain_len(self, start_col = 'start', end_col = 'end', col_name = '' , inplace = False):
        '''
        Calculate domain length
        ----------
        PARAMETERS:
        start_col: column with start position
        end_col:   column with end position
        col_name:  rename column
        inplace:   return to domain dataframe
        '''

        if not col_name:
            col_name = len(self.domain.columns)

        if inplace:
            self.domain[col_name] = self.domain[end_col].astype(int) - self.domain[start_col].astype(int) +1

        else:
            cop = self.domain.copy()
            cop[col_name] = cop[end_col].astype(int) - cop[start_col].astype(int) +1
            return self._copy(cop)


    def distribution(self, col = 'domain', col_name = 'count', merge = False, sort = True,
                     ascending = False, frequency = False, frequency_col_name = 'frequency'):
        '''
        Calculate domain distribution
        ----------
        PARAMETERS:
        col:                domain column
        col_name:           return column name
        merge:              merge to original dataframe (boolean: True/False)
        sort:               sort values
        ascending:          ascending values (lower to higher) (boolean: True/False)
        frequency:          Calculate frequency of each domain (boolean: True/False)
        frequency_col_name: Frequency column name (default: frequency)
        '''
        if merge:
            g = b.domain.groupcol(col).size().to_frame()
            g.columns = [col_name]
            g = g.reset_index()

            if frequency:
                total = g[col_name].sum()
                g[frequency_col_name] = g[col_name]/total

            self.domain = self.domain.merge(g, on = col)

        else:
            g = b.domain.groupcol(col).size().to_frame()
            g.columns = [col_name]
            self.domain_distribution = g.reset_index()

            if frequency:
                total = g[col_name].sum()
                self.domain_distribution[frequency_col_name] = self.domain_distribution[col_name]/total

            if sort:
                self.domain_distribution.sort_values(col = col, ascending = False)

            return self.domain_distribution


    def percentiles(self,col, percentiles = [10,90]):
        '''
        Find percentile of a size distribution
        ----------
        PARAMETERS:
        col:         column with size
        percentiles: a list of percentiles (default: [10, 90])
        ----------
        RETURN:
        a dictionary
        '''

        dc = {}
        for arg in percentiles:
            dc[arg] = np.percentile(self.domain[col], arg)

        return dc

    def add_sequence(self, *sequence_fasta, on ='ID', inplace = False, how = 'left'):
        '''
        Add fasta sequence to domain dataframe.
        Useful to slice fasta sequence.
        ----------
        PARAMETERS:
        sequence_fasta: input fasta file/sequence
        on:             merge sequence on domain (default = 'ID')
        inplace:        Inplace (boolean: True/False)
        how:            How to merge (default: 'left')
        '''

        seq_df = tb.fasta2df(rcf._fasta_ls(rcf._flatten(sequence_fasta)))
        for col in seq_df.columns:
            if col in self.domain.columns and col != on:
                x = 1

                while True:
                    if col+str(x) in self.domain.columns:
                        x +=1

                    else:
                        seq_df.rename(columns={col: col+(str(x))}, inplace=True)
                        break
        if inplace:
            self.domain = self.domain.merge(seq_df, on = on, how = how)

        else:
            cop = self.domain.copy()
            cop = cop.merge(seq_df, on = on, how = how)
            return self._copy(cop)


    def seq_len(self,col = 'Seq', col_name = '', inplace = False):
        '''
        Calculate sequence (str) length.
        ----------
        PARAMETERS:
        col:      Column with sequence string
        col_name: rename output column
        inplace:  Inplace (boolean: True/False)
        '''

        if not col_name:
            col_name = len(self.domain.columns)

        if inplace:
            self.domain[col_name] = self.domain[col_name]

    def _copy(self, df):
        return domain(df)

    def __repr__(self):
        pd.set_option('precision', 2)
        s = self.domain.to_string()
        return s

# fa1 = '''>AIF01865.1 Ammonium transporter (amt, AMT, MEP) [uncultured marine thaumarchaeote KM3_150_B03] >AIF13962.1 ammonium transporter (amt, AMT, MEP) [uncultured marine thaumarchaeote KM3_65_D04]
# MNSKNHKYALLLMAAVALTATGAMAQAYAQSTVDGMDGYQVGVPGAGIYTGNPNECWYDTDDDGTPDMYCFIDTGDTAWM
# LTASALVLFMTPGVAFLYGGLARSKNAVNTIGMTFIVIGLISVQWVLWGYSLAFGSVDNEANMFMGNLDYVGFNQVSHWA
# PLGEPGSCEGTWSDYYQMQQMKTTGYCSQGWPGTVPHQLFAMFQATFAIITPALIVGGLVDRMKFSALVIFILLWGTFVY
# DPIAHWVWGGGYIGNLDFDPDLSPSYGLDFAGGTVVHITSGFAALAAALVLGRRLGYGKVPMEPHNIPMVVLGASILWFG
# WFGFNAGSEVLADGITVSAWTVTNTATGMASVTWLLMSWGHTGRPSIAGAATGAVAGLVAITPASGWVGPMASIIIGVAA
# GTLCYGAVAFKNSRKWDDALDVWGVHGIGGFTGAVLTGTLASPHIWDTGDGIGAWTGTAEGFEQQAINIVAACMSVAYTF
# AVTIAILKIMDAIWPGGIRVTPREEEVGLDIAQNGERAYVFE'''.split('\n')
#
# fa2 = '''>AIE94919.1 ammonium transporter (amt, AMT, MEP) [uncultured marine thaumarchaeote AD1000_54_F06]
# MKSRNHKYALILMAAVAMTATGALSTAYAQQTSDGMDGYVVGDPDGGAGIYTGNPNECWYDTDDDGVPDMYCYVDTGDMA
# WMLTASSLVLFMTPGVAFFYGGLARSKNAVNTIGMVFIIMGLMSVQWVLWGYSLAFGGIDNDANMFMGNLDYVGFNQVSH
# WAPLGAPSPCEDTWAHYYQMQTMKEGDVCSDTWPGTVPHQLFAMFQATFAIITPALIVGGLVDRMKFSAIVVFVLLWGTF
# VYDPICHWVWGGGYIGGGSLDFNPDLSPSYALDFAGGTVVHISSGFTALAGALILGRRLGYGKVPMEPHNVPMVVLGASI
# LWFGWFGFNAGSEVMVDGITVSAWTVTNTATGMATVTWVLMSWAHTGKPSIVGAATGAVAGLVAITPASGWVGPMAAIII
# GIAAGTVCYGCVAFKNARKWDDALDVWGVHGMGGLTGAILTGTLASPHIWDTGDGIGAWTGTPEGYEQQAISIVAAGMSI
# AYAFGISIVIFKVMDAVWPGGIRVTPKEEEIGLDLSQNGERAYVNE'''
#
# fa3 = '''>WP_011523244.1 sensor signal transduction histidine kinase [Candidatus Koribacter versatilis]
# MISTGRVRIQPSLVVWFLLLISAGSSLFALNPDLSISQYAHSTWRVQDGAFRSAPNAVAQTKDGYLWIGTEGGLVHFDGV
# RFVPWVPPAGVKLLDPRIFSLMAASDGSLWIGTGYSISHWRRNELINYSQLSGRIEAIAEDHDGTVWFVRTQITDGGGPV
# CRITNDQPQCFGKADGIPFPIAVQLRVGNSGELWVGGYSELCRWKPASLSSDCFAKGSQVPETFASIKAIATGNDGTVWV
# ARERPGSFLQLERFAQEKWTTLSYPEIAINNSDVTTLFVDRDNTIWVGSANHGVFRIVGNTVRSFGRTDGLSSDAVGRFY
# QDVEGTVWVVTSAGIDNFRDLKVVTYSMREGLTAAGAGTVLGTRDGTVWIGNFHALDFMQGSKLSSIRAGNGLPGLYITT
# FFEDHAGRLWVGIDDGLWVYENQTFRPVRHADGSKLGIVFSITEDTLHNIWARAGKNLDRIADYRLQEETTSPQISTSYI
# LAASPQGGIYLGLVSGDLVQYDGGKSQTFASNEVGNTRQIRDLLVEPDGSVWGTTLDEIVRWKNGERKNLTTRNGLPCDE
# IFALVEDSRGSLWIESKCGVIEIERAQLDAWWEHPETVVKFGLLDGSDGMQAGLTPLKPQATRSADGKLWFVNGRILQML
# DPNHLQRNPVPPPVQIEEIVADHKSHSPQAGLRLPALTRDLEIDYTALSFVAPQKVQFRYMLEGRDTAWQETMTRRQAFY
# NNLGPGHYRFRVMASNNDGVWNEAGAYLDFSILPAYYQTVWFRLLCAIAFLMVLWSIFQIRVHQLRRQFEIGVEARVNER
# TRIARELHDTLLQTLHGLMFQFQAVRNLLPRRPEDAMRSLDDAIVETEKALAEGRNAIQGIRSESRDGDDLAEFLKNASK
# DFASTAKPGESLPTFDLIEEGQRRSVSSDVNNEVCRIALELLRNAFRHAQATRIEAEIRYDAQMLRLRIRDNGKGIDPVV
# LREGGVAGHWGLKGVRERAERIGAKIEFWSDVGLGTEIQVTVPAGVAYQAESEERLLQSNPRTKSRAKQS'''

# a = search(fa1, fa3, fa2,  methods = {'hmmscan': {}})
# a.parser(hmmscan = {'of':'domain'})
# b = domain(a.hmmscan_parsed)
# #
# b.add_sequence([fa1,fa2,fa3], inplace = True)
# #
