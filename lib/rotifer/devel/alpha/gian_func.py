from networkx.algorithms.bipartite.basic import color
from networkx.generators import small


def load_seq_scan(name, folder, haldane=False):
    '''
    load a seqscan result into a dataframe
    '''
    import pandas as pd
    pre = f'{folder}/{name}'
    info = pd.read_csv(
            f'{pre}.c100i100.tsv', sep='\t', names=['c100i100', 'pid']
            )
    info = info.merge(
            pd.read_csv(
                f'{pre}.c80i70.tsv',
                sep='\t',
                names=['c80i70', 'c100i100']
                ), how="left"
            )
    info = info.merge(
            pd.read_csv(
                f'{pre}.c80e3.tsv',
                sep='\t',
                names=['c80e3', 'c80i70']
                ),
            how="left")
    if haldane:
        info = info.merge(
                pd.read_csv(
                    f'{pre}.aravind.scan.arch',
                    sep='\t',
                    names=['c100i100', 'aravind'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
        info = info.merge(
                pd.read_csv(
                    f'{pre}.pfam.hmmscan.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
    else:
        info = info.merge(
                pd.read_csv(
                    f'{pre}.query.profiledb.rps.arch',
                    sep='\t',
                    names=['c100i100', 'aravind'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
        info = info.merge(
                pd.read_csv(
                    f'{pre}.query.pfam.rps.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]),
                how="left")
    return info


def cluster2aln(
        group_cluster,
        df,
        esl_index_file,
        grouper='c80e3',
        redundancy_cluster='c80i70',
        align_method='famsa',
        query=False,
        cpu=12
        ):
    import os
    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence
    with tempfile.TemporaryDirectory() as tmpdirname:
        if query:
            df.query(
                    query
                    )[redundancy_cluster].drop_duplicates().dropna().to_csv(
                            f'{tmpdirname}/accs',
                            index=None,
                            header=None
                            )
        else:
            df[df[grouper] == group_cluster][
                    redundancy_cluster
                    ].drop_duplicates(
                            ).dropna(
                                    ).to_csv(
                                            f'{tmpdirname}/accs',
                                            index=None,
                                            header=None
                                            )
        if not os.path.exists(esl_index_file+'.ssi'):
            Popen(
                    f'esl-sfetch --index {esl_index_file}',
                    stdout=PIPE,
                    shell=True
                    ).communicate()
        Popen(
                f'esl-sfetch -f {esl_index_file} {tmpdirname}/accs > \
                        {tmpdirname}/accs.fa',
                stdout=PIPE,
                shell=True
                ).communicate()
        b = sequence(
                f'{tmpdirname}/accs.fa'
                ).align(
                        method=align_method,
                        cpu=cpu)
        return b


def chunks(list_to_chunk, chunk_size):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(list_to_chunk), chunk_size):
        yield list_to_chunk[i:i + chunk_size]


def cluster_Co_occurrence(
        df,
        count='c80e3',
        freq_cutoff=0.3,
        only_query=True,
        annotation=False):
    ''' Function to count the co-occurence of clusters within
    the NeighborhoodDF.
    The count parameter shoul use to define the cluster
    that would be analysed and
    the freq_cutoff is to use to define a minimum cut_off
    to display
    '''
    x = df[['block_id', count]].merge(
            df[['block_id', count]],
            how='outer',
            on='block_id'
            )
    x.columns = [
            'block_id',
            'query_cluster',
            'neighbor_cluster'
            ]
    x = x[x.query_cluster != x.neighbor_cluster]
    if only_query:
        xx = x[
                x['query_cluster'].isin(
                    df.query('query ==1').pid.unique()
                    )
                ]
    else:
        xx = x

    xxx = xx.groupby(
            ['query_cluster', 'neighbor_cluster']
            ).block_id.nunique().reset_index()
    xxxx = xx.groupby(
            'query_cluster'
            ).block_id.nunique().rename(
                    'query_blocks'
                    ).reset_index().merge(
                            xxx,
                            how='left'
                            ).sort_values(
                                    ['query_blocks', 'block_id'],
                                    ascending=False
                                    )
    xxxx['query_freq'] = xxxx.block_id/xxxx.query_blocks
    if not annotation:
        return xxxx.query('query_freq >= @freq_cutoff')

    andf = df.groupby(count).agg(
            pfam=('pfam', count_series),
            aravind=('aravind', count_series)
            ).reset_index()
    xxxx = xxxx.merge(
            andf.rename({
                count: 'query_cluster',
                'pfam': 'query_pfam',
                'aravind': 'query_aravind'
                }, axis=1),
            how='left')
    xxxx = xxxx.merge(
            andf.rename({
                count: 'neighbor_cluster',
                'pfam': 'neighbor_pfam',
                'aravind': 'neighbor_aravind'
                }, axis=1),
            how='left')
    return xxxx.query('query_freq >= @freq_cutoff')


def count_series(
        series,
        normalize=False,
        cut_off=False,
        count='domain',
        as_list=False
        ):
    '''
    Function to flatten a pd.Series.value_counts
    and count the number of domains recognized in
    the whole proteins set of proteins
    Mind that by default, it count domain not architecture.
    If one wants to count architecture, it should change the count
    option to 'architecture'.
    The normalize options is to write the results as frequency.
    The cut_off options is to print results only above a given threashold.
    If the normalize function is True, you should give the cutoff value as pct.
    Usage:
    in: count_series(df.pfam, normalize=True, cut_off=10)
    out:'T6SS_HCP(15.15%), SP(11.36%)'
    in: df.groupby('c80e3').agg({'pfam':count_series})
                                                 pfam
    c80e3
    AAN66278.1
    ABO22067.1                           SP(1), TM(1)
    ABU79335.1  Gal_mutarotas_2(1), Glyco_hydro_31(1)
    ABV15410.1                           SP(1), TM(1)
    ABV15414.1                                  TM(1)

    '''
    flattened_list = []
    if count == 'architecture':
        s = series.dropna().value_counts(normalize=normalize)
    else:
        s = series.str.split('+').explode().dropna().value_counts(
                normalize=normalize
                )

    if cut_off:
        #cut_off = cut_off/100
        s = s.where(lambda x: x >= cut_off).dropna()
    for y, z in s.items():
        if normalize:
            flattened_list.append(f'{y}({100 * z:.2f}%)')
        else:
            flattened_list.append(f'{y}({z})')
    
    if as_list:
        if len(flattened_list) ==0:
            return list([''])
        else:
            return flattened_list

    return ', '.join(flattened_list)


def annotation(seqobj, coordinates, delimiter=True):
    '''
    Method that receives an sequence object and a list with tuples
    containing [(start, end annotation)]
    To add annotaion in aligment
    Example:
        gf.annotation(
        seqobject,
        [(35, 217,
        ' domain 1 (probability 90%, WP_xxxxxx)'),
        (240, 380, 'domain 2')
        ] )
    I sill have to finish it, but I  am almost done
    '''
    import pandas as pd
    s = seqobj.copy()
    t = pd.Series(list(s.df.iloc[0, 1]))

    t.iloc[:] = ' '
    for x in coordinates:
        if len(x) > 3:
            query_anchor = x[3]
            tq = pd.Series(
                    list(
                        s.df.query(
                            'id ==@query_anchor'
                            ).sequence.values[0]
                        )
                    )
            start = tq.where(
                    lambda x: x != '-'
                    ).dropna().iloc[x[0]:x[1]].index.min()
            end = tq.where(
                    lambda x: x != '-'
                    ).dropna().iloc[x[0]:x[1]].index.max()
        else:
            start = x[0]
            end = x[1]

        annotation = x[2]
        size = end - start
        size_an = len(annotation)
        if (size - 2) < size_an:
            x = size - size_an + - 2
            annotation = annotation[0:x]

        if delimiter:
            an = f'|{annotation}|'.center(size, '*')
        else:
            an = f'{annotation}'.center(size, '*')
        t.update(pd.Series(list(an), index=range(start, end)))
    s.df = pd.concat(
            [pd.DataFrame(
                [['annotation', ''.join(t.to_list()), 'annotation']],
                columns=['id', 'sequence', 'type']
                ),
             s.df
             ])
    return s


def hmmsearch_full2pandas(
        file,
        error_lines=True,
        keep_threshold=False
       ):
    """
    Function to load the result of hmmsearch or scan ang obtain
    the results for the full protein.
    If error_lines = False, it skip the results that gave error
    loading to the dataframe.
    Those errors only happens if any field present an string with
    more than two spaces in sequence (mostly found in  descrption field )
    """
    import pandas as pd
    import re
    from io import StringIO
    with open(file, 'r') as f:
        text = f.read()
    columns = ['Evalue',
               'score',
               'bias',
               'BD_Evalue',
               'BD_score',
               'BD_bias',
               'exp',
               'N',
               'Sequence',
               'Description'
               ]
    # u = '------- ------ -----    ------- ------ -----   ---- --  --------       -----------'
    u = '------- ------ -----    ------- ------ -----   ---- --  --------   -----------'
    string_tore = 'Domain annotation for each sequence'
    match = re.findall(f'{u}(.+?){string_tore}', text, re.DOTALL)[0].strip()
    m2 = re.sub(r'[^\S\r\n]{2,}', '\t', match)
    m3 = m2.replace('\n\t', '\n')
    df = pd.read_csv(
            StringIO(m3),
            sep="\t",
            names=columns,
            usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            error_bad_lines=error_lines)
    if keep_threshold:
        return df

    df = df.query('~Evalue.str.contains("threshold")')
    df.Evalue = df.Evalue.astype(float)

    return df


def hhr_to_aln(seqobj, hhr, database=False):

    """
    Function to add hhr info into sequence object aligment
    """
    import rotifer.interval.utils as ru
    from collections import OrderedDict
    import re
    import pandas as pd

    si = re.findall('No ([0-9]+)', hhr, re.DOTALL)
    si = [int(x) for x in si]
    query = re.findall(r'Query\s+(.*?)\s', hhr, re.MULTILINE)[0]
    tdf = pd.DataFrame()
    for x in range(len(si)):
        try:
            match = re.findall(
                f'No {si[x]}(.+?)No ',
                hhr,
                re.DOTALL
            )[0].strip()
        except IndexError:
            match = re.findall(
                f'No {si[x]}(.+)\n',
                hhr,
                re.DOTALL
            )[0].strip()
        if database == "pfam":
            try:
                idp = re.findall(
                    r'>.+?;(.+?);',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    r'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()
        else:
            try:
                idp = re.findall(
                    r'>(.+?)\.',
                    match,
                    re.DOTALL
                )[0].strip()
            except IndexError:
                idp = re.findall(
                    r'>(.+?)\n',
                    match,
                    re.DOTALL
                )[0].strip()

        values = match.split('\n')[1].split()
        values = pd.Series(values).str.split(
            '=',
            expand=True
        ).rename({
            0: 'col',
            1: 'val'},
                 axis=1)
        values['model'] = idp
        values = values.pivot(
            index='model',
            columns='col',
            values='val'
        )
        values['query'] = query
        values['start'] = re.findall(
            rf'{query}\s+([0-9]+) ',
            match,
            re.DOTALL
        )[0]
        values['end'] = re.findall(
            rf'{query}.+?([0-9]+?) \(',
            match,
            re.DOTALL)[-1]

        tdf = pd.concat([tdf, values])

    hhr_nolr = ru.filter_nonoverlapping_regions(
        tdf.astype({'start': int, 'end': int}).reset_index(),
        end='end',
        start='start',
        reference='query',
        criteria=OrderedDict([('Probab', False), ('region_length', True)]
                             ))
    # ISSUE !!!!!! When the rep seq do not have the first domains
    # it breaks the code because there is no way to map the domain.
    # This is an imporvisation that will overpass it,
    # but with error in annotation!!
    hhr_nolr = hhr_nolr[hhr_nolr.end > hhr_nolr.start]
    return hhr_nolr


def add_arch_to_df(seqobj, full=False, inplace=True):
    '''
    Add architecture and clusters info to a neighborhood df
    '''

    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    seqobjc = seqobj.copy()
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdirname:
        if hasattr(seqobjc, 'ndf'):
            df = seqobjc.ndf
        # temporary save fasta sequence file
            sequence(
                    df.query(
                        'type =="CDS"'
                        ).pid.dropna().drop_duplicates().to_list()
                    ).to_file(f'{tmpdirname}/seqfile')
        else:
            print('Seq object do not have neighborhood df, running annotation on aligment sequences')
            df = seqobjc.df
            sequence(
                    df.id.dropna().drop_duplicates().to_list()
                    ).to_file(f'{tmpdirname}/seqfile')

        os.chdir(tmpdirname)
        os.makedirs('tmpd')
        if full:
            script = 'scanseqs.nih.sh'
        else:
            script = 'scanseqs.short.nih.sh' 

        Popen(
                f'/home/nicastrogg/bin/{script} tmp seqfile',
                stdout=PIPE,
                shell=True
                ).communicate()
        info = pd.read_csv(
                'tmp.c100i100.tsv',
                sep='\t',
                names=['c100i100', 'pid']
                )
        info = info.merge(
                pd.read_csv(
                    'tmp.c80i70.tsv',
                    sep='\t',
                    names=['c80i70', 'c100i100']
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.c80e3.tsv',
                    sep='\t',
                    names=['c80e3', 'c80i70']),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.profiledb.rps.arch',
                    sep='\t',
                    names=['c100i100', 'profiledb'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.pfam.rps.arch',
                    sep='\t',
                    names=['c100i100', 'pfam'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        info = info.merge(
                pd.read_csv(
                    'tmp.query.both.rps.arch',
                    sep='\t',
                    names=['c100i100', 'arch'],
                    usecols=[0, 1],
                    skiprows=[0]
                    ),
                how="left")
        if inplace:
            if hasattr(seqobjc, 'ndf'):
                seqobj.ndf = df.merge(info, how='left')
                seqobj.ndf = seqobj.ndf.drop_duplicates().reset_index(drop=True)
            else:
                info = info.rename({'pid':'id'}, axis=1)
                seqobj.df = df.merge(info, how='left')

            seqobj.profiledb = pd.read_csv(
                    'tmp.query.profiledb.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'star', 'end'])
            with open('tmp.profiledb.rps.query.out') as f:
                seqobj.profiledbout = ''.join(f.readlines())
            seqobj.pfam = pd.read_csv(
                    'tmp.query.pfam.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'start', 'end'])
            with open('tmp.pfam.rps.query.out') as f:
                seqobj.pfamout = ''.join(f.readlines())
            seqobj.arch = pd.read_csv(
                    'tmp.query.both.rps.dom.tsv',
                    sep='\t', names=['pid', 'model', 'start', 'end'])
            os.chdir(cwd)
            return 'architecuture add to seqobj'

        ### The annotation of arch on sequence objct is not working if not inplace ########
        with open('tmp.pfam.rps.query.out') as f:
            seqobjc.pfamout = ''.join(f.readlines())

        seqobjc.pfam = pd.read_csv(
                'tmp.query.pfam.rps.dom.tsv',
                sep='\t')

        seqobjc.arch = pd.read_csv(
                'tmp.query.both.rps.dom.tsv',
                sep='\t')

        seqobjc.profiledb = pd.read_csv(
                'tmp.query.profiledb.rps.dom.tsv',
                sep='\t')
        with open('tmp.profiledb.rps.query.out') as f:
            seqobjc.profiledbout = ''.join(f.readlines())
        seqobjc.ndf = df.merge(info, how='left')
        seqobjc.ndf = seqobj.ndf.drop_duplicates().reset_index(drop=True)
        os.chdir(cwd)
        
    return seqobjc


def annotate_seqobj(seqobj,
                    ndf = False,
                    cnt='arch',
                    selected_collumns=False
                    ):
    '''
    Add annotation to seqobject tax,clusters,compact_neighborhood
    '''

    import pandas as pd

    if isinstance(ndf, pd.DataFrame):
        df = ndf.copy()
    else:
        df = seqobj.ndf.copy()
    seqobj = seqobj.copy()
    accs = seqobj.df.id.to_list()
    cn = df[
            df.block_id.isin(df.query('pid in @accs').block_id)
            ].compact_neighborhood(cnt)
    if selected_collumns:
        selected_collumns = selected_collumns
    else:
        selected_collumns = ['pid',
         'block_id',
         'profiledb',
         'pfam',
         'arch',                    
         'c80e3',
         'c80i70']

    seqobj.df = seqobj.df.merge(df.query('pid in @accs')[selected_collumns ].rename(
            {'pid': 'id'}, axis=1
            ).set_index('block_id').join(cn), on='id', how='left')

    return seqobj


def psiblast(acc,
             db='latest50',
             cpu=96,
             aln=True,
             max_out = 5000):
    '''
    Psiblast it can accept sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =[
            'hit',
            'query',
            'hitstart',
            'hitend',
            'evalue',
            'querycoverage',
            'querystart',
            'queryend',
            'iteration',
            'bitscore',
            'length'
            ]
    print(db)
    if db =='latest50':
        db = f'{os.environ["FADB"]}/nr/latest/nr.50.mmseqs.fa'
    elif db.startswith('nr'):
        db = f'{os.environ["FADB"]}/nr/{db}.mmseqs.fa'
    elif db == 'all.fa' :
        db = f'{os.environ["FADB"]}/all.fa'
    elif db == "prok" :
        db = f'{os.environ["FADB"]}/prok.fa '
    elif db == "euk":
        db = f'{os.environ["FADB"]}/euk.fa'
    elif db == "af":
        db = f'{os.environ["FADB"]}/uniref50.fa'
    else:
        db = db

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            acc.to_file(f'{tmpdirname}/seqfile') 
            if aln:
                Popen(f'splishpsi -a {cpu} -in_msa {tmpdirname}/seqfile -d {db} -num_alignments {max_out} -num_descriptions {max_out} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            else:
                Popen(f'splishpsi -a {cpu} -i {tmpdirname}/seqfile -d {db} -num_alignments {max_out} -num_descriptions {max_out} > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            Popen(f'/netmnt/vast01/cbb01/proteinworld/tools/bin/blast2table {tmpdirname}/out > {tmpdirname}/out.tsv',
                  stdout=PIPE,
                  shell=True).communicate()
            t = pd.read_csv(f'{tmpdirname}/out.tsv', sep="\t", names=cols)
            with open(f'{tmpdirname}/out') as f:
                    blast_r = f.read()
        
    return (t, blast_r) 


def search2aln(df, coverage=50, evalue=1e-3, id_with_coord = False, arch=None):
    import os
    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import pandas as pd
    if coverage > 1:
        coverage = coverage/100
    '''
    From a search tabular file, filter by cov, and evalue to make a aln object
    '''
    df = df.query(f'evalue <= {evalue} and querycoverage > {coverage}')
    seqobj = sequence(df.hit.drop_duplicates().to_list())
    if arch == 'profiledb':
        arch = ' '
    if arch:
        with tempfile.TemporaryDirectory() as tmpdirname:
            seqobj.to_file(f'{tmpdirname}/fasta')
            Popen(f'cat {tmpdirname}/fasta|splishrps {arch}|rps2arch -pcut 0.001 -overlap 0.03 > {tmpdirname}/arch', stdout=PIPE,shell=True).communicate()
            archdf = pd.read_csv(f'{tmpdirname}/arch',sep="\t", names=['id', 'arch', 'position'])

    seqobj.df = seqobj.df.merge(df[['hit','evalue','hitstart','hitend']].rename({'hit':'id'},axis=1), how='left')
    seqobj.df.sequence = seqobj.df.apply(lambda x: x.sequence[x.hitstart:x.hitend], axis=1)
    seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)
    evalues = seqobj.df.copy()[['id', 'evalue']]
    seqobj = seqobj.align()
    seqobj.df = seqobj.df.merge(evalues,on=['id'], how='left')
    seqobj = seqobj.sort(['evalue'])
    if arch:
        # We should fix the arch with coordinates similar to the evalue strategy
        print('arch option may contains bug!i')
        seqobj.df = seqobj.df.merge(archdf, how='left')

    #seqobj.df.id = seqobj.df.apply(lambda x: f'{x.id}/{x.hitstart}-{x.hitend}', axis=1)

    if id_with_coord:
        seqobj.df = seqobj.df.drop_duplicates(['id']).reset_index(drop=True)
        return seqobj
    seqobj.df[['hit_start', 'hit_end']] = seqobj.df.id.str.split('/', expand=True)[1].str.split('-', expand=True)
    seqobj.df.id = seqobj.df.id.str.split('/', expand=True)[0]
    seqobj.df = seqobj.df.drop_duplicates(['id', 'hit_start', 'hit_end']).reset_index(drop=True)
    return seqobj


def add_arch_to_seqobj(seqobj,db='profiledb', cpu=96):
    '''
    Psiblast it can accept sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =['id','arch','evalue']
    if db == 'profiledb':
        db = ' '

    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        acc = sequence(seqobj.df.id.to_list())
        acc.to_file(f'{tmpdirname}/seqfile') 
        Popen(f'cat {tmpdirname}/seqfile| splishrps -a {cpu} {db} > {tmpdirname}/out' , stdout=PIPE,shell=True).communicate()
        Popen(f'rps2arch {tmpdirname}/out > {tmpdirname}/out.tsv' , stdout=PIPE,shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/out.tsv', sep='\t', names=cols)
        seqobj.df = seqobj.df.merge(t, how='left', on='id') 
    return seqobj


def full_annotate(seqobj,
                  progress=True,
                  batch_size=8,
                  mirror="/am/ftp-genomes",
                  threads=None,
                  after=5,
                  before=5,
                  eukaryotes=False,
                  full=True):
    from rotifer.db import ncbi 
    from rotifer.devel.alpha import gian_func as gf
    gnc = ncbi.GeneNeighborhoodCursor(
            progress=progress,
            batch_size=batch_size,
            mirror=mirror,
            threads=threads,
            after=after,
            before=before,
            eukaryotes=eukaryotes)
    seqobj.ndf = gnc.fetchall(seqobj.df.id.to_list())
    gf.add_arch_to_df(seqobj, full=full)
    return seqobj

def padding_df(df, how='right'):
    cdf = df.copy()
    c = df.columns
    pad_col_name=[]
    for x in c:
        cdf[x] = cdf[x].fillna('').astype(str)
        w = cdf[x].str.len().max()
        cdf[x] = cdf[x].str.pad(width =w, side=how)
        if how == 'left':
            pad_col_name.append(x.rjust(w))
        else:
            pad_col_name.append(x.ljust(w))

    cdf.columns = pad_col_name    
    return cdf


def alnxaln(seqobj, clustercol = 'c50i0', minseq=10, aln_method='famsa'):
    """ IF aln_method other than FAMS, mafft linsi will be used"""
    import os
    import tempfile
    import pandas as pd
    from rotifer.io import hhsuite
    
    path = os.getcwd()
    result_table = pd.DataFrame()
    if clustercol in seqobj.df.columns:
        pass
    else:
       print(f'add the cluster{clustercol} columns to your seqobj')
       return

    with tempfile.TemporaryDirectory() as tmpdirname:
        os.chdir(tmpdirname)
        seqobj = seqobj.copy()
        alndict = dict()
        clist = list(seqobj.df[clustercol].value_counts().where(lambda x: x>= minseq).dropna().index)
        for x in clist:
            if aln_method == 'famsa':
                alndict[x] = seqobj.filter(f'{clustercol} ==\"{x}\"').align()
            else:    
                alndict[x] = seqobj.filter(f'{clustercol} ==\"{x}\"').align(method='linsi')
            alndict[x].to_file(f'{x}.aln') 
        for x in clist:
            os.system(f'for x in *.aln; do hhalign -i $x -t {x}.aln;done')
            os.system('cat *.hhr >> allhhr.txt')
            hhrs = hhsuite.read_hhr('./')
            result_table = pd.concat([result_table, hhrs])
        with open('./allhhr.txt', 'r') as f : allhhr = f.read()        
    os.chdir(path)
    return (result_table, allhhr, alndict) 


def split_by_model(seqobj, model, evalue=1e-6):
    '''
    Using a hmm model to split your sequence object to match only the model region.
    If more than one match in one protein, it will split the match in two sequence.
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    import rotifer
    from rotifer.core import loadpath
    
    hmmer2table_from_rotifer = f'{loadpath.path().globalpath}rotifer/bin/hmmer2table'

    nseqobj =seqobj.copy()
    with tempfile.TemporaryDirectory() as tmpdirname:
        nseqobj.df.sequence = nseqobj.df.sequence.str.replace('-','')
        nseqobj.to_file(f'{tmpdirname}/seqfile') 

        if isinstance(model, rotifer.devel.beta.sequence.sequence):
            model.to_file(f'{tmpdirname}/model')
            Popen(f'hmmbuild {tmpdirname}/model.hmm {tmpdirname}/model',
                  stdout=PIPE,
                  shell=True).communicate()
            model = f'{tmpdirname}/model.hmm'

        Popen(f'hmmsearch -E {evalue} {model} {tmpdirname}/seqfile > {tmpdirname}/hhsearch_r',
              stdout=PIPE,
              shell=True).communicate()
        Popen(f'{hmmer2table_from_rotifer} {tmpdirname}/hhsearch_r |domain2architecture |architecture2table > {tmpdirname}/hhsearch_table',
              stdout=PIPE,
              shell=True).communicate()
        t = pd.read_csv(f'{tmpdirname}/hhsearch_table', sep="\t")

    nseqobj.df = nseqobj.df.merge(t.rename({'ID':'id'},axis=1), how='left')
    nseqobj.df.end = nseqobj.df.end.fillna(nseqobj.df.length)
    nseqobj.df.start = nseqobj.df.start.fillna(1)
    nseqobj.df.sequence = nseqobj.df.apply(lambda x: x.sequence[int(x.start):int(x.end)], axis=1)
    nseqobj = nseqobj.filter('domain =="model"').align()
    return (nseqobj)

def remove_redundancy(seqobj, identity =80, coverage = 70):
    from rotifer.devel.beta.sequence import sequence as sequence
    import pandas as pd
    s = seqobj.copy()
    if len(s.df.id.value_counts().value_counts()) >1:
        print('multiples sequences with same name, using index to remove redundancy')

    s.add_cluster(coverage=coverage, identity=identity, inplace=True)
    s = s.filter(keep=s.df[f'c{coverage}i{identity}'].drop_duplicates().to_list())
    return(s)

def add_cordinates_to_aln(seqobj, full_length_seqobj=False):
    from rotifer.devel.beta.sequence import sequence as sequence
    c = seqobj.copy()
    if full_length_seqobj:
        cx = full_length_seqobj.copy()
    else:    
        cx = sequence(c.df.id.tolist())
    cx.df = cx.df.rename({'sequence':'full_sequence', 'length':'full_length'}, axis=1)
    c.df = c.df.merge(cx.df[['id','full_sequence','full_length']], how='left')
    c.df['start'] = c.df.apply(lambda x : 0 + x.full_sequence.find(x.sequence.replace('-','')), axis=1)
    c.df['end'] = c.df.start + c.df['length']
    c.df['C_term'] = c.df['full_length'] - c.df['end']
    c.df =  c.df.drop(['full_sequence'], axis=1)
    return c


def extend_aln(seqobj, n_terminal=50, c_terminal=50,c_term_list=[],n_term_list=[], full_length_seqobj=False):
    from rotifer.devel.beta import sequence as rdbs
    from rotifer.devel.alpha import gian_func as rdagf
    import numpy as np
    """
    Extend the N and/or C termianl parts of a given sequence alignment
    """
    if len(c_term_list) == 0:
        c_term_list = seqobj.df.id.tolist()
    if len(n_term_list) == 0:
        n_term_list = seqobj.df.id.tolist()
    tm = rdagf.add_cordinates_to_aln(seqobj, full_length_seqobj=full_length_seqobj)
    tm.df['previous_length'] = tm.df['length']
    if full_length_seqobj:
        tm2 = full_length_seqobj.copy()
    else:    
        tm2 = rdbs.sequence(seqobj.df.id.tolist())
    tm2.df = tm2.df.merge(tm.df[['id', 'full_length', 'start', 'end', 'C_term', 'previous_length']], how='left')
    tm2.df.end = np.where(tm2.df.id.isin(c_term_list),
                          np.where(
                              tm2.df.end + c_terminal <= tm2.df.full_length,
                              tm2.df.end + c_terminal,
                              tm2.df.full_length),
                          tm2.df.end)
    tm2.df.start = np.where(tm2.df.id.isin(n_term_list),
                           np.where(
                                tm2.df.start - n_terminal > 0,
                                tm2.df.start - n_terminal,
                                1),
                            tm2.df.start)
    if tm2.df['end'].isna().any():
        te = ', '.join(tm2.df.query("end.isna()").id.tolist())
        print(f'some pid is missing the coordinates  {te}')
        tm2.df = tm2.df.query('~end.isna()')
        tm2.df.start = tm2.df.start.astype(int)
        tm2.df.end = tm2.df.end.astype(int)
        tm2.df = tm2.df.drop_duplicates('id')
    tm2.df.sequence = tm2.df.apply(lambda x : x.sequence[x.start -1: x.end], axis=1)
    return tm2.align()

def get_correspondent_position(seqobj_source,seqobj_target, pid, position):
    """ Get the correspondent amminoacid position from one aligment/single sequence (rotifer seqobj) in a different sequence aligment
    seqobj_source: Sequence object where the user knows the desired position
    seqobt_target: Sequence object where the user need to find the correspondent position
    pid: Sequence used as anchor to find the correspondent position, it should be present in both source and target sequence object
    position: The postion in the pid wherethe user wants to find the correspondet position in the target sequence objct
    """
    from Bio import pairwise2
    import pandas as pd

    source =  seqobj_source.df.query('id == @pid').sequence.iloc[0].replace('-', '')
    target =  seqobj_target.df.query('id == @pid').sequence.iloc[0].replace('-','')
    aln = pairwise2.align.globalxx(source, target)[0]
    source = aln[0]
    target = aln[1]
    d = {'seq': list(source), 'trimmed' :list(target)}
    aln_df = pd.DataFrame(data=d).reset_index().rename({'index': 'pos'}, axis=1).query('trimmed != "-"').reset_index(drop=True)
    return int(aln_df.query('pos ==@position -1').index[0]) +1


def read_predicted_topologies(file):
    from rotifer.devel.beta import sequence as rdbs
    import pandas as pd
    import numpy as np
    def g (s):
        tomap = pd.Series(list(s.sequence))
        o = tomap.where(lambda x: x!="-").dropna().rename('aln').to_frame()
        o['top'] = list(s.topology)
        od = o.top.to_dict()
        tomap.update(od)
        return "".join(tomap.tolist())

    a = np.array(['header', 'sequence', 'topology'])
    df = pd.read_csv(file, names=['data'])
    df['data_type'] = np.resize(a, len(df))
    h = df.query('data_type =="header"').data.str.strip('>').rename("id").reset_index(drop=True)
    s = df.query('data_type =="sequence"').data.rename("sequence").reset_index(drop=True)
    t = df.query('data_type =="topology"').data.rename("topology").reset_index(drop=True)
    df2 = pd.concat([h,s,t], axis=1)
    aln = rdbs.sequence(df2[['id', 'sequence']]).align()
    aln.df['topology'] = df2.topology
    aln.df['aln_top'] = aln.df.apply(g, axis=1)
    aln.con = aln.add_consensus().df.iloc[0:4,0:2].rename({'sequence': 'aln_top'}, axis=1)
    aln.topology = pd.concat([aln.con,aln.df[['id', 'aln_top']]]).reset_index(drop=True)
    return aln

def uniref50_add_info(seqobj):
    seqobj = seqobj.copy()
    from rotifer.devel.beta.sequence import sequence
    from rotifer.db.local import ete3
    from rotifer.devel.alpha import gian_func as gf
    dd = sequence(seqobj.df.id.tolist())
    dd.add_cluster(coverage=0, identity=0, inplace=True)
    seqobj.df = seqobj.df.merge( dd.df[['id', 'length', 'c0i0']].rename({'length': 'plen'},axis=1))
    xx = ete3.TaxonomyCursor()
    seqobj.df[['description', 'cluster_size', 'Tax', 'TaxID', 'RepID']] = seqobj.df.description.str.split("n=|Tax=|TaxID=|RepID=", expand=True)
    seqobj = gf.add_cordinates_to_aln(seqobj)
    tid = xx.fetchall(seqobj.df.TaxID.drop_duplicates().tolist())
    tid.taxid = tid.taxid.astype(int)
    tid.rename({'taxid':'TaxID'}, axis=1, inplace=True)
    seqobj.df.TaxID = seqobj.df.TaxID.astype(int)
    seqobj.df = seqobj.df.merge(tid, how='left')
    return seqobj


def uniref50_to_ncbi(uniref50_strings):
     import pandas as pd
     import sqlite3
    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()
     pids = [uniprotid.split('_')[1] for uniprotid in uniref50_strings]

     placeholders = ', '.join(['?' for _ in pids])
     sql_query = f"""
         SELECT *
         FROM EMBL
         WHERE uniprotid IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=pids)
     uniref50DF['Uniref50'] = 'UniRef50_' + uniref50DF.uniprotid
     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return uniref50DF

def uniref50_to_clusters(uniref50_strings):
     import pandas as pd
     import sqlite3
    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()

     placeholders = ', '.join(['?' for _ in uniref50_strings])
     sql_query = f"""
         SELECT *
         FROM uniref50
         WHERE uniref50 IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=uniref50_strings)
     uniprotid = uniref50DF.uniprotid.unique().tolist() 

     placeholders2 = ', '.join(['?' for _ in uniprotid])
     # Construct the SQL query
     sql_query2 = f"""
     SELECT *
     FROM "EMBL"
     WHERE uniprotid IN ({placeholders2});
     """

     # Execute the SQL query
     finalDF = pd.read_sql_query(sql_query2, conn, params=uniprotid)
     finalDF = finalDF.merge(uniref50DF, how='outer').drop_duplicates()


     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return finalDF

def pid2uniref50(pids):
     import pandas as pd
     import sqlite3
     from rotifer.db import ncbi

     i = ncbi.IPGCursor()
     i = i.fetchall(pids)
     i = pd.concat(i)

    # Connect to the SQLite database
     conn = sqlite3.connect('/netmnt/vast01/cbb01/proteinworld/People/gian/data/idmapping_uniref50.new.db')
     cursor = conn.cursor()
     pids = i.pid.tolist()

     placeholders = ', '.join(['?' for _ in pids])
     sql_query = f"""
         SELECT *
         FROM EMBL
         WHERE [EMBL-CDS] IN ({placeholders});
     """
     uniref50DF = pd.read_sql_query(sql_query, conn, params=pids)
     uniref50DF['Uniref50'] = 'UniRef50_' + uniref50DF.uniprotid
     # Close the cursor and connection
     cursor.close()
     conn.close()

     # Return the result as a pandas Series
     return uniref50DF


def mview(seqobj,output='sequence.html', background='black', consensus = [100,90,80,70,60], organism=False, find=False):
    '''
    Using Mview to create a HTML file of the sequence algiment. 
    Need to add rule for replace the sequence object description column for organism name when organims is True
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cwd = os.getcwd()
    consensus = ','.join(map(str,consensus))
    f = ''
    if find:
        f = f'-find {find}'
    l2 = '-label2'
    if organism:
        l2 = ''
    if background =='black':
        b = 'black'
        t = 'white'
    else:
        b='white'
        t='black'

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        seqobj.to_file(f'{tmpdirname}/seqfile') 
        Popen(f'mview -in fasta -html head -consensus on -con_threshold {consensus} -label0 {l2} -label4 -label5 -pagecolor {b} -labcolor {t} -textcolor {t} -alncolor {b} -coloring any -con_coloring any {tmpdirname}/seqfile -css off -colormap CLUSTAL -bold -con_colormap CLUSTAL {f} > {cwd}/{output}',
              stdout=PIPE,
              shell=True
              ).communicate()
        
    return f'Consensus file saved on {cwd}/{output}' 

def rpsblast2table(psiblast_output, simple_profiledb=True):
    '''
    Function to read a raw string with a psi blast output and parse it to a table
    ''' 
    import pandas as pd
    import numpy as np

    l = pd.Series(psiblast_output.splitlines())
    f = l.str.startswith(('Query', 'Sbjct', ' Score','>', 'Query='))
    l = l.loc[f].to_frame()
    l.columns = ['text']
    l.loc[l.text.str.startswith("Query="), 'ID'] = l.query('text.str.startswith("Query=")').text.str.split(expand=True)[1]
    l.loc[l.text.str.startswith(">"), 'domain'] = l.query('text.str.startswith(">") == True').text.str.split(expand=True)[0].str.strip('>')
    l.loc[l.text.str.startswith(">PWD"), 'domain'] = l.query('text.str.startswith(">PWD") == True').text.str.split(expand=True)[2].str.strip(',')
    l.loc[l.text.str.startswith(">"), 'source'] = 'allprofiles'
    l.loc[l.text.str.startswith(">PWD"), 'source'] = 'Pfam'
    if simple_profiledb:
        l.loc[l.source=="allprofiles", 'domain'] = l.query('source =="allprofiles"').domain.str.split('.', expand=True)[0].str.split("_", expand=True)[0]
    l.loc[l.text.str.startswith(" Score"), 'evalue'] =  pd.to_numeric(l.query('text.str.startswith(" Score") == True').text.str.split(',',expand=True)[1].str.split('=', expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Query "), 'start'] = pd.to_numeric(l.query('text.str.startswith("Query ") == True').text.str.split(expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Query "), 'end'] = pd.to_numeric(l.query('text.str.startswith("Query ") == True').text.str.rsplit(expand=True, n=1)[1], errors='coerce')
    l.loc[l.text.str.startswith("Sbjct "), 'hit_start'] = pd.to_numeric(l.query('text.str.startswith("Sbjct ") == True').text.str.split(expand=True)[1], errors='coerce')
    l.loc[l.text.str.startswith("Sbjct "), 'hit_end'] = pd.to_numeric(l.query('text.str.startswith("Sbjct ") == True').text.str.rsplit(expand=True, n=1)[1], errors='coerce')
    l.evalue = l.evalue.ffill()
    l.domain = l.domain.ffill()
    l.ID = l.ID.ffill()
    l.source = l.source.ffill()
    l = l.query('text.str.startswith("Query ") or text.str.startswith("Sbjct ")')
    xx = l.groupby(['ID', 'domain', 'evalue', 'source']).agg({'hit_start':'min', 'hit_end':'max', 'start':'min', 'end':'max'}).reset_index()
    return xx

def phobius2table(phobius_output, short=True, add_evalue= 101e-4):
    '''
    Function to read a raw string with a psi blast output and parse it to a table
    Short : Ouputs a short df 
    add_evalue : add a dummy evalue to the short df to later be used in architecture functions.
    ''' 
    import pandas as pd
    l = pd.Series(phobius_output.splitlines())
    l = pd.Series(l)
    f = l.str.startswith(('FT', 'ID'))
    l = l.loc[f].to_frame()
    l.columns = ['text']
    l.loc[l.text.str.startswith("ID"), 'ID'] = l.query('text.str.startswith("ID")').text.str.split(expand=True)[1]
    l.ID = l.ID.ffill()
    l[['phobius', 'start', 'end','prediction']] = l.query('text.str.startswith("FT")').text.str.strip().str.split(expand=True)[[1,2,3,4]]
    l = l.query('text.str.startswith("FT")').iloc[:, 1:6]
    if short:
        l.phobius = l.phobius.replace({"SIGNAL": "SP", "TRANSMEM":"TM"})
        l.prediction = l['prediction'].fillna(l['phobius']).replace({"CYTOPLASMIC.":"IN", "NON":"OUT"})
        l = l.query('prediction in ["TM", "IN","OUT","SP"]')
        l = l[['ID','prediction','start', 'end']].rename({'prediction': 'phobius'}, axis=1)
        if add_evalue:
            l['evalue'] = add_evalue
    return l
def TMprediction(seqobj,
             predictior='phobius',
             cpu=96):
    '''
    Transmembrane predictor it accepts sequence object.
    Only working with phobius, more predictior to be added
    '''

    import tempfile
    from subprocess import Popen, PIPE
    from rotifer.devel.beta.sequence import sequence as sequence
    seqobj = seqobj.copy()

    # Shoul I remove gaps here??
    if predictior == 'phobius':
        with tempfile.TemporaryDirectory() as tmpdirname:
            # temporary save fasta sequence file
            if isinstance (seqobj, sequence):
                seqobj.to_file(f'{tmpdirname}/seqfile') 
                Popen(f'cat {tmpdirname}/seqfile | parallel -j {cpu} --pipe --recstart ">" phobius > {tmpdirname}/out 2> err',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            with open(f'{tmpdirname}/out') as f:
                    phobius_result = f.read()
        
    return (phobius_result) 

def rpsblast(seqobj,
             db=['allprofiles', 'pwld_new_pfam'],
             cpu=96):
    '''
    RPSblast it can accept sequence object. 
    DB can be suplied as python list, values accepted (NIH servers):
    allprofiles
    pwld_new_pfam
    pwld_pfam
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    import os
    if db == 'gian': #### Improvisation ################
        db =['/netmnt/vast01/cbb01/proteinworld/People/gian/projects/my_db/tmp/allprofiles', 'pwld_new_pfam']
    if isinstance(db,list):
        dbpath = f'{os.getenv("DATABASES")}/rpsdb/'
        with tempfile.TemporaryDirectory() as tmpdirname:
            for x in db:
                if "/" in x:
                    dbx = x
                else:    
                    dbx = f'{dbpath}{x}'
                # temporary save fasta sequence file
                if isinstance (seqobj, sequence):
                    seqobj =seqobj.copy()
                    seqobj.to_file(f'{tmpdirname}/seqfile') 
                    Popen(f'cat {tmpdirname}/seqfile | splishrps -d {dbx} -a {cpu} >> {tmpdirname}/out', stdout=PIPE, shell=True).communicate()
                    with open(f'{tmpdirname}/out') as f:
                            blast_r = f.read()
                else:
                    print('please supply a sequence obj')
        return blast_r 
    else:
        with tempfile.TemporaryDirectory() as tmpdirname:
            # temporary save fasta sequence file
            if isinstance (seqobj, sequence):
                seqobj =seqobj.copy()
                seqobj.to_file(f'{tmpdirname}/seqfile') 
                Popen(f'cat {tmpdirname}/seqfile | splishrps -d {db} -a {cpu} > {tmpdirname}/out', stdout=PIPE, shell=True).communicate()
                with open(f'{tmpdirname}/out') as f:
                        blast_r = f.read()
            else:
                print('please supply a sequence obj')

        return blast_r 

def get_plen(pidlist):
    """
    small function to get plen of proteins, can receive list, pd.Series or a sequence objct
    """
    import pandas as pd
    import rotifer.devel.beta.sequence as rdbs
    if isinstance(pidlist, pd.Series):
        pidlist = pidlist.drop_duplicates().tolist()
        
    if isinstance(pidlist, rdbs.sequence):
        results = pidlist
    else:
        results = rdbs.sequence(pidlist)
    
    results = results.df[['id','length']].rename({'id':'ID', 'length': 'full_length'}, axis=1)
    return results

def domtable(dom_table, seqobj=False):
    import pandas as pd
    import numpy as np
    from rotifer.devel.alpha import gian_func as gf
    copy = dom_table.copy()
    if seqobj:
        dom_table = dom_table.merge(gf.get_plen(seqobj))
    else:    
        dom_table = dom_table.merge(gf.get_plen(dom_table.ID.unique().tolist()))

    dom_table = dom_table.sort_values(['ID', 'start', 'end'])
    dom_table['region'] = np.where(dom_table.ID.shift(+1) == dom_table.ID, dom_table.start - dom_table.end.shift(+1),0)
    dom_table['region2'] =  (dom_table.start - dom_table.region) +1
    dom_table['region3'] =  (dom_table.start - dom_table.region2)
   # unk_central_dom = dom_table.query('region3 >=5')
   # unk_n_term_dom = dom_table.query('region3 ==-1 and region2 >=5')
    dom_table['region4'] = np.where(dom_table.ID.shift(-1) == dom_table.ID,0, dom_table.full_length - dom_table.end)
    dom_table['region5'] = dom_table.start.shift(-1)-1
    dom_table['region5'],dom_table['region6'] = np.where(dom_table.region4 > 0 ,[dom_table.end, dom_table.full_length], [ dom_table.end, dom_table.region5])
    dom_table['region7'] = dom_table.region6 - dom_table.region5
    dom_n = dom_table.drop_duplicates('ID')
    dom_n.end = dom_n.start
    dom_n.start = 1
    dom_n = dom_n[['ID', 'start', 'end']]
    dom_table = dom_table[['ID', 'region5' ,'region6', 'region7']]
    dom_table.columns = ['ID', 'start', 'end', 'rlen']
    dom_n['rlen'] = dom_n.end - dom_n.start
    dom_table = pd.concat([dom_table, dom_n]).sort_values(['ID', 'start', 'end'])
    xx = dom_table.query('rlen >=2')
    xx.domain =['unk']
    dom_file = copy[['ID', 'domain', 'start', 'end']]
    xx['domain'] = 'unk'
    final = pd.concat([dom_file,xx]).sort_values(['ID', 'start','end'])[['ID','domain', 'start','end']]

    return (final)

def draw_architecture(input_file,
                      file, id_list=True,
                      size_median=True,
                      not_spaces=False,
                      domain_rename=True,
                      db=['allprofiles', 'pwld_new_pfam'],
                      evalue=1e-3,
                      color_dict=None,
                      shape_dict = None,
                      round_dict = None):

    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import seaborn as sns
    import yaml
    import pandas as pd
    from rotifer.devel.beta.sequence import sequence
    import rotifer.devel.alpha.gian_func as gf
    from rotifer.interval import utils as riu
    from rotifer.core import functions as rcf
    from io import StringIO

    svg_dict = {'line_svg' : rcf.findDataFiles(":templates/line.svg"),
                'ANKs' : rcf.findDataFiles(":templates/Ank.svg"),
                'Bps' : rcf.findDataFiles(":templates/Bet_propelers.svg"),
                'HEATs' : rcf.findDataFiles(":templates/HEAT.svg"),
                'LRRs' : rcf.findDataFiles(":templates/LRR.svg"),
                'Pps' : rcf.findDataFiles(":templates/Pps.svg"),
                'Sel1' : rcf.findDataFiles(":templates/Sel1.svg"),
                'TPRs' : rcf.findDataFiles(":templates/TPR.svg")
                }

    #repeats = ['ANKs',  'Bps','HEATs',  'LRRs',  'Pps'  'Sel1',  'TPRs']
    # As reference Nature's standard figure sizes are 89 mm (3.50 inches) wide (single column) and 183 mm (7.20472 inches) wide (double column)
    # If a protein with 5000 aminoacids would fit a whole line of double colum figure, each aa would then use 0.001440944 inches.
    # If I use a PDF with 8 inches and 0.5 margin I would have 7 inches available
    # Therefor we could use this number as conversion factor to draw the domains shapes and rouglth estimate the size to write within it.
    #The font size is typically measured in points. One point is equal to 1/72 of an inch

    # Example : a domain of 300 aa would have a shape of 300 *0.001440944 = 0.4322832 inches
    # If I use font of size 5, each character would 
    # Font size (in inches)= 72/5 ≈0.069 inches
    #For proportional fonts like Arial, the average width of a character is typically about 0.6 times the font size.
    #Average character width (in inches)=0.069×0.6≈0.041 inches

    #Therefore to use calculate the number of characters that fits 300 aa domain:
    # 0.4322832/0.041 == ~10


    # Figure to scale (reccomended to keep 1 to produce final figure in adobe (best scale for publication, but needs to be proper adjusted the names)
    ##### To create PDF, the 4 works preatty good
    if id_list:
        id_list=input_file
        seqobj = sequence(id_list)
        seqobj.arch = gf.rpsblast(seqobj, db=db)
        tsv =  gf.rpsblast2table(seqobj.arch)
        tsv = tsv.query('evalue <= 0.1')
        seqobj.phobius = gf.TMprediction(seqobj)
        pt = gf.phobius2table(seqobj.phobius).rename({"phobius":"domain"}, axis=1)
        tt = pd.concat([tsv,pt])
        tt.start = tt.start.astype(int)
        tt.end = tt.end.astype(int)
        no = riu.filter_nonoverlapping_regions(tt.query('domain not in ["IN", "OUT"]').query(f'domain == "TM" or evalue <={evalue}')  ,**{**riu.config['archtsv'], 'maximum_overlap':0.4})
    else:
        if isinstance(input_file,pd.DataFrame):
            dt = input_file
        else:    
            dt = pd.read_csv(input_file, sep="\t", names=['ID','arch','dinfo'])
        dt.dinfo = dt.dinfo.str.split(',')
        dt = dt.explode('dinfo')
        dt[['start', 'end', 'domain']] = dt.dinfo.str.split('\.\.|\&', expand=True)
        dt = dt.sort_values(['ID', 'start', 'end'])
        dt.end = dt.end.astype(int)
        dt.start = dt.start.astype(int)
        no = dt
        seqobj = sequence(dt.ID.unique().tolist())
        id_list=False
 
        
    # get the organism name from the fasta description:
    seqobj.df['organism'] = seqobj.df.description.apply(gf.extract_organism_from_description)
    organism_dict = seqobj.df.set_index("id").organism.to_dict()


    scale_figure = 2
    domain_height = scale_figure/10
    font_point = float(0.013837) * 0.6  # 1/72 (size of one point character in inches) * actual size (avg 0.6 the real size)
    aa_scale = float(0.001440944) 

    shapes = ['box',
              'octagon',
    #          'circle',
              'hexagon',
              'ellipse',
              'diamond']
    #          'doubleoctagon',
    #          'tripleoctagon']
    shapes = pd.DataFrame(shapes, columns=['shapes']).reset_index().rename({'index':'shape_rank'}, axis=1)
    colors = sns.color_palette("pastel").as_hex() +sns.color_palette("deep").as_hex()  + sns.color_palette('muted').as_hex() + sns.color_palette('colorblind').as_hex()
    colors = pd.DataFrame(colors, columns=['colors']).reset_index().rename({'index':'color_rank'}, axis=1)
    color_code = shapes.join(colors, how='cross').join(pd.DataFrame(['rounded,filled', 'filled']), how='cross').rename({0:'rounded'}, axis=1).sort_values(['rounded','shape_rank','color_rank'], ascending=[False, True, True]).reset_index(drop=True)
    color_code = color_code.query('~(shapes =="ellipse" and rounded =="rounded,filled")')


    font_size= 4
    d = gf.domtable(no, seqobj=seqobj)
    d = d.rename({"ID":'pid'}, axis=1)

    if domain_rename:
        domain_dict = rcf.loadConfig(rcf.findDataFiles(":data/domain_rename.yaml"))
        domain_dict =  {v:k for k in domain_dict.keys() for v in domain_dict[k] }
        d['domain'] = d['domain'].replace(domain_dict)
    #due Some strange behaivour I will strip the columns domain:
    d['domain'] = d['domain'].str.strip()
    d['size']  = d.end - d.start
    d['esc'] = d['size'] * (aa_scale*scale_figure)
    #d.domain = d.domain.replace(domain_dict['Domain'])
    d['dl'] = d.domain.str.len()
    d['f_space']= round(d.esc/(font_point*font_size))
    d = d[d.f_space >= 1]
    d.f_space = d.f_space.astype(int)
    d['pid_order'] = d.pid.map(d.pid.drop_duplicates().reset_index(drop=True).reset_index().set_index('pid')['index'].to_dict())
    d = d.reset_index(drop=True).reset_index()

    if size_median:
        d.loc[d.domain !='unk', 'esc'] = d.loc[d.domain !='unk'].groupby('domain').esc.transform('median')

    # Adding custoum color and shape for TMs
    ################################################
    ll = d[~d.domain.isin(['unk', 'TM'])]
    ll = ll[~ll.domain.isin(svg_dict.keys())].domain.value_counts().reset_index()
    
    oo = ll.join(color_code)
    ## Adding the TM color and shape
    oo = pd.concat([pd.DataFrame(['TM', 0, 1, 'box', 0, '#FF912B', 'rounded,filled'], index=oo.columns).T, oo])
    ## fillling color and shape for domains not much represented on the ds
    oo.colors = oo.colors.fillna("#808080")
    oo.shapes = oo.shapes.fillna("ellipse")
    oo.rounded = oo.rounded.fillna("filled")
    sd = oo[['domain', 'shapes']].set_index('domain').shapes.to_dict()
    if shape_dict:
        sd = {**sd, **shape_dict}
    cd = oo[['domain', 'colors']].set_index('domain').colors.to_dict()
    if color_dict:
        cd = {**cd, **color_dict}
    rd = oo[['domain', 'rounded']].set_index('domain').rounded.to_dict()
    if round_dict:
        rd = {**rd, **round_dict}

    ######################################
    #Setting font size based on the shapes
    d['font_size'] = font_size
    d['domain_char'] = d.domain.str.len()
    d['write_size'] = d.domain_char * font_point * font_size
    d['over_write'] = d.esc / d['write_size']
    d.loc[d.over_write < 1, 'font_size'] =   d.loc[d.over_write < 1].font_size * d.loc[d.over_write < 1].over_write
    fsd = d.drop_duplicates('domain').set_index('domain').font_size.to_dict()


    draw_scale=True
    scale_size = 100 * (aa_scale*scale_figure)
    if not_spaces:
        d = d.query('domain !="unk"')
    gb = d.groupby('pid_order')
    li =[]
    A = pgv.AGraph()
    #Setting the size of output, based in the number of proteins:
    A.graph_attr['size'] = f'12,{d.pid.nunique()}'
    A.graph_attr['margin'] = '0.5'
    if draw_scale:
        l =['scale_text', 'scale']
        A.add_node('scale_text',
                   label='scale',
                   shape='none',
                   fontsize=font_size +1,
                   fontname='Arial')
        A.add_node('scale',
                   label='100',
                   width=scale_size,
                   fixedsize='true',
                   nodesep=0,
                   fontsize=font_size,
                   #shape='box',
                   image=svg_dict['line_svg'],
                   #style = 'filled',
                   height='0.00001',
                   #fillcolor='black',
                   penwidth=0,
                   labelloc='t',
                   #color='white',
                   labeldistance=1.5)
        li.append(l)

    for x in gb:
        l =[]
        for y,z in x[1].iterrows():
            if len(l) < 1:
                l.append(z['pid'])
                A.add_node(z['pid'],
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="5"></TD></TR>'
                               f'<TR><TD HEIGHT="5"></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="5">{x[1].pid.unique()[0]}</FONT></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial italic" POINT-SIZE="4">{organism_dict[x[1].pid.unique()[0]]}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=font_size +1,
                           fontname='Arial',
                           labelloc="t")

            if z.domain =='unk':
                A.add_node(z['index'],
                           label='',
                           width=z['esc'],
                           fixedsize='true',
                           nodesep=0,
#                           shape='box',
                           image=svg_dict['line_svg'],
                           height='domain_height',
                           penwidth=0)
                l.append(z['index'])
            elif z.domain in svg_dict.keys():
                A.add_node(z['index'],
                           label=z.domain,
                           width=z['esc'],
                           fixedsize='true',
                           nodesep=0,
                           image=svg_dict[z.domain],
                           fontsize=fsd[z['domain']],
                           height='domain_height',
                           penwidth=0,
                           imagescale="true",
                           fontname='Arial')
                l.append(z['index'])
            else:
                A.add_node(z['index'],
                           label=z.domain,
                           style=rd[z['domain']],
                           width=z['esc'],
                           fixedsize='true',
                           height=domain_height,
                           color=cd[z['domain']],
                           fontsize=fsd[z['domain']],
                           fillcolor=cd[z['domain']],
                           shape=sd[z['domain']],
                           nodesep=0,
                           penwidth=1,
                           fontname='Arial')
                l.append(z['index'])
        li.append(l)
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")

    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.02)
    A.graph_attr["fontnames"] = "svg"
    A.draw(file, prog="dot")


def extract_organism_from_description(text):
    import re
    match = re.search(r'\[(.*?)\]', text)
    if match:
        content = match.group(1)
        words = content.split()
        return ' '.join(words[:2])
    return None

import numpy as np
import pandas as pd
from seaborn import light_palette

def insert_na_rows_at_indexes(df, indexes, insert_position='upstream', fill_rules=None):
    """
    Inserts rows with NaN values into a DataFrame at the specified indexes
    and fills NaN values in the specified columns according to provided rules.

    Parameters:
    df : pandas.DataFrame
        The DataFrame to modify.
    indexes : list of int
        The indexes where rows will be inserted.
    insert_position : str, optional
        Determines where the row will be inserted relative to the index ('upstream' or 'downstream').
    fill_rules : dict, optional
        A dictionary where keys are column names and values are either 'forward', 'backward',
        or a specific string to fill NaN values with.

    Returns:
    new_df : pandas.DataFrame
        The modified DataFrame with inserted rows and filled values.
    """
    
    sorted_indexes = sorted(indexes, reverse=True)

    # Create a DataFrame for the NaN row
    na_row_df = pd.DataFrame([{col: np.nan for col in df.columns}], columns=df.columns)

    # Convert DataFrame to a list of dictionaries for easier manipulation
    df_list = df.to_dict(orient='records')

    for index in sorted_indexes:
        # Adjust index based on the 'insert_position'
        if insert_position == 'downstream':
            index += 1  # Insert after the index

        # Insert NaN row into the list at the adjusted index
        df_list.insert(index, na_row_df.iloc[0].to_dict())

    # Convert the list of dictionaries back to DataFrame
    new_df = pd.DataFrame(df_list)

    # Apply the fill rules based on the fill_rules dictionary
    if fill_rules:
        for col, rule in fill_rules.items():
            if rule == 'forward':
                # Fill NaN with the value from the previous row (shift forward)
                new_df[col] = new_df[col].fillna(method='ffill')
            elif rule == 'backward':
                # Fill NaN with the value from the next row (shift backward)
                new_df[col] = new_df[col].fillna(method='bfill')
            else:
                # Fill NaN with the specified string
                new_df[col] = new_df[col].fillna(rule)

    return new_df


def insert_na_rows_at_indexes_2(df, indexes, columns_to_fill=None, fill_direction=None, insert_position='upstream'):

    import pandas as pd
    import numpy as np

    """
    Inserts rows with NaN values into the DataFrame at specified indexes and fills NaN values 
    based on the given columns and directions.

    Parameters:
    df (pd.DataFrame): The original DataFrame.
    indexes (list of int): A list of indexes at which to insert new rows with NaN values.
    columns_to_fill (list of str): List of columns where NaN values should be filled.
    fill_direction (list of str): List of directions ('backward' or 'forward') for each column 
                                  in 'columns_to_fill'. Determines how to fill NaN values.
    insert_position (str): Specifies where to insert the NaN row ('upstream' for before the index,
                           'downstream' for after the index).

    Returns:
    pd.DataFrame: A new DataFrame with the NaN rows inserted and optionally filled.
    """
    # Ensure indexes are sorted in descending order to handle multiple insertions correctly
    sorted_indexes = sorted(indexes, reverse=True)

    # Create a DataFrame for the NaN row
    na_row_df = pd.DataFrame([{col: np.nan for col in df.columns}], columns=df.columns)

    # Convert DataFrame to a list of dictionaries for easier manipulation
    df_list = df.to_dict(orient='records')

    for index in sorted_indexes:
        # Adjust index based on the 'insert_position'
        if insert_position == 'downstream':
            index += 1  # Insert after the index

        # Insert NaN row into the list at the adjusted index
        df_list.insert(index, na_row_df.iloc[0].to_dict())

    # Convert the list of dictionaries back to DataFrame
    new_df = pd.DataFrame(df_list)

    # If columns to fill and directions are provided
    if columns_to_fill and fill_direction:
        for col, direction in zip(columns_to_fill, fill_direction):
            if direction == 'forward':
                # Fill NaN with the value from the previous row (shift forward)
                new_df[col] = new_df[col].fillna(method='ffill')
            elif direction == 'backward':
                # Fill NaN with the value from the next row (shift backward)
                new_df[col] = new_df[col].fillna(method='bfill')

    return new_df

def pick_colors(df,column='dom',top_domains=20, pallets = ['pastel','deep','muted','colorblind'], own_dict=False, domain_rename=True, second_rename=False):
    import seaborn as sns
    import pandas as pd
    from rotifer.core import functions as rcf
    from packaging import version
    df = df.copy()
    ### Creating a pd.Series with  hex colors from sns pallet function:
    colors = []
    for pallet in pallets:
        colors = colors +  sns.color_palette(pallet).as_hex() 
    colors = pd.DataFrame(colors, columns=['colors']).reset_index().rename({'index':'rank'}, axis=1)


    if domain_rename:
        domain_dict = rcf.loadConfig(rcf.findDataFiles(":data/domain_rename.yaml"))
        domain_dict =  {v:k for k in domain_dict.keys() for v in domain_dict[k] }
        df[column] = df[column].replace(domain_dict)
    if second_rename:
        df[column] = df['column'].replace(second_rename)

        
    ### Using the dataframe to rank the most commons domains to be colored
    if version.parse(pd.__version__) >= version.parse("2.0.0"):
        domain_rank = df[column].value_counts().reset_index().rename({column:'domain'}, axis=1).reset_index().query( 'domain not in [ "-", "?", "|", "", " ", "LIPO", "TM", "SIG", "SP"]').rename({'index':'rank'}, axis=1)[['rank', 'domain']]
    else:
        domain_rank = df[column].value_counts().reset_index().rename({'index':'domain'}, axis=1).reset_index().query( 'domain not in [ "-", "?", "|", "", " ", "LIPO", "TM", "SIG", "SP"]').rename({'index':'rank'}, axis=1)[['rank', 'domain']]

    domain_rank = domain_rank.merge(colors, how='left')
    domain_rank['colors'] = domain_rank['colors'].apply(lambda x: x if not pd.isna(x) else domain_rank['colors'].dropna().sample(1).values[0])

    if top_domains:
        domain_rank.loc[top_domains:, 'colors'] = "#D3D3D3"
    #Creating the color_dict:
    domain_rank = domain_rank.merge(colors, how='left')
    color_dict = domain_rank.set_index('domain').colors.to_dict()
    if own_dict:
        color_dict = {**color_dict, **own_dict}
    return color_dict


def split_domain(df,
                 column='arch',
                 fill ='?',
                 domain_rename=True,
                 remove_tm=False,
                 query_same_direction=False,
                 check_duplicates=False,
                 reverse_annotaion=False):
    '''
    It will use the query as anchor to delect the given number of after and before genes to
    create a dataframe of domains.
    One can use a dcitionary of domain to select the ones that would be kept and rename it
    '''
    import numpy as np
    import pandas as pd
    from rotifer.core import functions as rcf
    from rotifer.devel.alpha import gian_func as gf

    df = df.copy()
    df.reset_index(drop=True, inplace=True) ### Add recently 20241202
    # Asuring to not have duplicated protein in the input DF:
    if check_duplicates:
        domdf = df.drop_duplicates(['nucleotide','start','end']).reset_index(drop=True).query('type == "CDS"').copy()
    else:
        domdf = df.copy()

    # Creating an order for block_id domain representation:
    domdf['blockp'] = pd.Series(np.where(domdf.block_id == domdf.block_id.shift(1), 1,0))
    domdf.blockp = domdf.groupby(['block_id'])['blockp'].cumsum()

    if query_same_direction:
        bi = domdf.query('query ==1').drop_duplicates('block_id').query('strand == -1').block_id
        mask = domdf['block_id'].isin(bi)
        domdf.loc[mask, 'strand']  = domdf.loc[mask, 'strand'] * -1
        domdf.loc[mask, 'blockp'] = domdf.loc[mask].groupby('block_id')['blockp'].transform(lambda x: x[::-1].values)
        domdf = domdf.sort_values(['block_id','blockp']).reset_index(drop=True)

    domdf['dom'] = domdf[column].str.split('+')
    if reverse_annotaion:
        pass
    else:
        domdf.loc[domdf.strand  == -1, 'dom'] = domdf.loc[(domdf.strand  == -1)].dom.apply(lambda x: x[::-1])
    domdf = domdf.explode('dom').reset_index(drop=True)
    domdf['domp'] = pd.Series(np.where(domdf.pid == domdf.pid.shift(1), 1,0))
    domdf.domp = domdf.groupby(['pid','block_id'])['domp'].cumsum()
    domdf['blockp'] = pd.Series(np.where(domdf.block_id == domdf.block_id.shift(1), 1,0))
    domdf.blockp = domdf.groupby(['block_id'])['blockp'].cumsum()


    if fill:
        domdf.dom = domdf.dom.fillna(fill)

    if domain_rename:
        domain_dict = rcf.loadConfig(rcf.findDataFiles(":data/domain_rename.yaml"))
        domain_dict =  {v:k for k in domain_dict.keys() for v in domain_dict[k] }
        if isinstance(domain_rename, dict):
            domain_dict = {**domain_dict, **domain_rename}
        domdf['dom'] = domdf['dom'].replace(domain_dict)


    # removing tm sig
    if remove_tm:
        domdf['tmc'] = np.where(domdf.dom.isin(["TM","SIG"]),0,1)
        only_tm_sig = list(domdf.groupby('pid').tmc.sum().where(lambda x : x==0).dropna().index)
        domdf.dom = np.where(domdf.pid.isin(only_tm_sig), fill, domdf.dom)
        domdf = domdf.query('dom no in ["TM", "SIG"]')
     
    #reorder to keep the same input order:
    domdf = domdf.sort_values(['block_id','blockp','domp']) # Remove pid fomr the oder


    return domdf


def operon_fig_bkp(df,
                domain_rename=True,
                output_file='operon_fig_out.pdf',
                domain_column='arch',
                query_same_direction=True,
                color_dict=None,
                top_domains=10,
                sort_by='classification',
                height=0.28,
                f=2,
                fontsize=4,
                query_asterix=False,
                query_color='red',
                check_duplicates=True,
                light_palette='90',
                reverse_annotaion=False):

    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import pandas as pd
    from rotifer.devel.alpha import gian_func as gf
    from rotifer.core import functions as rcf
    import seaborn as sns

    if check_duplicates:
        te = df.query('type =="CDS"').copy()
        te[domain_column] =  te[domain_column].fillna('?')
    else:
        te = df.copy()

    te = gf.split_domain(te,
                         domain_rename=domain_rename,
                         column=domain_column,
                         query_same_direction=query_same_direction,
                         check_duplicates=check_duplicates,
                         reverse_annotaion=reverse_annotaion)
    te['shape'] = 'rectangle'
    # Removing duplicates sequencial domain in same protein
    #te = te.loc[~((te.pid.shift() == te.pid) & (te.dom.shift() == te.dom))].reset_index(drop=True)
    te = te.reset_index(drop=True)

    l = te.drop_duplicates(['block_id','pid'], keep='first').query('strand ==-1').index
    r = te.drop_duplicates(['block_id','pid'], keep='last').query('strand ==1').index
    te.loc[r,'shape'] = 'rarrow'
    te.loc[l,'shape'] = 'larrow'

    te['height'] = np.where(te['shape'].isin(['larrow', 'rarrow']),  height, height/f)
    te['style'] = 'filled' 
    # insert an space between oposite strands.
    add_rows = te[(te.strand > te.strand.shift(1) )& (te.block_id == te.block_id.shift(1))].index.tolist()
    te = gf.insert_na_rows_at_indexes(te, add_rows, fill_rules={'pid':'forward','block_id':'forward'})
    te.loc[te.nucleotide.isna(),['dom','domp','shape','style', 'height', 'width']] = ['', 0, 'rectangle','', 0.2, 2] 
    te = te.reset_index(drop=True).reset_index()
    if color_dict == None:
        color_dict = gf.pick_colors(te, column = 'dom', top_domains=top_domains)

    if light_palette:
         color_dict = {key: value + light_palette for key, value in color_dict.items()}


    color_dict = {**color_dict,
                  '' :"white",
                  "-" :"#D3D3D3",
                  " " :"#D3D3D3",
                  "?":"#D3D3D3",
                  "TM" :"#D5B60A",
                  "LIPO":"Blue",
                  "LP":"Blue",
                  "SP":"Red",
                  "SIG":"Red"}
    
    te['color'] = te.dom.replace(color_dict)
#    te.dom = te.dom.replace({'TM':'', 'LP':'', 'LIPO':'', 'SIG':'', 'SP': ''})
#    te.loc[(te.dom=='') & (te['shape'].isin(['larrow','rarrow'])), 'dom'] = ' '
    te.loc[(te.dom=='') & (te['shape'].isin(['larrow','rarrow'])), 'color'] = '#D3D3D3'
####### Adding new row to represent empty proteins with TM, SIG or LIPO
##### A lot of code to make sure the order, color and shapes  of the added rows are correct.
    ra = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['rarrow']))].index
    te = gf.insert_na_rows_at_indexes(
       te,
       ra,
       fill_rules={
           'nucleotide':'backward',
           'pid':'backward',
           'block_id':'backward',
           'arch':'backward',
           'dom':'backward',
           'color':'backward',
           'shape':'rectangle',
           'height':'backward',
           'strand':'backward',
           'style':'backward'})
#    te = te.reset_index(drop=True) 
    ra = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['rarrow']))].index
    te.loc[ra,['shape','color']] =[['rarrow', '#D3D3D3']]
    la = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['larrow']))].index
    te = gf.insert_na_rows_at_indexes(
       te,
       la,
       insert_position='downstream',
       fill_rules={
           'nucleotide':'forward',
           'pid':'forward',
           'block_id':'forward',
           'arch':'forward',
           'dom':'forward',
           'color':'forward',
           'shape':'rectangle',
           'height':'forward',
           'strand':'forward',
           'style':'forward'})

#    te = te.reset_index(drop=True) 
    #Replacing the - or ? marker for empty spaces to decrease the noise on the figure.
    te['bcolor'] = te.color
    te['penwidth'] = 1
    te['width'] = te.dom.str.len()/40
    te.loc[te.dom =="TM", ['bcolor', 'width', 'penwidth']] =['black', 0.03, 0.5]
    te.loc[te.dom.isin(["SIG", "SP"]), ['bcolor', 'width', 'penwidth']] =['black', 0.03, 0.5]
    te.loc[te.dom =="PSE", ['color', 'bcolor']] =['#D3D3D340', '#D3D3D3']

    te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['larrow','rarrow'])), ['dom', 'color', 'bcolor']] = [' ', '#D3D3D3', '#D3D3D3']

    if query_asterix:
        to_a = te.query('query ==1').query("dom not in ['TM','LP','LIPO','SIG','SP']").drop_duplicates(subset=['pid'], keep='last').index ###
        te.loc[to_a, 'dom'] = te.loc[to_a, 'dom'] +'*' ###
    if query_color:    
        to_a2 = te.query('query ==1').query("dom not in ['TM','LP','LIPO','SIG','SP']").index ###
        te.loc[to_a2, 'bcolor'] = query_color ###
    te.dom = te.dom.replace({'TM':'', 'LP':'', 'LIPO':'', 'SIG':'', 'SP': ''})
    te.dom = te.dom.replace({'?':' ', '-':' '})
    te.loc[te.dom == "", ['height']] = 0.185
    ####Creating a column with the pid of the query seed to later be used in given a header of the figure.
    te['query_pid'] = te.groupby('block_id')['pid'].transform(lambda x: x[te['query'] == 1].iloc[0] if any(te['query'] == 1) else np.nan)   
    ##### Reseting the index to use the new index as name of each node in the figure
    del te['index']
    te = te.reset_index()

    if sort_by:
        ordered_list = te.sort_values([sort_by,'block_id','blockp']).block_id.drop_duplicates().tolist() 
    else:
        ordered_list = te.block_id.drop_duplicates().tolist() 
    

    gb = te.groupby('block_id')
    li=[]
    A = pgv.AGraph()
   # return gb
    for x in ordered_list:
        block_id_df = gb.get_group(x)
        l =[]
        for y,z in block_id_df.iterrows():
            if len(l) < 1:
                l.append(z['block_id'])
                A.add_node(z['block_id'],
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="4">{block_id_df.query_pid.unique()[0]}</FONT></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="4">{block_id_df.block_id.unique()[0]}</FONT></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial italic" POINT-SIZE="4">{block_id_df.organism.unique()[0]}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Arial')
            A.add_node(z['index'],
                           label=z.dom,
                           shape=z['shape'],
                           width=z['width'],
                           style=z['style'],
                           height=z['height'],
                           color=z.bcolor,
                           penwidth=z['penwidth'],
                           fillcolor=z.color,
                           labeldistance=0.5,
                           fontsize=fontsize)
            l.append(z['index'])

        li.append(l)
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")


    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.02)
    A.draw(output_file, prog="dot")
    return te


def compact_to_df(compact,sep='Tab', columns= ['pid', 'compact', 'organism', 'pid_compact']):
    """
    convert TASS annotated compact represantation to DF
    """
    import pandas as pd
    if sep=='Tab':
        x = pd.read_csv(compact, sep="\t", names=columns)
    else:
        x = pd.read_csv(compact, delimiter=r"\s\s+", names=columns)
    x = x.drop_duplicates(subset=['pid'])
    x['arch'] = x.compact.str.split('->')

    x = x.explode('arch')
    x2 = x[['arch','organism']]
    x2['strand'] = -1
    x2.loc[~x2.arch.fillna('').str.contains('<-'), 'strand'] = 1
    x2.arch = x2.arch.str.split('<-')
    x2 = x2.explode('arch')
    x2 = x2.query('arch not in ["", "||"]')
    x2.arch = x2.arch.str.replace('\|\|', '|--')
    x2.arch = x2.arch.str.split('|')
    x2 = x2.explode('arch')
    x2.loc[x2.arch.fillna(" ").str.startswith('--'), 'strand'] = 1
    x2.arch = x2.arch.str.strip('--')
    x2 = x2.reset_index().rename({'index':'block_id'}, axis=1)
    x2['block_id'] = x2.block_id.replace(x.pid.to_dict())
    x2 = x2.reset_index().rename({'index':'pid'}, axis=1)
    x2['query'] = 0
    x2.loc[x2.arch.fillna(' ').str.endswith('*'), 'query'] = 1
    x2.arch = x2.arch.str.strip('*')
    x2['nucleotide'] = 'dummy_collumn' 
    return x2


def check_spacing(text):
    import re
    """
    Check if the string contains a tab or two or more spaces immediately
    after the first word in the input text.

    Parameters:
    -----------
    text : str
        The input string to be checked.

    Returns:
    --------
    str
        A message indicating whether a tab, two or more spaces, or neither 
        were found after the first word:
        - "Tab found" if a tab ('\t') follows the first word.
        - "Two or more spaces found" if two or more spaces follow the first word.
        - "Neither found" if neither is present after the first word.
    """
    # Check for a tab directly after the first word
    if re.match(r'^\S+\t', text):
        return "Tab"
    
    # Check for two or more spaces directly after the first word
    if re.match(r'^\S+  +', text):
        return "Two_spaces"
    
    return "No tabular format found(with spaces/Tab), check your input"
def pid2tax(pidlist,full=False):
    import pandas as pd
    from rotifer.db import ncbi
    ic = ncbi.IPGCursor()
    if isinstance(pidlist,pd.Series):
        pidlist = pidlist.unique().tolist()
    ipgs = ic.fetchall(pidlist)
    #ipgs = pd.concat(ipgs) Removing this line because it seems no reason to be here and is breaking the code
    assemblies_id = ipgs.assembly.unique().tolist()
    a = ncbi.assemblies(baseurl="/am/ftp-genomes/ASSEMBLY_REPORTS", taxonomy=False)
    a.rename({'#assembly_accession': 'assembly_accession'}, axis=1, inplace=True)
    a = a.query('assembly_accession in @assemblies_id')
    tc =   ncbi.TaxonomyCursor()
    taxonomies = tc.fetchall(a.taxid.tolist())
    a = a[['assembly_accession','taxid']]
    result = a.merge(taxonomies,
                     how='left').merge(ipgs[['pid','nucleotide','start','stop', 'assembly']],
                                       left_on='assembly_accession',
                                       right_on='assembly',
                                       how='left'
                                       )
    result['position'] = result.nucleotide + ":" + result.start.astype(str) + '-' +  + result.stop.astype(str)                 

    return result

def operon_fig2(df,
                domain_rename=True,
                output_file='operon_fig_out.pdf',
                domain_column='arch',
                query_same_direction=True,
                color_dict=None,
                shape_dict=None,
                top_domains=10,
                sort_by='pid',
                height=0.35,
                f=2,
                fontsize=10,
                query_asterix=False,
                query_color='red',
                highlight_domain=False,
                check_duplicates=True,
                light_palette='90',
                labeldistance=0.5,
                margin=0.03,
                unknow='  ?  ',
                id_position='side',
                comments=False,
                to_string=False,
                reverse_annotaion=False,
                lineage=False,
                yaml_file=False,
                architecture=False,
                genome_coordinates=True
                ):
    """
    Highlight domain should be dictnionary
    """
    import pygraphviz as pgv
    import numpy as np
    import pandas as pd
    from rotifer.devel.alpha import gian_func as gf
    import textwrap
    import io


    penwidth = 1
    TM_width = 0.03 
    space_width = 0.1 
    if fontsize == 4:
        labeldistance = 0.5
        margin = 0.015
        height = 0.15
        TM_width = 0.02
        space_width = 0.025
        penwidth = 0.4

    if yaml_file:
        yaml_df =gf.yamldomain2df(yaml_file)
        yaml_df.Include = yaml_df.Include.astype(int)
        yaml_df.loc[yaml_df.Include ==0, 'Display_name'] = '-'
        yaml_rename_dict = yaml_df.Display_name.to_dict()
        yaml_df.loc[yaml_df.Include ==0, 'color'] = '#D3D3D3'
        yaml_color_dict = yaml_df.color.to_dict()
        if isinstance(domain_rename, dict):
            domain_rename = {**domain_rename, **yaml_rename_dict}
        else:
            domain_rename = yaml_rename_dict


    if check_duplicates:
        te = df.query('type =="CDS"').copy()
        te[domain_column] =  te[domain_column].fillna('?')
    else:
        te = df.copy()

    te = gf.split_domain(te,
                         domain_rename=domain_rename,
                         column=domain_column,
                         query_same_direction=query_same_direction,
                         check_duplicates=check_duplicates,
                         reverse_annotaion=reverse_annotaion)
    te['shape'] = 'rectangle'
    # Removing duplicates sequencial domain in same protein
    #te = te.loc[~((te.pid.shift() == te.pid) & (te.dom.shift() == te.dom))].reset_index(drop=True)
    te = te.reset_index(drop=True)
    if architecture:
        te['block_id'] = te.apply(lambda x: x.block_id.split(":")[0]+":"+ str(x.start)+"-"+ str(x.end), axis=1)
        te['query'] = 1
        TM_width = TM_width *2.5
        if shape_dict:
            # Setting a default shape"
            te['shape'] = "oval"
            shape_dict = {**te.set_index('dom')['shape'].to_dict(), **shape_dict}
            te['shape'] = te.dom.map(shape_dict)
    else:    
        l = te.drop_duplicates(['block_id','pid'], keep='first').query('strand ==-1').index
        r = te.drop_duplicates(['block_id','pid'], keep='last').query('strand ==1').index
        te.loc[r,'shape'] = 'rarrow'
        te.loc[l,'shape'] = 'larrow'
        te['height'] = np.where(te['shape'].isin(['larrow', 'rarrow']),  height, height/f)

    te['style'] = 'filled' 
    # insert an space between oposite strands.
    add_rows = te[(te.strand > te.strand.shift(1) )& (te.block_id == te.block_id.shift(1))].index.tolist()
    te = gf.insert_na_rows_at_indexes(te, add_rows, fill_rules={'pid':'forward','block_id':'forward'})
    te.loc[te.nucleotide.isna(),['dom','domp','shape','style', 'height', 'width']] = ['white_space', 0, 'rectangle','', height/1.4, 2] 
    te = te.reset_index(drop=True).reset_index()
    color_dict2 = gf.pick_colors(te, column = 'dom', top_domains=top_domains)
    if color_dict == None:
        color_dict = color_dict2
    else:    
        color_dict = {**color_dict2, **color_dict}


    color_dict = {**color_dict,
                  'white_space' :"#FFFFFF",
                  "-" :"#D3D3D3",
                  " " :"#D3D3D3",
                  "?":"#D3D3D3",
                  "TM" :"#D5B60A",
                  "LIPO":"#0000FF",
                  "LP":"#0000FF",
                  "SP":"#FF0000",
                  "SIG":"#FF0000"}
    if yaml_file:
        color_dict = {**color_dict, **yaml_color_dict}

    if light_palette:
         color_dict = {key: value + light_palette for key, value in color_dict.items()}
    
    te['color'] = te.dom.replace(color_dict)
#    te.dom = te.dom.replace({'TM':'', 'LP':'', 'LIPO':'', 'SIG':'', 'SP': ''})
#    te.loc[(te.dom=='') & (te['shape'].isin(['larrow','rarrow'])), 'dom'] = ' '
#    te.loc[(te.dom=='') & (te['shape'].isin(['larrow','rarrow'])), 'color'] = color_dict[" "] 
####### Adding new row to represent empty proteins with TM, SIG or LIPO
##### A lot of code to make sure the order, color and shapes  of the added rows are correct.
    ra = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['rarrow']))].index
    te = gf.insert_na_rows_at_indexes(
       te,
       ra,
       fill_rules={
           'nucleotide':'backward',
           'pid':'backward',
           'block_id':'backward',
           'organism':'backward',
           'arch':'backward',
           'dom':'backward',
           'color':'backward',
           'shape':'rectangle',
           'height':'backward',
           'strand':'backward',
           'style':'backward'})
#    te = te.reset_index(drop=True) 
    ra = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['rarrow']))].index
    te.loc[ra,['shape', 'dom', 'height', 'width']] =[['rarrow', ' ', height * 0.95, TM_width * 3]] # remove color from the colors to change  te.loc[ra,['shape','color', 'dom']] =[['rarrow', color_dict[" "], unknow]]
    la = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['larrow']))].index
    te = gf.insert_na_rows_at_indexes(
       te,
       la,
       insert_position='downstream',
       fill_rules={
           'nucleotide':'forward',
           'pid':'forward',
           'block_id':'forward',
           'organism':'forward',
           'arch':'forward',
           'dom':'forward',
           'color':'forward',
           'shape':'rectangle',
           'height':'forward',
           'strand':'forward',
           'style':'forward'})

#    te = te.reset_index(drop=True) 
    #Replacing the - or ? marker for empty spaces to decrease the noise on the figure.
    la = te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP'])) & (te['shape'].isin(['larrow']))].index
    te.loc[la,['shape', 'dom', 'height', 'width']] =[['larrow', ' ', height * 0.95, TM_width * 3]] # samething as above  remove color from the colors to change  te.loc[ra,['shape','color', 'dom']] =[['rarrow', color_dict[" "], unknow]]
    te['bcolor'] = te.color
    if highlight_domain:
        mask = te.dom.isin(list(highlight_domain.keys()))
        te.loc[mask, 'bcolor'] = te.loc[mask, 'dom'].replace(highlight_domain)

    te['penwidth'] = penwidth
    te['width'] = te.dom.str.len()/40
    te.loc[te.dom =="white_space", 'width'] = space_width
    if architecture:
        te.loc[te.dom.isin(["TM", "SIG", "SP","LIPO","LP"]), ['bcolor', 'width', 'penwidth', 'shape']] =['black', TM_width, 0.5, 'oval']
    else:    
        te.loc[te.dom.isin(["TM", "SIG", "SP","LIPO","LP"]), ['bcolor', 'width', 'penwidth']] =['black', TM_width, 0.5]

    te.loc[te.dom =="PSE", ['color', 'bcolor']] =['#D3D3D340', '#D3D3D3']

    te.loc[(te.dom.isin(['TM','LP','LIPO','SIG','SP', " "])) & (te['shape'].isin(['larrow','rarrow'])), ['dom', 'color', 'bcolor']] = [' ', color_dict[" "], color_dict[" "]]

    if query_asterix:
        to_a = te.query('query ==1').query("dom not in ['TM','LP','LIPO','SIG','SP']").drop_duplicates(subset=['pid'], keep='last').index ###
        te.loc[to_a, 'dom'] = te.loc[to_a, 'dom'] +'*' ###
        te.loc[te.dom == " *", 'dom'] = "*"
    if query_color:    
        to_a2 = te.query('query ==1').query("dom not in ['TM','LP','LIPO','SIG','SP']").index ###
        te.loc[to_a2, 'bcolor'] = query_color ###
    te.dom = te.dom.replace({'TM':'', 'LP':'', 'LIPO':'', 'SIG':'', 'SP': '', 'white_space': ''})
    te.dom = te.dom.replace({'?':unknow, '-':unknow})
    te.loc[te.dom == "", ['height']] = height/1.5
    ####Creating a column with the pid of the query seed to later be used in given a header of the figure.
    te['query_pid'] = te.groupby('block_id')['pid'].transform(lambda x: x[te['query'] == 1].iloc[0] if any(te['query'] == 1) else np.nan)   
    ##### Reseting the index to use the new index as name of each node in the figure
    del te['index']
    te = te.reset_index()

    if sort_by:
        ordered_list = te.sort_values([sort_by,'block_id','blockp']).block_id.drop_duplicates().tolist() 
    else:
        ordered_list = te.block_id.drop_duplicates().tolist() 
    
    te.organism = te.organism.fillna(' ')
    #Getting the mas string size to padd the image headers 
    w1 = te.block_id.str.len().max()
    w2 = te.organism.str.len().max()
    max_width = max(w1, w2)
    if lineage:
        te['lineage'] = te.block_id.map(df[['block_id', 'lineage']].drop_duplicates('block_id').set_index('block_id').lineage.to_dict())
        if genome_coordinates:
            te['to_bellow'] = te.block_id.fillna('-') +'/' + te.organism.fillna('-') + ' (' + te.lineage.fillna('-') + ')'
        else:
            te['to_bellow'] = te.organism.fillna('-') + ' (' + te.lineage.fillna('-') + ')'

        te['to_bellow'] = te.to_bellow.str.replace('>', '&#62;')
    else:    
        if genome_coordinates:
            te['to_bellow'] = te.block_id.fillna('-') +'/' + te.organism.fillna('-')
        else:
            te['to_bellow'] = te.organism.fillna('-')

    ## Adjusting the height os arrows accourding to the string size
    te['dom_string_size'] = te.dom.str.len()
    te.loc[(te.dom_string_size ==1) & (te['shape'].isin(['rarrow', 'larrow'])), 'height'] = height * 0.80
    te.loc[(te.dom_string_size ==2) & (te['shape'].isin(['rarrow', 'larrow'])), 'height'] = height * 0.90
    te.loc[(te.dom_string_size ==3) & (te['shape'].isin(['rarrow', 'larrow'])), 'height'] = height * 0.99
    bellow_max = te.to_bellow.str.len().max()
    gb = te.groupby('block_id')
    li=[]
    label_list=[]
    A = pgv.AGraph()
    if id_position =='bellow':
        bellow = True
    else:
        bellow = False

   # return gb
    for group in ordered_list:
        block_id_df = gb.get_group(group)
        if isinstance( comments, pd.DataFrame):
            if group in comments.index:
                l = []
                node_fake = f'{group}_fake_comment'
                l.append(node_fake)
                A.add_node(node_fake,
                           label='',
                           shape='none',
                           fontsize=fontsize,
                           fontname='Arial',
                           margin=margin,
                           height=height)
                node_name = f'{group}_comments'
                comment_text = comments.loc[group][0].replace("&", "&amp;")
                comment_text = comment_text.replace("\n", "<br/>")
                comment_text = textwrap.wrap(comment_text, width=100)
                comment_text = [line.ljust(100) for line in comment_text]
                comment_text = '<br/>'.join(comment_text)
                
                l.append(node_name)
                A.add_node(node_name,
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{fontsize +2}"><B>{comment_text:<{max_width}}</B></FONT></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Arial',
                           margin=margin,
                           height=height)

                li.append(l)

        if bellow:
            l =[]
            node_fake = f'{group}_fake_text_3'
            l.append(node_fake)
            A.add_node(node_fake,
                       label='',
                       shape='none',
                       fontsize=fontsize,
                       fontname='Arial',
                       margin=margin,
                       height=height)
        else:
            l =[]
            node_fake = f'{group}_fake_text_2'
            l.append(node_fake)
            A.add_node(node_fake,
                       label='',
                       shape='none',
                       fontsize=fontsize,
                       fontname='Arial',
                       margin=margin,
                       height=height)
            l.append(group)
            label_list.append(group)
            if genome_coordinates:
                print('oo')
                A.add_node(group,
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="{fontsize}">{block_id_df.query_pid.unique()[0]}</FONT></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{fontsize}">{block_id_df.block_id.unique()[0]:<{max_width}}</FONT></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas italic" POINT-SIZE="{fontsize}">{block_id_df.organism.unique()[0]:<{max_width}}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Consolas',
                           margin=margin,
                           height=height)
            else:
                print('ooooi')
                A.add_node(group,
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                               f'<TR><TD HEIGHT="4"></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"><FONT FACE="Arial" POINT-SIZE="{fontsize}">{block_id_df.query_pid.unique()[0]}</FONT></TD></TR>'
                              f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas italic" POINT-SIZE="{fontsize}">{block_id_df.organism.unique()[0]:<{max_width}}</FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Consolas',
                           margin=margin,
                           height=height)

        for y,z in block_id_df.iterrows():

            A.add_node(z['index'],
                           label=z.dom,
                           shape=z['shape'],
                           width=z['width'],
                           style=z['style'],
                           height=z['height'],
                           color=z.bcolor,
                           penwidth=z['penwidth'],
                           fillcolor=z.color,
                           labeldistance=labeldistance,
                           margin=margin,
                           fontname="Consolas",
                           fontsize=fontsize)
            l.append(z['index'])
    
        if bellow:
            li.append(l)
            l =[]
            node_fake = f'{group}_fake_text'
            l.append(node_fake)
            A.add_node(node_fake,
                       label='',
                       shape='none',
                       fontsize=fontsize,
                       fontname='Arial',
                       margin=margin,
                       height=height)
            l.append(group)
            label_list.append(group)
            if genome_coordinates:
                A.add_node(group,
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{fontsize}">{block_id_df.query_pid.unique()[0]}/<I>{block_id_df.to_bellow.unique()[0]:<{bellow_max}}</I></FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Consolas',
                           margin=margin,
                           height=height)
            else:
                A.add_node(group,
                           label=(
                               f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
                               f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{fontsize}">{block_id_df.query_pid.unique()[0]}/<I>{block_id_df.to_bellow.unique()[0]:<{bellow_max}}</I></FONT></TD></TR>'
                               '</TABLE>>'),
                           shape='none',
                           fontsize=fontsize,
                           fontname='Consolas',
                           margin=margin,
                           height=height)
            li.append(l)
        else:
            li.append(l)
    # Drawing the nodes using the li list to make the loci in the same rank        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")
    # Creating invisible Edges using the list of first nodes tok align it to the left.    
    #for i in range(len(label_list) - 1):
    #    A.add_edge(label_list[i], label_list[i + 1], style='invis')    


    v = [li[x][0] for x in range(len(li))]
   # [A.add_edge(f"{v.pop(0)}:w", f"{v[0]}:w", penwidth=0) for x in range(len(v)-1)]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.02)
    if to_string:
        buffer = io.BytesIO()
        A.draw(buffer, format='svg', prog='dot')
        svg_output = buffer.getvalue().decode('utf-8')
        buffer.close()
        return svg_output

    
    A.graph_attr["page"] = "8.5,11" 
    A.draw(output_file, prog="dot")
#    with open("Output.txt", "w") as text_file:  In case we ant to ahve the .dot file
#        text_file.write(A.string())
    return te


def compact_to_df(compact,sep='Tab', columns= ['pid', 'compact', 'organism', 'pid_compact']):
    """
    convert TASS annotated compact represantation to DF
    """
    import pandas as pd
    if sep=='Tab':
        x = pd.read_csv(compact, sep="\t", names=columns)
    else:
        x = pd.read_csv(compact, delimiter=r"\s\s+", names=columns)
    x = x.drop_duplicates(subset=['pid'])
    x['arch'] = x.compact.str.split('->')

    x = x.explode('arch')
    x2 = x[['arch','organism']]
    x2['strand'] = -1
    x2.loc[~x2.arch.fillna('').str.contains('<-'), 'strand'] = 1
    x2.arch = x2.arch.str.split('<-')
    x2 = x2.explode('arch')
    x2 = x2.query('arch not in ["", "||"]')
    x2.arch = x2.arch.str.replace('\|\|', '|--')
    x2.arch = x2.arch.str.split('|')
    x2 = x2.explode('arch')
    x2.loc[x2.arch.fillna(" ").str.startswith('--'), 'strand'] = 1
    x2.arch = x2.arch.str.strip('--')
    x2 = x2.reset_index().rename({'index':'block_id'}, axis=1)
    x2['block_id'] = x2.block_id.replace(x.pid.to_dict())
    x2 = x2.reset_index().rename({'index':'pid'}, axis=1)
    x2['query'] = 0
    x2.loc[x2.arch.fillna(' ').str.endswith('*'), 'query'] = 1
    x2.arch = x2.arch.str.strip('*')
    x2['nucleotide'] = 'dummy_collumn' 
    return x2


def check_spacing(text):
    import re
    """
    Check if the string contains a tab or two or more spaces immediately
    after the first word in the input text.

    Parameters:
    -----------
    text : str
        The input string to be checked.

    Returns:
    --------
    str
        A message indicating whether a tab, two or more spaces, or neither 
        were found after the first word:
        - "Tab found" if a tab ('\t') follows the first word.
        - "Two or more spaces found" if two or more spaces follow the first word.
        - "Neither found" if neither is present after the first word.
    """
    # Check for a tab directly after the first word
    if re.match(r'^\S+\t', text):
        return "Tab"
    
    # Check for two or more spaces directly after the first word
    if re.match(r'^\S+  +', text):
        return "Two_spaces"
    
    return "No tabular format found(with spaces/Tab), check your input"
def pid2tax(pidlist,full=False):
    import pandas as pd
    from rotifer.db import ncbi
    ic = ncbi.IPGCursor()
    if isinstance(pidlist,pd.Series):
        pidlist = pidlist.unique().tolist()
    ipgs = ic.fetchall(pidlist)
    # ipgs = pd.concat(ipgs) removing this line because there ar eno reason to have it here, also check duplicate of this function?? what happend with my codes!!!
    assemblies_id = ipgs.assembly.unique().tolist()
    a = ncbi.assemblies(taxonomy=False)
    a.rename({'#assembly_accession': 'assembly_accession'}, axis=1, inplace=True)
    a = a.query('assembly_accession in @assemblies_id')
    tc =   ncbi.TaxonomyCursor()
    taxonomies = tc.fetchall(a.taxid.tolist())
    a = a[['assembly_accession','taxid']]
    result = a.merge(taxonomies,
                     how='left').merge(ipgs[['pid','nucleotide','start','stop', 'assembly']],
                                      left_on='assembly_accession',
                                       right_on='assembly',
                  
                                       how='left'
                                       )
    result['position'] = result.nucleotide + ":" + result.start.astype(str) + '-'  + result.stop.astype(str)                 

    return result


def compact_to_df2(compact,
                   columns= ['pid', 'compact', 'organism'],
                   columns_to_keep=[0,1,2],
                   only_parser=False):
    """
    convert TASS annotated compact represantation to DF
    If only parser, it will retunr the dataframe before transform it
    Your file must have at leat block identification and compact representation colum

    """
    import pandas as pd
    from rotifer.devel.alpha import gian_func as gf

    with open(compact, 'r') as f:
        file = f.readlines()
    file = pd.Series(file, name='lines').to_frame()
    file.lines = file.lines.str.strip('\n')
    file.lines = file.lines.str.strip()
    file.lines = file.lines.str.replace('\t', '  ')

    file = file.query('lines !=""').reset_index(drop=True)
    to_read = file.query('~lines.str.startswith("#")')
    # Auto checking for the type of separator \t or two or more spaces:
    sep=r"\s\s+"

            
    comments_lines = file.query('lines.str.startswith("#")')
    comments_lines.lines = comments_lines.lines.str.strip()

    if comments_lines.empty:
        pass
    else:
        c_ref = gf.find_next_non_existing(comments_lines.index.tolist())
        pid_ref = to_read.loc[c_ref].lines.str.split(sep, expand=True)[0].tolist()
        comments_lines['pid_ref'] = pid_ref
        comments_lines = comments_lines.groupby('pid_ref').agg(l =('lines' ,lambda x : '; '.join(x.tolist())))

    
    x = to_read.lines.str.split(sep, expand=True)
    if only_parser:
        return x
    
    if x.shape[1] > 2:
        x = x[columns_to_keep]
        x.columns = columns
    else:
        x = x[columns_to_keep[:-1]]
        x.columns = columns[:-1]
        x['organism'] = 'Dummy_organism'


    x = x.drop_duplicates(subset=['pid'])
    x['arch'] = x.compact.str.split('->')

    x = x.explode('arch')
    x2 = x[['arch','organism']]
    x2['strand'] = -1
    x2.loc[~x2.arch.fillna('').str.contains('<-'), 'strand'] = 1
    x2.arch = x2.arch.str.split('<-')
    x2 = x2.explode('arch')
    x2 = x2.query('arch not in ["", "||"]')
    x2.arch = x2.arch.str.replace('\|\|', '|--')
    x2.arch = x2.arch.str.split('|')
    x2 = x2.explode('arch')
    x2.loc[x2.arch.fillna(" ").str.startswith('--'), 'strand'] = 1
    x2.arch = x2.arch.str.strip('--')
    x2 = x2.reset_index().rename({'index':'block_id'}, axis=1)
    x2['block_id'] = x2.block_id.replace(x.pid.to_dict())
    x2 = x2.reset_index().rename({'index':'pid'}, axis=1)
    x2['query'] = 0
    x2.loc[x2.arch.fillna(' ').str.endswith('*'), 'query'] = 1
    x2.arch = x2.arch.str.strip('*')
    x2['nucleotide'] = 'dummy_collumn' 
    return x2, comments_lines

def find_next_non_existing(lst):
    result = []
    seen = set(lst)  # Store the original elements in a set for quick lookup
    
    for num in lst:
        next_num = num + 1  # Start with the next integer
        # Increment until we find a unique number not in the original list
        while next_num in seen:
            next_num += 1
        result.append(next_num)
    
    return result


def get_relations(domdf,
                  dom='full',
                  min_connections = 5,
                  self_relation = False,
                  to_yaml=False,
                  iteration=1,
                  filter_rename_yaml  = False):
    from rotifer.devel.alpha import gian_func as gf
    """ 
    domdf to network, it can create dictionary to create a YAML file or a dataframe.
    dom: list of domains to be used as query in the network, if full all domdf will be used
    min_connections: Minimum number of relationship (Operon or fusion) anode should have to be presente in the network.
    yaml: If decided to transform the network in dictionary to latter be exported as YAML file.
    Iteration: If a list of domains was supplied, the number of iteration in the network to collect other nodes.
    """
    domdf.dom = domdf.dom.str.strip()
    if filter_rename_yaml:
        import yaml
        with open(filter_rename_yaml, "r") as file:
            data = yaml.safe_load(file)
            data = pd.DataFrame(data).T
            data.loc[data.Display_name =="", 'Display_name'] = data.loc[data.Display_name ==""].index.tolist()
            rename_dict = data.Display_name.to_dict()
            domdf.dom = domdf.dom.replace(rename_dict)
            to_display = data.query('Include ==1').Display_name.unique().tolist()
            domdf = domdf.query('dom in @to_display')

    if isinstance(dom, str):
        dom = [dom]
    if dom[0] =='full':
        block_id = domdf.block_id.tolist()
    else:
        for x in range(iteration):
            block_id = domdf.query('dom in @dom').block_id.tolist()
            dom = domdf.query('block_id in @block_id').dom.unique().tolist()
    to_net = domdf.query('block_id in @block_id')[['pid', 'block_id', 'dom', 'strand']]
    to_net['block_pos'] = 1
    to_net['block_pos'] = to_net.groupby(['block_id'])['block_pos'].cumsum()
    to_net = to_net.merge(to_net, on ='block_id')
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Downstream'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Upstream'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Nterminal'
    to_net.loc[((to_net.strand_x == 1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Cterminal'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Downstream'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x != to_net.pid_y)), 'edge_type'] = 'Upstream'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x > to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Cterminal'
    to_net.loc[((to_net.strand_x == -1 ) & (to_net.block_pos_x < to_net.block_pos_y) & (to_net.pid_x == to_net.pid_y)), 'edge_type'] = 'Nterminal'

    to_net= to_net.query('edge_type.notna()')
    to_net = to_net[['block_id','dom_x','dom_y','edge_type']]
    to_net.columns = ['block_id','source','target','edge_type']
    if to_yaml:
        to_net['edge_count'] = to_net.groupby(['source', 'target']).edge_type.transform('count')
        to_net = to_net.query('edge_count >= @ min_connections').drop(['edge_count'], axis=1)
        yamldf = pd.pivot_table(to_net,
                            index='source',
                            columns='edge_type',
                            values = ['target', 'block_id'],
                            aggfunc={
                                'target':lambda x : gf.count_series(x, as_list=True),
                                'block_id': 'nunique'}
                            )
        yamldf = yamldf.loc[:,'block_id'].sum(axis=1).rename('total').to_frame().join(yamldf.loc[:,'target']).sort_values('total')
        yamldf.fillna('', inplace=True)
        yamldf.Nterminal = yamldf.Nterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Cterminal = yamldf.Cterminal.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Downstream = yamldf.Downstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf.Upstream = yamldf.Upstream.apply(lambda x: x.split(' ') if isinstance(x, str) else x)
        yamldf['L'] = yamldf[['Cterminal', 'Downstream', 'Nterminal', 'Upstream']].apply(
                lambda row: [item for sublist in row for item in sublist if item != ''],
                axis=1
                )
        yamldf.L = yamldf.L.apply(len)
        yamldf = yamldf.query('L >0')
        yamldf = yamldf.drop(['L'], axis=1)

        return(yamldf.to_dict(orient='index'))
    else:
        to_net['edge_total'] = to_net.groupby(['source', 'target']).edge_type.transform('count')
        to_net = to_net.groupby(['source', 'target', 'edge_type']).agg(
                {'block_id' :'count',
                 'edge_total' : 'first'}).reset_index().sort_values('block_id')
        to_net = to_net.rename({'block_id':'edge_count'},axis=1)        
        to_net = to_net.query('edge_total >= @min_connections').drop(['edge_total'], axis=1)
        if self_relation:
            pass
        else:
            to_net = to_net[to_net.source != to_net.target]
        return(to_net)

def yaml2net(dict_from_yaml):
    """
    Transform the dictionary from a yaml network file back to the network
    """
    dfc  = pd.DataFrame(dict_from_yaml).T
    dfc.index.name = 'source'
    a = dfc.Cterminal.explode().dropna().str.extract(r'(.+?)(?:\W*\((\d+)\))?$').rename({0: 'target', 1:'cterminal'}, axis=1)
    b = dfc.Downstream.explode().dropna().str.extract(r'(.+?)(?:\W*\((\d+)\))?$').rename({0: 'target', 1:'dowstream'}, axis=1)
    c = dfc.Nterminal.explode().dropna().str.extract(r'(.+?)(?:\W*\((\d+)\))?$').rename({0: 'target', 1:'Nterminal'}, axis=1)
    d = dfc.Upstream.explode().dropna().str.extract(r'(.+?)(?:\W*\((\d+)\))?$').rename({0: 'target', 1:'Upstream'}, axis=1)
    te = pd.concat([a,b,c,d]).reset_index().melt(id_vars=['source', 'target']).dropna()
    te.columns = ['source', 'target', 'edge_type', 'edge_count']
    te.edge_count = te.edge_count.astype(int)
    return te.sort_values('edge_count')


def plot_network(networkdf,
                 community_to_color='leidein',
                 outputfile='net.svg',
                 view=False,
                 node_size='freq',
                 color_edge=True,
                 highlight_query=[],
                 position='spring'):
    """
    PLot the network in svg
    """
    import community
    import seaborn as sns
    import leidenalg as la
    import igraph as ig
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    import networkx as nx

    G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count'])
    #Community detection by Louvain
    partition = community.best_partition(G,weight='edge_count')
    #Df to easily map community to node:
    c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
       {'index': 'cluster', 0: 'Louvain'}, axis=1)
    #Leiden partition:
    part = la.find_partition(ig.Graph.from_networkx(G), la.ModularityVertexPartition, weights='edge_count')
    c['leidein'] = pd.Series(part.membership)
    #Leiden partition playing with resolution parameter:
    part2 = la.find_partition(ig.Graph.from_networkx(G), la.CPMVertexPartition, weights='edge_count', resolution_parameter = 2.5)
    c['leidein_2'] = pd.Series(part2.membership)


    #Adding a color column for each community method used:
    c['Louvain_color'] = c.Louvain.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_color'] = c.leidein.map(pd.Series(sns.color_palette('pastel',c.leidein.nunique()).as_hex()).to_dict())
    c['leidein_2_color'] = c.leidein_2.map(pd.Series(sns.color_palette('pastel',c.leidein_2.nunique()).as_hex()).to_dict())

    #Creating a list to create relative node size based in frequency
    if node_size=='freq':
        Node_size_dict =  (np.sqrt(pd.concat([networkdf.source, networkdf.target]).value_counts()) *30).to_dict()
        node_size = [Node_size_dict[x] for x in G.nodes]
    else:
        node_size = node_size

    # getting specific colors and for edges:
    d = c.set_index(community_to_color).fillna('black')[f'{community_to_color}_color'].to_dict()

    #Colors for Nodes
    nc = c.set_index('cluster')[f'{community_to_color}_color'].fillna('black').to_dict()
    Node_colors = [nc[x] for x in G.nodes]

    Node_edge = ['red' if x in highlight_query else nc [x] for x in G.nodes ]

    #Edge Colors
    if color_edge:
        community_dict = c.set_index('cluster')[community_to_color].to_dict()
        edge_colors=[]
        ew =[]
        for x,y in G.edges:
            if community_dict[x] == community_dict[y]:
                edge_colors.append(d[community_dict[x]])
            else:
                edge_colors.append('black')


    #Drawing the network:
    fig, ax = plt.subplots(figsize=(20,14 ))
    if position=='spring':
        pos = nx.spring_layout(G)
    else:    
        pos = nx.kamada_kawai_layout(G)
    nx.draw_networkx_nodes(G, pos, c.cluster, node_shape="o",node_size=node_size, node_color=Node_colors, edgecolors=Node_edge)
    nx.draw_networkx_edges(G, pos, alpha=0.5, edge_color=edge_colors)
    label_pos = {node: (x, y - 0.03) for node, (x, y) in pos.items()}
    nx.draw_networkx_labels(G,label_pos, font_size= 8)
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if view==True:
        plt.plot()
    else:
        plt.savefig(outputfile, format='svg', bbox_inches='tight')

def domain2compact(domdf, domain_rename=False):
    if domain_rename:
        if isinstance(domain_rename, dict):
            domdf.dom = domdf.dom.replace(domain_rename)
        else:
           import yaml
           with open(domain_rename, "r") as file:
               data = yaml.safe_load(file)
               s = pd.DataFrame(data).T
               tex = s.query('Include ==0').index.tolist()
               s.loc[tex, 'Display_name'] = '?'
               domain_rename = s.Display_name.where(lambda x : x!="").dropna().to_dict()
               domdf.dom = domdf.dom.replace(domain_rename)

    g = domdf.copy().groupby(['pid'], sort=False).agg(
            block_id = ('block_id', 'first'),
            organism =  ('organism' ,'first'),
           strand = ('strand', 'first'),
            query = ('query', 'first'),
            arch = ('dom', lambda x: '+'.join(list(x))))
    if domain_rename:
        ##Tidyng multiples ? from the excluded node list (removing tralling or ? from proteins with annotations)
        g['arch'] = g['arch'].str.replace(r'^\?\+|\+\?', '', regex=True)

    g.loc[g['query'] ==1, 'arch'] = g.loc[g['query'] ==1,'arch'] +'*'
    g.loc[g['strand'] ==1, 'arch'] = g.loc[g['strand'] ==1,'arch'] + '->'
    g.loc[g['strand'] == -1, 'arch'] = '<-' + g.loc[g['strand'] == -1, 'arch'].astype(str)
    m1  = (g.strand ==  1) & (g.block_id == g.block_id.shift(-1)) & (g.strand != g.strand.shift(-1))
    m2  = (g.strand == -1) & (g.block_id == g.block_id.shift(-1)) & (g.strand != g.strand.shift(-1))
    m3 =  (g.strand ==  1) & (g.block_id == g.block_id.shift(1)) & (g.strand != g.strand.shift(1))
    m4 =  (g.strand == -1) & (g.block_id == g.block_id.shift(1)) & (g.strand != g.strand.shift(1))
    g.loc[m1, 'arch'] = g.loc[m1,'arch'].astype(str) + '|'
    g.loc[m2, 'arch'] = g.loc[m2,'arch'].astype(str) + '|'
    g.loc[m3, 'arch'] = '|' + g.loc[m3,'arch'].astype(str)
    g.loc[m4, 'arch'] = '|' + g.loc[m4,'arch'].astype(str)
    g = g.groupby(['block_id'], sort=False).agg(
            compact = ('arch', lambda x: ''.join(list(x))),
            organism =  ('organism' ,'first'))
    return g

def domtable_to_yaml_names(domtable, cutoff=0, to_file=False, color=False):
    def safe_for_yaml(x):
        try:
            # First handle np.nan, None
            if pd.isna(x):
                return ''
            # Now handle string "nan" (case insensitive)
            if isinstance(x, (str, np.str_)):
                if x.strip().lower() == 'nan':
                    return ''
                else:
                    return str(x)
            # Safe types
            elif isinstance(x, (int, float, bool)):
                return x
            else:
                # fallback: convert to string
                return str(x)
        except Exception:
            # In case x is weird object or array, fallback
            return ''


    to_remove = ["TM", "SP", "LP", "LIPO", '-', '?', 'SIG', 'tRNA', 'PSE']
    to_remove = domtable.query('dom in @to_remove').dom.unique().tolist()
    to_dict = domtable.dom.value_counts().to_frame()
    if pd.__version__.startswith('2.'):
        to_dict.rename({'count': 'dom'}, axis=1, inplace=True) # Add to finx compatibility
    to_dict.loc[to_dict.dom > cutoff, 'Display_name'] = to_dict.loc[to_dict.dom > cutoff].index.tolist()
    to_dict['Include']= 0
    to_dict.loc[to_dict.dom > cutoff, 'Include'] = 1
    to_dict.loc[to_remove , 'Display_name'] = '-'
    to_dict['Display_name'] = to_dict['Display_name'].fillna('-')
    to_dict.loc[to_dict.Display_name =='', 'Include'] = 0
    to_dict['Function'] = [['', ''] for x in range(to_dict.shape[0])]
    to_dict['Notes'] = [['', ''] for x in range(to_dict.shape[0])]
    if color:
        if isinstance(color,dict):
            to_dict['color'] = to_dict.index.map(color).fillna("#D3D3D3")
            to_dict.loc[to_dict.dom <= cutoff, 'color'] = '#D3D3D3'
        else:
            print('color must be a domain color dict (gf.pickcolor on domtable)') 

    to_dict = to_dict.drop('dom', axis=1)
    to_dict = to_dict.applymap(safe_for_yaml) 
    to_dict = to_dict.to_dict(orient='index')
    if to_file:
        import yaml
        with open(to_file, 'w') as file:
            yaml.dump(to_dict, file, default_flow_style=False, sort_keys=False)
    else:
        return to_dict
def yamldomain2df(file):
    """
    Load yaml file created by the domtable_to_yaml_names into a pandas dataframe
    """
    import yaml
    with open(file, "r") as file:
        data = yaml.safe_load(file)
        data = pd.DataFrame(data).T
        data.Display_name =  data.Display_name.fillna('')
     
    return data 


def display_html_popup_from_file(file_path, title="Popout Window", use_button=True):
    from IPython.display import display, HTML
    import html
    # Read the content of the HTML or SVG file
    with open(file_path, 'r', encoding='utf-8') as file:
        html_content = file.read()
    
    # Properly escape the content for safe inclusion in the HTML
    escaped_content = html.escape(html_content).replace("\n", "&#10;")

    if use_button:
        # Generate the pop-out script for a new window
        script = f"""
        <script>
            function openPopup() {{
                var newWindow = window.open("", "{title}", "width=800,height=600");
                newWindow.document.write(`
                    <html>
                        <head>
                            <title>{title}</title>
                        </head>
                        <body>
                            {html_content}
                        </body>
                    </html>
                `);
                newWindow.document.close();
            }}
        </script>
        <button onclick="openPopup()">Open Content in Popout</button>
        """
    else:
        # Directly display the content embedded in an iframe
        script = f"""
        <iframe style="border:none; width:100%; height:600px;" srcdoc="{escaped_content}"></iframe>
        """

    display(HTML(script))

def is_running_in_jupyter():
    try:
        from IPython import get_ipython
        if 'IPKernelApp' in get_ipython().config:
            return True
        else:
            return False
    except (ImportError, AttributeError):
        return False    


def find_potential_typos(df, column, threshold=2):
    import textdistance
    results = []
    for i, word1 in enumerate(df[column]):
        for j, word2 in enumerate(df[column]):
            if i < j:  # To avoid duplicate and self-comparisons
                distance = textdistance.levenshtein(word1, word2)
                if distance <= threshold:  # Change threshold as needed
                    results.append((word1, word2, distance))
    return results


def plot_network2(networkdf,
                 community_to_color='leidein',
                 outputfile='net.svg',
                 view=False,
                 node_size='freq',
                 color_edge=True,
                 query_list=[],
                 position='Kamada-kawai',
                 legend=False,
                 legend_title='Functional theme', 
                 highlight_color='green',
                 highlight_list =[], 
                 return_position=False,
                 label=True,
                 font_size=4,
                 fig_size = (11.7,8.3)):
    """
    PLot the network in svg
    fig_size 11.7, 8.3 is a full size landscape in a standard journal
    """
    import community
    import seaborn as sns
    import leidenalg as la
    import igraph as ig
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    from rotifer.devel.alpha import gian_func as gf
    import networkx as nx
    plt.rcParams['svg.fonttype'] = 'none'

    G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count'])
    #Community detection by Louvain
    partition = community.best_partition(G,weight='edge_count')
    #Df to easily map community to node:
    c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
       {'index': 'cluster', 0: 'Louvain'}, axis=1)
    #Leiden partition:
    part = la.find_partition(ig.Graph.from_networkx(G), la.ModularityVertexPartition, weights='edge_count')
    c['leidein'] = pd.Series(part.membership)
    #Leiden partition playing with resolution parameter:
    part2 = la.find_partition(ig.Graph.from_networkx(G), la.CPMVertexPartition, weights='edge_count', resolution_parameter = 2.5)
    c['leidein_2'] = pd.Series(part2.membership)


    #Adding a color column for each community method used:
    c['Louvain_color'] = c.Louvain.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_color'] = c.leidein.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_2_color'] = c.leidein_2.map(pd.Series(sns.color_palette('pastel',c.leidein_2.nunique()).as_hex()).to_dict())


    #Colors for Nodes
    if isinstance(community_to_color, str):
        nc = gf.get_network_community_color(networkdf, community_to_color = community_to_color)
    else:
        nc = community_to_color

    Node_colors = [nc[x] for x in G.nodes]

    #Edge Colors
    if color_edge:
        edge_colors=[]
        edge_alpha=[]
        edge_styles=[]
        for x,y in G.edges:
            if nc[x] == nc[y]:
                edge_colors.append(nc[x])
                edge_styles.append('-')
                edge_alpha.append(0.7)
            else:
                edge_colors.append('black')
                edge_styles.append(':')
                edge_alpha.append(0.2)

    # Creating the G fot the queries, and the G tor the other nodes:
    G_queries = G.copy()
    G_queries = G_queries.subgraph(query_list)
    Node_query_colors = [nc[x] for x in G_queries.nodes]
    G_others = G.copy()
    G_others.remove_nodes_from(query_list)
    Node_others_colors = [nc[x] for x in G_others.nodes]

    #Creating a list to create relative node size based in frequency
    if node_size=='freq':
        Node_size_dict =  (np.sqrt(pd.concat([networkdf.source, networkdf.target]).value_counts()) *30).to_dict()
        node_size = [Node_size_dict[x] for x in G.nodes if x not in query_list]
        node_query_size = [Node_size_dict[x] for x in G_queries.nodes]
        node_others_size = [Node_size_dict[x] for x in G_others.nodes]
    else:
        node_size = node_size

    #Drawing the network:
    fig, ax = plt.subplots(figsize=(fig_size))
    if isinstance(position, dict):
        pos = position
        print ('Dict position')
    else:
        if position=='spring':
            pos = nx.spring_layout(G)
            print('Spring')
        else:    
            pos = nx.kamada_kawai_layout(G)
            print('Kamada')
    if len(query_list) > 0 :
        nx.draw_networkx_nodes(G_queries,
                               pos,
                               node_shape="h",
                               node_size=node_query_size,
                               node_color=Node_query_colors,
                               edgecolors=[highlight_color if node in highlight_list else nc[node] for node in G_queries.nodes])
        nx.draw_networkx_nodes(G_others,
                               pos,
                               node_shape="o",
                               node_size=node_others_size,
                               node_color=Node_others_colors,
                               edgecolors=[highlight_color if node in highlight_list else nc[node] for node in G_others.nodes])
    else:
        nx.draw_networkx_nodes(G,
                               pos,
                               node_shape="o",
                               node_size=node_size,
                               node_color=Node_colors,
                               edgecolors=[highlight_color if node in highlight_list else nc[node] for node in G.nodes])
        
    nx.draw_networkx_edges(G, pos, alpha=edge_alpha, edge_color=edge_colors, style=edge_styles)
    label_pos = {node: (x, y - 0.03) for node, (x, y) in pos.items()}
    if label:
        nx.draw_networkx_labels(G,label_pos, font_size= font_size)
    plt.tight_layout()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    if legend:
        import matplotlib.patches as mpatches
        legend_handles = [mpatches.Patch(color=color, label=group) for group, color in legend.items()]
        ax.legend(handles=legend_handles, title=legend_title, fontsize=font_size)
    if view==True:
        if return_position:
            plt.plot()
            return(pos)
        plt.plot()
    else:
        if return_position:
            plt.savefig(outputfile, format='svg', bbox_inches='tight')
            return(pos)

        plt.savefig(outputfile, format='svg', bbox_inches='tight')

def get_network_community_color(networkdf,
                 community_to_color='leidein'):
    """
    PLot the network in svg
    """
    import community
    import seaborn as sns
    import leidenalg as la
    import igraph as ig
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    import networkx as nx

    G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count'])
    #Community detection by Louvain
    partition = community.best_partition(G,weight='edge_count')
    #Df to easily map community to node:
    c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
       {'index': 'cluster', 0: 'Louvain'}, axis=1)
    #Leiden partition:
    part = la.find_partition(ig.Graph.from_networkx(G), la.ModularityVertexPartition, weights='edge_count')
    c['leidein'] = pd.Series(part.membership)
    #Leiden partition playing with resolution parameter:
    part2 = la.find_partition(ig.Graph.from_networkx(G), la.CPMVertexPartition, weights='edge_count', resolution_parameter = 2.5)
    c['leidein_2'] = pd.Series(part2.membership)


    #Adding a color column for each community method used:
    c['Louvain_color'] = c.Louvain.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_color'] = c.leidein.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_2_color'] = c.leidein_2.map(pd.Series(sns.color_palette('pastel',c.leidein_2.nunique()).as_hex()).to_dict())


    # getting specific colors and for edges:
    d = c.set_index(community_to_color).fillna('black')[f'{community_to_color}_color'].to_dict()

    #Colors for Nodes
    nc = c.set_index('cluster')[f'{community_to_color}_color'].fillna('black').to_dict()

    return nc

def reposition_nodes_keep(G, pos, x_correction=0, y_correction=0):
    """
    Keeps the specified nodes in the graph, removes all other nodes.
    The remaining nodes are repositioned to maintain the same aspect ratio and fill the space.

    Parameters:
    G (networkx.Graph): The graph from which nodes will be kept.
    nodes_to_keep (list): A list of nodes to keep in the graph.

    Returns:
    dict: A dictionary of positions for the remaining nodes after scaling.
    """

    # Get the positions of the remaining nodes
    remaining_nodes = list(G.nodes)
    remaining_pos = {node: pos[node] for node in remaining_nodes}

    # Calculate the bounding box for the remaining nodes
    x_vals = [p[0] for p in remaining_pos.values()]
    y_vals = [p[1] for p in remaining_pos.values()]

    # Find the scaling factors based on the bounding box of the remaining nodes
    x_range = max(x_vals) - min(x_vals)
    y_range = max(y_vals) - min(y_vals)

    # Avoid division by zero by defaulting scale factors to 1
    scale_x = 1 / x_range if x_range != 0 else 1
    scale_y = 1 / y_range if y_range != 0 else 1

    # Scale positions to maintain aspect ratio
    scaled_pos = {node: ((pos[0] * scale_x) + x_correction, (pos[1] * scale_y) + y_correction) for node, pos in remaining_pos.items()}

    # Return the scaled positions
    return scaled_pos

def remote_blast(acc,
             db='nr_cluster_seq',
             aln=False,
             max_out = 5000):
    '''
    Remote_blast accepts sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    cols =[
            'hit',
            'query',
            'hitstart',
            'hitend',
            'evalue',
            'querycoverage',
            'querystart',
            'queryend',
            'iteration',
            'bitscore',
            'length'
            ]
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            acc.to_file(f'{tmpdirname}/seqfile') 
            if aln:
                Popen(f'psiblast -in_msa {tmpdirname}/seqfile -db {db} -num_alignments {max_out} -num_descriptions {max_out}  > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            else:
                Popen(f'blastp -query {tmpdirname}/seqfile -db {db} -num_alignments {max_out} -num_descriptions {max_out}  > {tmpdirname}/out',
                      stdout=PIPE,
                      shell=True
                      ).communicate()
            #Popen(f'blast2table {tmpdirname}/out > {tmpdirname}/out.tsv',
            #      stdout=PIPE,
            #      shell=True).communicate()
            #t = pd.read_csv(f'{tmpdirname}/out.tsv', sep="\t", names=cols)
            with open(f'{tmpdirname}/out') as f:
                    blast_r = f.read()
        
    return  blast_r 


def dash_aligner_view(seqobj):
    import dash_bio as dashbio
    from dash import Dash, html, Input, Output, callback
    import io
    # Convert the FASTA string to a format compatible with Dash Bio
    fasta_string = seqobj.to_string()
    fasta_data = io.StringIO(fasta_string)

    app = Dash(__name__)

    app.layout = html.Div([
        dashbio.AlignmentChart(
            id='my-default-alignment-viewer',
            data=fasta_data.read(),
            height=900,
            tilewidth=30,
        ),
        html.Div(id='default-alignment-viewer-output')
    ])

    @callback(
        Output('default-alignment-viewer-output', 'children'),
        Input('my-default-alignment-viewer', 'eventDatum')
    )
    def update_output(value):
        if value is None:
            return 'No data.'
        return str(value)

    if __name__ == '__main__':
        app.run(debug=True)
        
def wrap2html(match, background='black'):
    from rotifer.devel.alpha import gian_func as gf
    import re
    if background =='black':
        foreground = '#FFFFFF'
    else:
        foreground = 'black'
    
    first_dict = gf.colors_for_html_function(foreground=foreground)[0]

    def wrap_match(match):
        return f'<SPAN class={first_dict[match.group(0)[0]]}>{match.group(0)}</SPAN>'
    #gconsensus = {'h':r'[AILMFWVh]', 'p':r"[KR\+]", 'n' :r'[ED\-]', 'po' : r"[NQST]", "C":r"C", "pro":r"P", "G":r"G", "Ali":r"[HY]", '-':r'[\-xX\.]'}
    gp = {'h':r'[AILMFWVh]', 'p':r"[KR]", 'n' :r'[ED]', 'po' : r"[NQST]", "C":r"C", "pro":r"P", "G":r"G", "Ali":r"[HY]", '-':r'[\-xX]'}
    combined_pattern = '|'.join(f"({pattern}+)" for pattern in gp.values())
    result = re.sub(combined_pattern, wrap_match, match)
    return result

def to_html2(seqobj, output, background = 'black', columns =[], consensus=(50,60,70,80,90)):
    from rotifer.devel.alpha import gian_func as gf
    if background =='black':
        foreground = '#FFFFFF'
    else:
        foreground = 'black'
   ### creating the CSS style block with  the color dict 
    cd = gf.colors_for_html_function(foreground=foreground)[1]
    css_rules = ""
    for key, value in cd.items():
        css_rule = f"    .{key} {{ color: {value}; }}\n"
        css_rules += css_rule
    # Wrap CSS rules in <style> tag
    css_style_block = f"<style>\n{css_rules}</style>"

    s = seqobj.copy()
    #Getting the length of the aligment:
    l = len(s.df.iloc[0,1])
    s = s.add_consensus(cutoffs=consensus)
    s.df = pd.concat([s._scale_bar(l), s.df])

    s.df['html'] = s.df.sequence.apply(gf.wrap2html)
    if len(columns) <1:
        columns = ['id','html']
    else:
        columns = columns.append(['id','html'] )
    my_text = gf.padding_df(s.df[columns]).to_csv(sep="\t", index=None, header = None)
    """
    Save the MSA as HTML file
    
    Parameters:
    - file_name: str, name of the file to save (e.g., "output.html").
    - Background color: str, [white or black].
    """
    html_content = f"""<!DOCTYPE html>
<HTML>
<HEAD>
<META http-equiv="Content-Type" content="text/html; charset=utf-8"/>
<TITLE>Rotifer HTML view</TITLE>
{css_style_block}

</HEAD>
<BODY style="background-color:black;color:white">
<TABLE style="border:0px; background-color:black; color:white; a:link:blue; a:active:red; a:visited:purple;">
<TR><TD><PRE>
</PRE></TD></TR>
<TR><TD><PRE>
{my_text}
</PRE></TD></TR>
</TABLE>
<P><SMALL><A HREF="my page/">rotifer</A> 2024 <A HREF="mailto:ggnicastr@gmail.com">Gianlucca G. Nicastro</A></SMALL></P>
</BODY>
</HTML>
"""
    # Write the HTML content to the file
    with open(output, "w", encoding="utf-8") as file:
        file.write(html_content)
    print(f"HTML file '{output}' has been successfully created.")

def yaml2df(file, transpose=True):
    from rotifer.core.functions import loadConfig
    import pandas as pd
    yaml_dict = loadConfig(file)
    df = pd.DataFrame(yaml_dict)
    if transpose:
        return df.T
    else:
        return df

def yamldf2yamlfile(yamldf, output_file):
    import yaml
    """
    Function to convert yaml table to yaml file
    """
    to_file = yamldf.to_dict(orient='index')
    with open(output_file, 'w') as file:
        yaml.dump(to_file, file, default_flow_style=False, sort_keys=False)
    return print(f'dataframe saved as {output_file}')    


def plot_network_dash(networkdf,
                 community_to_color='leidein',
                 outputfile='net.svg',
                 view=False,
                 node_size='freq',
                 color_edge=True,
                 query_list=[],
                 position='Kamada-kawai',
                 legend=False,
                 legend_title='Functional theme', 
                 highlight_color='green',
                 highlight_list =[], 
                 return_position=False,
                 label=True,
                 font_size=4,
                 fig_size = (11.7,8.3)):
    """
    PLot the network in svg
    fig_size 11.7, 8.3 is a full size landscape in a standard journal
    """
    import community
    import seaborn as sns
    import leidenalg as la
    import igraph as ig
    import matplotlib.pyplot as plt
    from matplotlib import cm as cm
    import numpy as np
    from rotifer.devel.alpha import gian_func as gf
    import networkx as nx
    plt.rcParams['svg.fonttype'] = 'none'

    G = nx.from_pandas_edgelist(networkdf, edge_attr=['edge_type', 'edge_count'])
    #Community detection by Louvain
    partition = community.best_partition(G,weight='edge_count')
    #Df to easily map community to node:
    c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
       {'index': 'cluster', 0: 'Louvain'}, axis=1)
    #Leiden partition:
    part = la.find_partition(ig.Graph.from_networkx(G), la.ModularityVertexPartition, weights='edge_count')
    c['leidein'] = pd.Series(part.membership)
    #Leiden partition playing with resolution parameter:
    part2 = la.find_partition(ig.Graph.from_networkx(G), la.CPMVertexPartition, weights='edge_count', resolution_parameter = 2.5)
    c['leidein_2'] = pd.Series(part2.membership)


    #Adding a color column for each community method used:
    c['Louvain_color'] = c.Louvain.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_color'] = c.leidein.map(pd.Series(sns.color_palette('pastel',c.Louvain.nunique()).as_hex()).to_dict())
    c['leidein_2_color'] = c.leidein_2.map(pd.Series(sns.color_palette('pastel',c.leidein_2.nunique()).as_hex()).to_dict())


    #Colors for Nodes
    if isinstance(community_to_color, str):
        nc = gf.get_network_community_color(networkdf, community_to_color = community_to_color)
    else:
        nc = community_to_color

    Node_colors = [nc[x] for x in G.nodes]

    #Edge Colors
    if color_edge:
        edge_colors=[]
        edge_alpha=[]
        edge_styles=[]
        for x,y in G.edges:
            if nc[x] == nc[y]:
                edge_colors.append(nc[x])
                edge_styles.append('-')
                edge_alpha.append(0.7)
            else:
                edge_colors.append('gray')
                edge_styles.append(':')
                edge_alpha.append(0.2)

    # Creating the G fot the queries, and the G tor the other nodes:
    G_queries = G.copy()
    G_queries = G_queries.subgraph(query_list)
    Node_query_colors = [nc[x] for x in G_queries.nodes]
    G_others = G.copy()
    G_others.remove_nodes_from(query_list)
    Node_others_colors = [nc[x] for x in G_others.nodes]
    pos = nx.kamada_kawai_layout(G)


    #Creating a list to create relative node size based in frequency
    if node_size=='freq':
        Node_size_dict =  (np.sqrt(pd.concat([networkdf.source, networkdf.target]).value_counts()) *3).to_dict()
        node_size = [Node_size_dict[x] for x in G.nodes]
        node_query_size = [Node_size_dict[x] for x in G_queries.nodes]
        node_others_size = [Node_size_dict[x] for x in G_others.nodes]
    else:
        node_size = node_size
    
    return (G,Node_colors,node_size, nc, edge_colors,edge_styles, edge_alpha, pos)


def get_pubmed_bibtex(pmid):
    from Bio import Entrez
    from rotifer.db.ncbi import NcbiConfig as config

    Entrez.email = config['email']  # Provide your email for Entrez requests
    try:
        # Fetch the PubMed article information in XML format
        handle = Entrez.esummary(db="pubmed", id=pmid, retmode="xml")
        records = Entrez.read(handle)
        if not records:
            print(f"Error: No records found for PMID {pmid}")
            return None

        docsum = records[0]  # Get the first record

        # Extract information (safe checking for missing elements)
        title = docsum.get("Title", "No title available")
        authors = docsum.get("AuthorList", ["No authors available"])
        source = docsum.get("Source", "No source available")
        pubdate = docsum.get("PubDate", "No publication date available")
        pmid_value = docsum.get("Id", "No PMID available")

        # Create BibTeX entry
        bibtex_entry = f"@article{{{pmid_value},\n"
        bibtex_entry += f"  author = {{{', '.join(authors)}}},\n"
        bibtex_entry += f"  title = {{{title}}},\n"
        bibtex_entry += f"  journal = {{{source}}},\n"
        bibtex_entry += f"  year = {{{pubdate}}},\n"
        bibtex_entry += f"  publisher = {{{source}}},\n"
        bibtex_entry += f"  PMID = {{{pmid_value}}}\n"
        bibtex_entry += f"}}\n"

        return bibtex_entry
    except Exception as e:
        print(f"Error fetching data for PMID {pmid}: {e}")
        return None

def generate_bibtex(pmid_list, filename):
    from rotifer.devel.alpha import gian_func as gf
    with open(filename, "w") as bibfile:
        for pmid in pmid_list:
            # Fetch the BibTeX entry for each PMID
            bibtex_entry = gf.get_pubmed_bibtex(pmid)
            if bibtex_entry:
                bibfile.write(bibtex_entry)


def extract_citations_from_latex(latex_file):
    import re
    with open(latex_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    citations = re.findall(r'\\cite\{(.*?)\}', content)
    
    # Flatten list if there are multiple citations inside a single \cite{a,b,c}
    citation_list = [cite.strip() for group in citations for cite in group.split(',')]
    
    return set(citation_list)  # Use set to remove duplicates

def create_models_dict_db(suffix='pkl'):
    """ 
    Load multiple seqobject pickle files in a dictionary
    Is importante to have the .name agrument from qhe seqobject as a proper name and not the name of the first sequence
    """

    from glob import glob
    files = glob(f'*.{suffix}')
    alndb = {}
    for x in files:
        m = pd.read_pickle(x)
        alndb[m.name] = m

    return alndb    

def create_HHMdb(model_dict_db, hhm_db_name='hmmdb'):
    """
    Using the model_dict db to create an HHM database, it uses the names of the models that cames from seqobjc.name
    """

    import tempfile
    import subprocess
    import rotifer
    from rotifer.devel.alpha import gian_func as gf
    import os
    import shutil

    hhdb_script  = f"{rotifer.config['base']}/bin/build_hhmDB.sh"
    alndb = model_dict_db
    with tempfile.TemporaryDirectory() as temp_dir:
        os.makedirs(f'{temp_dir}/alns', exist_ok=True)
        for x in alndb.keys():
            aln_file = os.path.join(temp_dir, 'alns', f'{x}.fa')
            alndb[x].to_file(aln_file)
            with open(aln_file, "r") as f:
                    lines = f.readlines()
            lines.insert(0, f"#{x}\n")
            with open(aln_file, "w") as f:
                    f.writelines(lines)
        pwd = os.getcwd()
        os.chdir(temp_dir)
        cmd = [hhdb_script,  'alns', 'hhdb_test']
        subprocess.run(cmd, check=True)
        shutil.copytree(temp_dir, os.path.join(pwd, hhm_db_name))
        os.chdir(pwd)
    return ('hmmdb produce')

def hmm_models_vs_sequence_db(models_dic_db, seqdb):
    import tempfile
    import os
    import subprocess
    import rotifer
    hmmer2table_from_rotifer  = f"{rotifer.config['base']}/bin/hmmer2table"
    with tempfile.TemporaryDirectory() as temp_dir:
        sequence_db = os.path.join(temp_dir, 'sequenceDB.fa')
        output_raw = os.path.join(temp_dir, 'output.txt')
        tbl_output = os.path.join(temp_dir, 'tbl.tsv')

        seqdb.to_file(sequence_db)
        for x in models_dic_db.keys():
            hmm_file = os.path.join(temp_dir, f'{x}.hmm')
            models_dic_db[x].to_file(hmm_file)
            result = subprocess.run(['hmmsearch', hmm_file,  sequence_db], capture_output=True, text=True)
            with open(output_raw, 'a') as file:
                file.write(result.stdout)
        tbl_result = subprocess.run([hmmer2table_from_rotifer, '-y', output_raw], capture_output=True, text=True)        
        with open(tbl_output, 'a') as file:
            file.write(tbl_result.stdout)
            
        with open(output_raw, 'r') as file:
            output = file.read()
        tbl_o = pd.read_csv(tbl_output, sep="\t")    
            
    return [output, tbl_o]            


def pid_ncbi_to_uniprotkb(ncbi_ids, max_retries=3, delay=20):
    from unipressed import IdMappingClient
    import pandas as pd
    import time
    from rotifer.devel.beta.sequence import sequence
    if isinstance(ncbi_ids, list):
        checker = 'list'
        pass
    elif isinstance(ncbi_ids,sequence):
        seqobj = ncbi_ids.copy()
        checker = 'sequenceobj'
        ncbi_ids = ncbi_ids.df.id.dropna().unique().tolist()
    # Create an ID mapping client
    request = IdMappingClient.submit(
        source="RefSeq_Protein",
        dest="UniProtKB",
        ids=set(ncbi_ids)
    )
    time.sleep(delay)
    for attempt in range(1, max_retries + 1):
        try:
            results = list(request.each_result())  # Attempt fetching results
            results =  pd.DataFrame(results).rename({"from":"id", "to": "AFDB_id"}, axis=1)
            if checker == "sequenceobj":
                seqobj.df = seqobj.df.merge(results, how='left')
                return seqobj
            else:
                return results
  # Convert to DataFrame if successful
        except Exception  as e:
            if attempt < max_retries:
                print(f"Attempt {attempt} failed: {e}. Retrying in {delay} seconds...")
                time.sleep(delay)  # Wait before retrying
            else:
                print(f"Final attempt {attempt} failed. Raising error.")
                raise  # Raise error after max retries

def uniprot2fasta(uniprot_ids):
    import requests
    from rotifer.devel.beta import sequence as rdbs
    """
    Retrieve FASTA sequences from UniProtKB for given UniProt IDs.
    
    Args:
        uniprot_ids (str/list): Single UniProt ID or list of IDs
        
    Returns:
        str: FASTA formatted sequences if successful, None otherwise
    """
    # Handle single ID input
    if isinstance(uniprot_ids, str):
        uniprot_ids = [uniprot_ids]
        
    # Build API query
    url = "https://rest.uniprot.org/uniprotkb/stream"
    params = {
        "format": "fasta",
        "query": " OR ".join(f'accession:{id}' for id in uniprot_ids)
    }
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()  # Raise exception for HTTP errors
        return rdbs.sequence(response.text)
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving sequences: {e}")
        return None
def domain_seg(seqobj,
               ted_path= "/netmnt/vast01/cbb01/proteinworld/People/gian/projects/TED/TED_database.db",
               duplicates=False):
     '''
     The function is implementend only for high confidence domain prediction, don't think it worth use lowr than that
     '''
     from rotifer.devel.alpha import gian_func as gf
     import numpy as np
     import sqlite3

     ###inner function
     def parse_seg_per_conf(segdf, conf_level):
         segdf = segdf.dropna(subset=[f'{conf_level}_c'])
         segdf[f'{conf_level}_c'] = segdf[f'{conf_level}_c'].str.split(',')
         segdf = segdf.explode(f'{conf_level}_c')
         segdf = segdf[['pid','uniprot',f'{conf_level}_c']].rename({f'{conf_level}_c':'coordinates'},axis=1)
         segdf['domain_name'] = segdf.apply(lambda x : f'{conf_level}/{x.uniprot}/{x.coordinates}', axis=1)
         segdf.coordinates = segdf.coordinates.str.split('_')
         segdf = segdf.explode(['coordinates'])
         segdf[['start', 'end']] = segdf.coordinates.str.split('-', expand=True)
         segdf['confidence'] = conf_level

         return segdf

     uniref_id = gf.pid_ncbi_to_uniprotkb(seqobj.df.id.unique().tolist())
     conn = sqlite3.connect(ted_path)
     cursor = conn.cursor()
     value_to_search = [f"AF-{id}-F1-model_v4" for id in uniref_id.to.tolist()]
     placeholders = ','.join(['?' for _ in range(len(value_to_search))])
     query = f"SELECT * FROM ted WHERE AFDB_model_ID IN ({placeholders})"

     # Define the value to search for


     # Execute the query
     cursor.execute(query, value_to_search)

     # Fetch all matching rows
     results = cursor.fetchall()
     output = pd.DataFrame(results, columns=['to', 'length', 'high' ,'medium' ,'low' ,'high_c', 'medium_c', 'low_c'])
     output.to = output.to.str.extract(r'AF-([A-Z0-9]+)-F1-model_v4')
     output = uniref_id.merge(output, how='left', on='to').rename({'from':'pid', 'to':'uniprot'}, axis=1)
     output = output.replace('na', np.nan)
     output = output.dropna(subset = ['length'])
     if not duplicates:
         output = output.drop_duplicates(['pid'])
     out2 = pd.DataFrame()
     for con in ['high', 'medium', 'low']:
         out2 = pd.concat([out2, parse_seg_per_conf(output, con)])
     seg = out2.domain_name.value_counts().where(lambda x: x>1).dropna().index.tolist()
     out2['segmented'] = np.where(out2.domain_name.isin(seg), True, False)
     to_merge = seqobj.df[['id','sequence']].rename({'id':'pid'}, axis=1).copy()
     out2 = out2.merge(to_merge,on='pid', how='left') 
     out2['aln_start'] = out2.apply(lambda x: pd.Series(list(x.sequence)).rename('seq').to_frame().query('seq !="-"').reset_index().loc[int(x.start)]['index'], axis=1)
     out2['aln_end'] = out2.apply(lambda x: pd.Series(list(x.sequence)).rename('seq').to_frame().query('seq !="-"').reset_index().loc[int(x.end) -1]['index'], axis=1)
     out2 = out2.drop(['sequence'], axis=1)

     return out2

def segment_seqobj(seqobj, confidence='medium'):
    from rotifer.devel.alpha import gian_func as gf
    import rotifer.interval.utils as riu
    segtable = gf.domain_seg(seqobj)
    segtable['fake_name'] = 'x'
    if confidence =='low':
        pass
    elif confidence =="medium":
        segtable = segtable.query('confidence !="low"')
    else:
        segtable = segtable.query('confidence =="high"')

    dt = riu.filter_nonoverlapping_regions(segtable, reference='fake_name', start='aln_start', end='aln_end', criteria={"region_length":False})
    dom_dict={}
    for x,y in dt.iterrows():
        dom_dict[y.domain_name] = seqobj.slice((y.aln_start, y.aln_end))
    return [dom_dict, segtable]

def colors_for_html_function(foreground):
    from rotifer.devel.beta.sequence import config
    import yaml
    """ Function should be improved to receive any yaml file with color mapping, at the moment it is harded coded
    it should be ease to just change the cd loading
    """
    cd = {**yaml.load(open(config['html_colors']), Loader=yaml.Loader), '-':foreground, 'x':foreground, 'X':foreground}
    #Converting the color dictionary to 2 different dictinary to apply CSS style to improve page size and loading
    color_to_key = {}
    key_to_color = {}
    counter = 1

    for color in set(cd.values()):
        new_key = f'c{counter}'
        color_to_key[color] = new_key
        key_to_color[new_key] = color
        counter += 1

    # Step 2: Generate the first dictionary
    first_dict = {k: color_to_key[v] for k, v in cd.items()}

    # Step 3: Generate the second dictionary
    second_dict = key_to_color
    return (first_dict, second_dict)
def single_tax2ndf(ndf):
    r = ndf.copy()
    from ete4 import NCBITaxa
    import rotifer
    ncbi = NCBITaxa()
    tax = r.taxid.unique().tolist()
    lineages = [ncbi.get_lineage(taxid) for taxid in tax]
    names = [ncbi.get_taxid_translator(lineage) for lineage in lineages]
    z = pd.read_csv(f"{rotifer.config['base']}/share/rotifer/taxonomy/tass2.txt", sep="\t", comment='#', names=['tax']).reset_index()
    z.rename({'index':'tax_level'}, axis=1, inplace=True)
    tt = pd.DataFrame()
    for taxid, lineage, name in zip(tax, lineages, names):
         z2 = pd.Series(name).reset_index().rename({'index':'taxid', 0:'taxname'}, axis=1)
         z2['real_taxid'] = taxid
         z3 = z.merge(z2, left_on=z['tax'].str.casefold(), right_on=z2['taxname'].str.casefold())
         tt = pd.concat([tt, z3])

    tt = tt.sort_values('tax_level').drop_duplicates('real_taxid')
    r['single_tax_name'] = r.taxid.map(tt.set_index('real_taxid').taxname.to_dict())
    r['single_tax_taxid'] = r.taxid.map(tt.set_index('real_taxid').taxid.to_dict())
    return r

def mmseqs_easy_search(acc,
             cpu=96,
             aln=True,
             max_out = 10000,
             sensitivity=3,
             iterations=1,          
             db = '/netmnt/vast01/cbb01/proteinworld/data/fadb/tmp/nr.50.mmseqs.db',
             columns = "query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln,theader"):
    '''
    MMseqs easy-search that accepts sequence object. 
    '''

    import tempfile
    import subprocess
    from subprocess import Popen, PIPE, STDOUT
    from rotifer.devel.beta.sequence import sequence as sequence
    import os
    import pandas as pd
    import time

    max_out = str(max_out)
    sensitivity = str(sensitivity)
    iterations = str(iterations)
    cwd = os.getcwd()

    with tempfile.TemporaryDirectory() as tmpdirname:
        # temporary save fasta sequence file
        if isinstance (acc, sequence):
            command = ["mmseqs",
                       "easy-search",
                       f'{tmpdirname}/seqfile',
                       db,
                       f'{tmpdirname}/out.tsv',
                       f'{tmpdirname}/tmp',
                       '--db-load-mode',
                       '2',
                       '--max-seqs',
                       max_out,
                       '-s',
                       sensitivity,
                       '--num-iterations',
                       iterations,
                       '--format-output',
                       columns]
            if aln:
                acc.to_file(f'{tmpdirname}/seqfile',  output_format='stockholm') 
                Popen(f'mmseqs convertmsa {tmpdirname}/seqfile {tmpdirname}/msaconverted', 
                  stdout=PIPE,
                  shell=True
                  ).communicate()
                Popen(f'mmseqs msa2profile {tmpdirname}/msaconverted {tmpdirname}/converted_fasta', 
                  stdout=PIPE,
                  shell=True
                  ).communicate()
                Popen(f'mmseqs search {tmpdirname}/converted_fasta {db} {tmpdirname}/resultDB {tmpdirname}/tmp --db-load-mode 2 -a -s {sensitivity} --max-seqs {max_out} --num-iterations {iterations}', 
                  stdout=PIPE,
                  shell=True
                  ).communicate()
                Popen(f'mmseqs convertalis {tmpdirname}/converted_fasta {db} {tmpdirname}/resultDB {tmpdirname}/out.tsv --db-load-mode 2 --format-output {columns}', 
                  stdout=PIPE,
                  shell=True
                  ).communicate()
            else:
                acc.to_file(f'{tmpdirname}/seqfile') 
                Popen(" ".join(command),
                      stdout=PIPE,
                      shell=True
                      ).communicate()

            t = pd.read_csv(f'{tmpdirname}/out.tsv', sep="\t", names=columns.split(','))
            t2 = sequence(t.rename({'target':'id', 'taln':'sequence'}, axis=1))
            with open(f'{tmpdirname}/out.tsv') as f:
                    blast_r = f.read()
        
    return (t, blast_r, t2) 

import os
import requests
import tempfile
from Bio.PDB import PDBParser, DSSP

def get_alphafold_dssp(afid,  sequence_object=False, id_to_guide=False, chain_id="A"):
    """
    Download AlphaFold PDB, run DSSP, return sequence and secondary structure string.

    Parameters:
    - afid: str, UniProt ID (e.g. "P0DTC2")
    - chain_id: str, chain to extract (AlphaFold models always use 'A')

    Returns:
    - sequence: str, one-letter amino acids
    - ss_string: str, secondary structure (H, E, C)
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{afid}-F1-model_v4.pdb"
    filename = f"AF-{afid}-F1-model_v4.pdb"

    with tempfile.TemporaryDirectory() as tmpdir:
        pdb_path = os.path.join(tmpdir, filename)

        # Download PDB
        r = requests.get(url)
        if r.status_code != 200:
            raise ValueError(f"Failed to download AlphaFold structure for {afid}: HTTP {r.status_code}")
        with open(pdb_path, "wb") as f:
            f.write(r.content)

        # Parse structure
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(afid, pdb_path)
        model = structure[0]

        # Run DSSP
        dssp = DSSP(model, pdb_path)

        # Build sequence and secondary structure string
        sequence = ""
        ss_string = ""
        for key in dssp.keys():
            chain, res_id = key
            if chain != chain_id:
                continue
            aa = dssp[key][1]
            raw_ss = dssp[key][2]
            ss = (
                "H" if raw_ss in ("H", "G", "I") else
                "E" if raw_ss in ("E", "B") else
                "C"
            )
            sequence += aa
            ss_string += ss
        if sequence_object:
            sequence_object = sequence_object.add_pdb(id_to_guide, pdb_file=pdb_path)
            return sequence_object

    return sequence, ss_string

def polish_2_residues_df(polish_file):
    import re
    from collections import Counter

    def drop_leading_empty_columns(df):
        while True:
            df.columns = range(df.shape[1])  # Reindex columns as 0, 1, 2, ...
            if df.empty or df.shape[1] == 0:
                break
            first_col = df.columns[0]
            # Check if entire first column is empty (NaN or blank string)
            if df[first_col].apply(lambda x: pd.isna(x) or (isinstance(x, str) and x.strip() == '')).all():
                df = df.drop(columns=first_col)
            else:
                break
        return df

    def insert_tabs(line, positions):
        chars = list(line)

        # Shift mode_second to AFTER the space
        adjusted_positions = []
        for i, pos in enumerate(positions):
            if pos is None:
                continue
            # If it's the second position (index 1), insert after the space
            adjusted_pos = pos + 1 if i == 1 else pos
            if 0 <= adjusted_pos <= len(line):
                adjusted_positions.append(adjusted_pos)

        # Insert tabs in reverse order to avoid shifting issues
        for pos in sorted(set(adjusted_positions), reverse=True):
            chars.insert(pos, '\t')

        return ''.join(chars)

    # Step 1: Read lines
    with open(polish_file, "r") as f:
        lines = [line.rstrip('\n') for line in f]

    # Step 2: Extract space-after-word positions per line
    results = [[m.start(2) for m in re.finditer(r'(\w)( )', line)] for line in lines]

    # Step 3: Compute stats across all lines
    first_elements = [lst[0] for lst in results if len(lst) >= 1]
    second_elements = [lst[1] for lst in results if len(lst) >= 2]
    last_elements = [lst[-1] for lst in results if lst]

    max_first = max(first_elements) if first_elements else None
    mode_second = Counter(second_elements).most_common(1)
    mode_second = mode_second[0][0] if mode_second else None
    mode_last = Counter(last_elements).most_common(1)
    mode_last = mode_last[0][0] if mode_last else None

    # Step 4: Function to insert tabs at multiple positions

    # Step 5: Apply to each line
    new_lines = [insert_tabs(line, [max_first, mode_second, mode_last]) for line in lines]
    tsv_string = '\n'.join(new_lines)
    # Step 6: Save result in pd.DataFrame
    from io import StringIO
    df = pd.read_csv(StringIO(tsv_string), sep='\t', header=None)
    rdf = df.set_index([0])[2].str.split("", expand=True)
    rdf = drop_leading_empty_columns(rdf)
    rdf.insert(0,'start',df[1].str.strip().values)
    rdf['end'] = df[3].str.strip().fillna('').values
    rdf.columns = list((range(1,rdf.shape[1]+1)))
    rdf.index = rdf.index.str.strip()
    rdf = rdf.replace('', ' ')
    rdf.index.name =None

    return rdf


def html_highlight_aln(s):
    from rotifer.core.functions import loadConfig
    from rotifer.core  import config as CoreConfig
    import numpy as np
    cd = loadConfig(
            ':colors.html_aa_colors',
            system_path=CoreConfig['baseDataDirectory'])
    ### getting the consensus value to map the colors filling na with "  " to color white 
    d = cd[s.fillna('_').iloc[-1]]
    #d = aa_groups_colors[s.fillna('  ').iloc[-1]]
    return np.where(
        s == '  ',
        'color:"";background-color:""',
        np.where(
            s == s.iloc[-1],
            f'color:{d["fcolor"]};background-color:{d["color"]}',
            np.where(
                s.isin(d['residues']),
                f'color:{d["fcolor"]};background-color:{d["color"]}',
                f'color:black;background-color:')))

def html_highlight_consensus(s):
    from rotifer.core.functions import loadConfig
    from rotifer.core  import config as CoreConfig
    from rotifer.devel.alpha import gian_func as gf
    import numpy as np
    cd = loadConfig(
            ':colors.html_aa_colors',
            system_path=CoreConfig['baseDataDirectory'])
    d = cd[s.fillna('_').iloc[-1]]
    """TODO: Docstring for highlight_consensus.

    :arg1: TODO
    :returns: TODO

    """
    return np.where(
        s.isin(cd["ALL"]["residues"]),
        f'color:{d["fcolor"]};background-color:{d["color"]}',
        f'color:{d["fcolor"]};background-color:{d["color"]}',
        )

def veremos(aln_r,
            consensus=False,
            annotations=False,
            font_size=6,
            output='figure.pdf',
            landscape=False,
            aln_length=70
             ):
    from rotifer.core.functions import chunks
    from rotifer.devel.alpha import gian_func as gf
    from weasyprint import HTML
    import sys
    if not consensus:
        consensus = aln_r[aln_r.index.str.contains("Consensus", case=False)].index[0]
    if not annotations:
        try:
            annotations = aln_r[aln_r.index.str.startswith("#")].index.str.lstrip('#').tolist()
            aln_r.index = aln_r.index.str.lstrip('#')
        except:
            annotations = False 

    if landscape:
        landscape='landscape'
    else:
        landscape =""


    headers = {
        'selector': 'th:not(.index_name)',
        'props': f'''font-size: {font_size}px;
        text-align: left;
        font-family:"Lucida Console", Monaco, monospace;
        color:black;'''
    }
    slices = chunks(aln_r.columns, aln_length)
    to_merge = []
    for x in slices:
        sliced = aln_r.loc[:,x].copy()
        slice_consensus = ([consensus], sliced.columns)
        if annotations:
            slice_annotaion = (annotations, sliced.columns)
            con_ann = annotations + [consensus]
        else:
            con_ann = [consensus]

        slice_sequences = sliced.loc[~sliced.index.isin(con_ann)].index.tolist() 
        slice_sequences = slice_sequences + [consensus]
        slice_sequences = (slice_sequences, sliced.columns)
        if sys.version_info.minor > 8:
            df_style = sliced.loc[:,x].style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(gf.html_highlight_aln, axis=0, subset=slice_sequences).hide(axis='columns').apply(
                gf.html_highlight_consensus, subset=slice_consensus
            ).set_table_styles(
                [headers]
            )
        else:
            df_style = sliced.loc[:,x].style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(gf.html_highlight_aln, axis=0, subset=slice_sequences).hide_columns().apply(
                gf.html_highlight_consensus, subset=slice_consensus
            ).set_table_styles(
                [headers]
            )
        df_style = df_style.to_html(table_attributes='cellspacing=0, cellpadding=0')
        to_merge.append(df_style)
    ht = ''.join(to_merge)        
    #if whant to send to latex, replace set_stick... to:to_latex(environment='longtable', convert_css=True)
    html_full = f"""
    <html>
    <head>
      <meta charset="utf-8">
      <style>
      @page {{
      size: A4 {landscape};
      margin: 1cm;
        }}

        table, th, td {{
          border: none !important;
          border-collapse: collapse;
        }}
      </style>
    </head>
    <body>
      {ht}
    </body>
    </html>
    """
    HTML(string=html_full).write_pdf(f"{output}")


def add_afdb2seqobj(seqobj, afid, model_sequence):
    from rotifer.devel.alpha.gian_func import get_alphafold_dssp
    seq_obj_copy = seqobj.sort(by=['list'], id_list = [model_sequence]).copy()
    afseq, afdssp = get_alphafold_dssp(afid)
    model_sequence_aa = seqobj.filter(keep=[model_sequence]).df.sequence[0]
    model_sequence_Series = pd.Series(list(model_sequence_aa))
    start_index = afseq.find(model_sequence_aa.replace('-',''))
    end_index =  len(model_sequence_aa.replace('-','')) + start_index
    aa = afseq[start_index:end_index]
    dssp = afdssp[start_index:end_index]
    aa = pd.Series(list(aa))
    aa.index = model_sequence_Series.where(lambda x : x!='-').dropna().index
    dssp = pd.Series(list(dssp))
    dssp.index = model_sequence_Series.where(lambda x : x!='-').dropna().index
    model_sequence_Series.update(dssp)
    dssp_str = ''.join(list(model_sequence_Series))
    data = {'id':'SS', 'sequence':dssp_str, 'type':'DSSP_from_AF', 'description':f'AF_model:{afid} ncbipd:{model_sequence}'}
    seq_obj_copy.df = pd.concat([pd.DataFrame([data]),seq_obj_copy.df])
    return seq_obj_copy
def neighborhood_DF_max_distance(df, max_distance):
    keep_rows = set()
    for block_id, block in df.groupby("block_id"):
        for strand in [+1, -1]:
            sub = block[block['strand'] == strand].copy()
            if sub.empty:
                continue
            # Sort for biological order:
            sub = sub.sort_values('start' if strand == 1 else 'end', ascending=(strand == 1))
            query_idxs = sub.index[sub['query'] == 1].tolist()
            for query_idx in query_idxs:
                # Always keep the query
                keep_rows.add(query_idx)
                curr_idx = sub.index.get_loc(query_idx)

                # Walk backwards (upstream for strand)
                prev_idx = curr_idx
                while prev_idx > 0:
                    i = prev_idx
                    j = i - 1
                    if strand == 1:
                        dist = sub.iloc[i]['start'] - sub.iloc[j]['end']
                    else:
                        dist = sub.iloc[j]['start'] - sub.iloc[i]['end']
                    if dist <= max_distance:
                        keep_rows.add(sub.index[j])
                        prev_idx = j
                    else:
                        break

                # Walk forwards (downstream for strand)
                next_idx = curr_idx
                while next_idx < len(sub) - 1:
                    i = next_idx
                    j = i + 1
                    if strand == 1:
                        dist = sub.iloc[j]['start'] - sub.iloc[i]['end']
                    else:
                        dist = sub.iloc[i]['start'] - sub.iloc[j]['end']
                    if dist <= max_distance:
                        keep_rows.add(sub.index[j])
                        next_idx = j
                    else:
                        break
    # Return only kept rows (in original DataFrame order)
    return df.loc[sorted(keep_rows)].copy()

# pip install "https://github.com/etetoolkit/ete/archive/ete4.zip"  # ETE4
# (or) pip install ete3


def taxids_to_rank_dataframe(
    taxids,
    ranks=("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"),
    include_taxid_cols=False,
    rank_synonyms=None,
    ncbi_instance=None,
):
    from collections import defaultdict
    import pandas as pd

    try:
        # Prefer ETE4
        from ete4.ncbi_taxonomy import NCBITaxa
    except Exception:
        # Fallback to ETE3
        from ete3 import NCBITaxa
    """
    Build a DataFrame (indexed by NCBI taxid) with one column per requested taxonomic rank.

    Parameters
    ----------
    taxids : Iterable[int]
        List/iterable of NCBI taxonomy IDs.
    ranks : Iterable[str], optional
        Which ranks to include as columns. Default:
        ('superkingdom','kingdom','phylum','class','order','family','genus','species').
    include_taxid_cols : bool, optional
        If True, besides the taxon name also includes columns '{rank}_taxid'
        with the corresponding NCBI taxid. Default: False.
    rank_synonyms : dict[str,str], optional
        Mapping from rank synonyms to canonical names
        (e.g., {'division':'phylum','domain':'superkingdom'}). If None, a default map is used.
    ncbi_instance : NCBITaxa, optional
        An existing NCBITaxa instance (to reuse for speed). If None, one is created.

    Returns
    -------
    pandas.DataFrame
        DataFrame indexed by 'taxid' with one column per requested rank.
        If include_taxid_cols=True, also contains '{rank}_taxid'.

    Notes on Taxonomic Ranks (NCBI Taxonomy)
    ----------------------------------------
    The NCBI hierarchy may include many intermediate levels. Common ones are:

    - superkingdom (≈domain; synonym 'domain') → e.g. Bacteria, Archaea, Eukaryota
    - kingdom
    - subkingdom
    - superphylum
    - phylum (sometimes called 'division' in plants/algae; synonym normalized here)
    - subphylum
    - superclass
    - class
    - subclass
    - infraclass
    - superorder
    - order
    - suborder
    - infraorder
    - parvorder
    - superfamily
    - family
    - subfamily
    - tribe
    - subtribe
    - genus
    - subgenus
    - species group
    - species subgroup
    - species
    - subspecies
    - variety
    - form

    Depending on the lineage, some or none of these may be present. The function
    normalizes 'division' → 'phylum' and 'domain' → 'superkingdom'.

    Disclaimer
    ----------
    This function (including docstring) was generated with assistance from GPT-5.

    """

    ncbi = ncbi_instance or NCBITaxa()

    if rank_synonyms is None:
        rank_synonyms = {
            "division": "phylum",
            "domain": "superkingdom",
        }

    def _canon_rank(rank: str):
        if not rank:
            return None
        r = rank.strip().lower()
        return rank_synonyms.get(r, r)

    ranks = tuple(ranks)
    taxids = [int(t) for t in taxids]

    # Collect lineages
    lineage_map = {}
    all_lineage_taxids = set()
    for t in taxids:
        try:
            lin = ncbi.get_lineage(t) or []
        except Exception:
            lin = []
        lineage_map[t] = lin
        all_lineage_taxids.update(lin)

    # Bulk fetch ranks and names
    ranks_all = ncbi.get_rank(list(all_lineage_taxids)) if all_lineage_taxids else {}
    names_all = ncbi.get_taxid_translator(list(all_lineage_taxids)) if all_lineage_taxids else {}

    rows = []
    for t in taxids:
        lin = lineage_map.get(t, [])
        per_rank_tid = {}
        per_rank_name = {}

        for tid in lin:
            rk = _canon_rank(ranks_all.get(tid))
            if rk in ranks and rk not in per_rank_tid:
                per_rank_tid[rk] = tid
                per_rank_name[rk] = names_all.get(tid)

        # If the queried node itself is at a requested rank
        self_rank = _canon_rank(ranks_all.get(t))
        if self_rank in ranks and self_rank not in per_rank_tid:
            per_rank_tid[self_rank] = t
            per_rank_name[self_rank] = names_all.get(t) or ncbi.get_taxid_translator([t]).get(t)

        row = {rk: per_rank_name.get(rk) for rk in ranks}
        if include_taxid_cols:
            for rk in ranks:
                row[f"{rk}_taxid"] = per_rank_tid.get(rk)
        row["taxid"] = t
        rows.append(row)

    return pd.DataFrame(rows).set_index("taxid")


def dict_to_yaml(data: dict, filename: str):
    import yaml
    
    """
    Save a Python dictionary to a YAML file.

    Args:
        data (dict): The dictionary to be saved.
        filename (str): Path to the output YAML file.
    """
    with open(filename, "w", encoding="utf-8") as f:
        yaml.dump(
            data,
            f,
            default_flow_style=False,
            sort_keys=False,
            allow_unicode=True
        )


from typing import Any, Union
import yaml
from pathlib import Path
def yaml2dict(path: Union[str, Path]) -> Any:
    """
    Lê um arquivo YAML e retorna seu conteúdo como dict (ou list).

    Parameters
    ----------
    path : str | Path
        Caminho para o arquivo YAML.

    Returns
    -------
    dict | list
        Conteúdo do YAML carregado.
    """
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Arquivo não encontrado: {path}")

    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    return data
