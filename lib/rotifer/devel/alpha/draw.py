def arch_2_svg(df, length=False, font_size=8):
    import yaml
    import numpy as np
    import seaborn as sns
    import pandas as pd
    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import rotifer.interval.utils as riu
    from rotifer.devel.alpha.sequence import complete

    #loading the domain dict to select the domains that shoul be kept and rename it (if nescssary)
    #with open('../network/laks_domains.yaml', 'r') as f:
    #    domain_dict = yaml.load(f, Loader=yaml.SafeLoader)

    # As reference Nature's standard figure sizes are 89 mm (3.50 inches) wide (single column) and 183 mm (7.20472 inches) wide (double column) 
    # If a protein with 5000 aminoacids would fit a whole line of doublw colum figure, each aa would then use 0.001440944 inches.
    # Therefor we could use this number as conversion factor to draw the domains shapes.
    # Example : a domain of 300 aa would have a shape of 300 *0.001440944 = 0.4322832 inches
    #Therefore to use calculate the number of characters that fit a shape of given size we have to use the followin forumule:
    # [size of shape]/(font_point * [character point numeber])
    # Using the full_space_character one can easyly check 

    # Size in inches of each point in font size, to calculate how much fit n the shapes:
    font_point = float(0.013837)
    full_space_character = '━'
    aa_scale = float(0.001440944) # This scale would be enough for a protein of 5000 aa on a double column figure (seems better to scalonate *4)

    colors = sns.color_palette("pastel").as_hex() +sns.color_palette("deep").as_hex()  + sns.color_palette('muted').as_hex() + sns.color_palette('colorblind')
    shapes = ['box',
              'octagon',
              'circle',
              'hexagon',
              'ellipse',
              'diamond',
              'doubleoctagon',
              'tripleoctagon']
    color_code = pd.DataFrame(shapes, columns=['shapes']).join(pd.DataFrame(colors, columns=['colors']), how='cross')

    a = df.copy() 
    if add_length and 'length' not in df.columns:
        import rotifer.devel.beta.sequence as rdbs
        a['length'] = a.ID.map(rdbs.sequence(a.ID.drop_duplicates().tolist()).df.set_index('id').length.to_dict())


    a = complete.complete(a, values={'domain':'unk'})
    d = complete.complete(a, values={'domain':'unk'})

    a = a.reset_index()
    #a.rename({'index':'pid_order'}, axis=1, inplace=True)
    a['pid_order'] = a.ID.map({ k:v for (k,v) in zip(a.ID,range(0,len(a)))})
    a = a.reset_index(drop=True).reset_index()

    #Selecting the colors
    b = a.domain.value_counts().reset_index().join(color_code, how='left')
    cc = b.set_index('index').colors.to_dict()
    sc = b.set_index('index').shapes.to_dict()

    #Preparing to create different lines

    #m = a.groupby('pid')['index'].apply(lambda x: min(x))
    #li = list(a.groupby('pid')['index'].apply(lambda x: list(x)))
    #for x in range(len(li)):
    #    li[x].insert(0,a.iloc[li[x]].pid.unique()[0])

    #Creating the graph.
    A = pgv.AGraph()
    # adding nodes  some edges

    li =[]
    gb = a.groupby('pid_order')
    for x  in gb:
        l=[]
        for y,z in x[1].iterrows():
            A.add_node(z['index'],  label=z.domain, shape=sc[z['domain']], style='round,filled',  color='lightgray', fillcolor=cc[z['domain']])
            l.append(z['index'])    
        li.append(l)
        li.append(x[1].ID.drop_duplicates().tolist())
        A.add_node(x[1].ID.unique()[0],  label=x[1].ID.unique()[0], shape='none')
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")

    # v = list(m)
    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.05)
    A.draw("subgraph7.svg", prog="dot")

    d = d.query('domain !="unk"')
    d['size']  = d.end - d.start
    d['esc'] = d['size'] * (aa_scale*4)
    #d.domain = d.domain.replace(domain_dict['Domain'])
    d['dl'] = d.domain.str.len()
    d['f_space']= round(d.esc/(font_point*font_size))
    d = d[d.f_space >= 1]
    d.f_space = d.f_space.astype(int)
    d.loc[d.domain =='unk', 'domain'] = d.loc[d.domain =='unk', 'f_space'].apply(lambda x: x * full_space_character)
    d['pid_order'] = d.ID.map(d.ID.drop_duplicates().reset_index(drop=True).reset_index().set_index('ID')['index'].to_dict())
    d = d.reset_index(drop=True).reset_index()
    #if size_median:
    #     d['esc'] = d.groupby('domain').esc.transform('median')

    gb = d.groupby('pid_order')
    li =[]
    A = pgv.AGraph()
    for x in gb:
        l =[]
        for y,z in x[1].iterrows():
            if z.domain.__contains__(full_space_character):
                A.add_node(z['index'],  label=z.domain, width=z['esc'],fixedsize='true', shape='none',height='0.4', fontsize=font_size)
                l.append(z['index'])
            else:
                A.add_node(z['index'],  label=z.domain,style='rounded,filled', width=z['esc'],fixedsize='true',height='0.4', color='lightgray', fontsize=font_size)
                l.append(z['index'])
        li.append(l)
        li.append(x[1].ID.drop_duplicates().tolist())
        A.add_node(x[1].ID.unique()[0],  label=x[1].ID.unique()[0], shape='none')
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")
        
    v = [li[x][0] for x in range(len(li))]
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.graph_attr.update(ranksep= 0.05)
    A.draw("subgraph10.svg", prog="dot")

def operon_fig(df, domain_dict=False):
    ''' saves a svg figure of neighbourhood dataframe
    '''

    import pygraphviz as pgv
    from pygraphviz.agraph import re
    import yaml
    import numpy as np
    import pandas as pd
    fontsize = 8    
    height=0.4

    def split_domain(df, column='pfam',fill ='unk', domain_dict=domain_dict, strand=False, after=10, before=10, remove_tm=False):
        '''
        It will use the query as anchor to delect the given number of after and before genes to
        create a dataframe of domains.
        One can use a dcitionary of domain to select the ones that would be kept and rename it
        '''
        def reverse_domains(pid):
            def fix_direction(df):
                df = df.copy()
                q = df.query('query ==1').query('strand ==-1').block_id
                df.strand = np.where(df.block_id.isin(q), df.strand * -1, df.strand)
                return df
            pid = fix_direction(pid)
            strand = pid.strand.unique()[0]
            if strand == -1:
                pid.domp = pid.domp.iloc[::-1].to_list()
            return pid.reset_index(drop=True)

        domdf = df.reset_index(drop=True).select_neighbors(strand=strand, after=after, before=before).copy().reset_index(drop=True).query('type == "CDS"')
        domdf['dom'] = domdf[column].str.split('+')
        domdf = domdf.explode('dom')
        domdf['domp'] = pd.Series(np.where(domdf.pid == domdf.pid.shift(1), 1,0))
        domdf.domp = domdf.groupby(['pid','block_id'])['domp'].cumsum()

        if fill:
            domdf.dom = domdf.dom.fillna(fill)

        # removing tm sig
        if remove_tm:
            domdf['tmc'] = np.where(domdf.dom.isin(["TM","SIG"]),0,1)
            only_tm_sig = list(domdf.groupby('pid').tmc.sum().where(lambda x : x==0).dropna().index)
            domdf.dom = np.where(domdf.pid.isin(only_tm_sig), fill, domdf.dom)
            domdf = domdf.query('dom no in ["TM", "SIG"]')

        # Iterate instead of groupby().apply(): pandas 3 drops the 'pid'
        # grouping column from the frames apply() returns.
        domdf = pd.concat(
                [reverse_domains(g) for _, g in domdf.groupby('pid')]
                ).reset_index(
                        drop=True
                        ).sort_values([
            'nucleotide',
            'block_id',
            'start',
            'end',
            'domp'])

        if not domain_dict:
            return domdf
            
        #domdf = domdf[domdf.dom.isin(list(domain_dict.keys()))]
        domdf.dom = domdf.dom.replace(domain_dict)


        return domdf

    te = df.query('type =="CDS"').copy()
    te = split_domain(te)
    te['shape'] = 'rectangle'
    # Removing duplicates sequencial domain in same protein
    te = te.loc[~((te.pid.shift() == te.pid) & (te.dom.shift() == te.dom))].reset_index(drop=True)

    l = te.drop_duplicates(['block_id','pid'], keep='first').query('strand ==-1').index
    r = te.drop_duplicates(['block_id','pid'], keep='last').query('strand ==1').index
    te.loc[r,'shape'] = 'rarrow'
    te.loc[l,'shape'] = 'larrow'
    te['height'] = np.where(te['shape'].isin(['larrow', 'rarrow']),  height, None)
    te = te.reset_index(drop=True).reset_index()

    m = te.groupby('block_id')['index'].apply(lambda x: min(x))
    li = list(te.groupby('block_id')['index'].apply(lambda x: list(x)))

    A = pgv.AGraph()
    # add some edges
    for x,y in te.iterrows():
        A.add_node(y['index'],  label=y.dom, shape=y['shape'],style='filled',height=y['height'],  color='lightgray', fillcolor='lightgray', fontsize=fontsize)
        
    for x in range(len(li)):
        A.add_subgraph(li[x], rank="same")
        
    v = list(m)
    [A.add_edge(v.pop(0), v[0], penwidth = 0) for x in range(len(v)-1)]
    A.graph_attr.update(nodesep= 0.02)
    A.draw("subgraph3.svg", prog="dot")
    return print(A.string())

def _get_df_to_img(self, consensus, annotations=False, remove_gaps=False, adjust_coordinates=False):
    """
    Generate a DataFrame suitable for matplotlib image generation from a sequence alignment object.

    :param consensus: The consensus threshold used to color the alignment.
    :param annotations: List or str of annotation row IDs to include in the output.
    :param remove_gaps: ID of the sequence to use as a model for gap removal.
    :param adjust_coordinates: Whether to adjust coordinates when removing gaps.
    :return: A pandas DataFrame with sequences + annotations + consensus, gaps handled if requested.
    """

    import sys
    import pandas as pd
    from rotifer.devel.beta.sequence import sequence
    import numpy as np

    aln = self.copy()

    if remove_gaps:
        gaps, gdf = aln.gaps_to_numbers(remove_gaps, adjust_coordinates=adjust_coordinates)
        gdf.index = aln.df.query('type == "sequence"').id

    aromatic = ['F', 'Y', 'W', 'H']
    alifatic = ['I', 'V', 'L']
    hydrophobic = alifatic + ['A', 'C', 'F', 'M', 'W', 'Y']
    positive = ['H', 'K', 'R']
    negative = ['D', 'E']
    charged = positive + negative
    polar = charged + ['p', 'Q', 'N', 'S', 'T', 'C']
    alcohol = ['S', 'T']
    tiny = ['G', 'A', 'S']
    small = tiny + ['V', 'T', 'D', 'N', 'P', 'C']
    big = ['K', 'F', 'I', 'L', 'M', 'Q', 'R', 'W', 'Y', 'E']
    all_aa = ['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'H', 'C', 'P', 'K', 'R', 'D', 'E', 'Q', 'N', 'S', 'T']

    aa_groups_colors = {
        'a': [aromatic, '#2C68F3'],
        'l': [alifatic, '#2CF3EA'],
        'h': [hydrophobic, '#F3E42C90'],
        '+': [positive, '#2C68F3'],
        '-': [negative, '#F50EF195'],
        'c': [charged, '#38F50E'],
        'p': [polar, '#0EF5A150'],
        'o': [alcohol, '#AE5BF8'],
        'u': [tiny, '#EE9C0C'],
        's': [small, '#DA147750'],
        'b': [big, '#A28694'],
        '.': [all_aa, 'white'],
        'G': [all_aa, 'white'], 'A': [all_aa, 'white'], 'V': [all_aa, 'white'], 'I': [all_aa, 'white'],
        'L': [all_aa, 'white'], 'M': [all_aa, 'white'], 'F': [all_aa, 'white'], 'Y': [all_aa, 'white'],
        'W': [all_aa, 'white'], 'H': [all_aa, 'white'], 'C': [all_aa, 'white'], 'P': [all_aa, 'white'],
        'K': [all_aa, 'white'], 'R': [all_aa, 'white'], 'D': [all_aa, 'white'], 'E': [all_aa, 'white'],
        'Q': [all_aa, 'white'], 'N': [all_aa, 'white'], 'S': [all_aa, 'white'], 'T': [all_aa, 'white'],
        ' ': [all_aa, 'white'], '  ': [all_aa, 'white'], '_': [all_aa, 'white']
    }

    # Getting the residues table
    aln_r = aln.residues
    con = pd.Series(list(aln.consensus(consensus)))
    aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
    con.index += 1
    aln_r = pd.concat([aln_r, con.rename('consensus').to_frame().T], axis=0)

    # Adding annotations if present
    if annotations:
        if isinstance(annotations, str):
            ann = pd.Series(list(aln.df.query('id == @annotations').sequence.iloc[0]))
            ann.index += 1
            aln_r = pd.concat([ann.rename(annotations).to_frame().T, aln_r])
        else:
            for x in annotations:
                ann = pd.Series(list(aln.df.query('id == @x').sequence.iloc[0]))
                ann.index += 1
                aln_r = pd.concat([ann.rename(x).to_frame().T, aln_r])

    # Removing gaps if specified
    if remove_gaps:
        aln_r = aln_r.drop(gaps, axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0', '  ')

    return aln_r


def aln2img(seqobj, outfile, consensus=60, remove_gaps=False, annotations=False):
    import matplotlib.pyplot as plt
    import numpy as np

    aromatic = ['F', 'Y', 'W', 'H']
    alifatic = ['I', 'V', 'L']
    hydrophobic = alifatic + ['A', 'C', 'F', 'M', 'W', 'Y']
    positive = ['H', 'K', 'R']
    negative = ['D', 'E']
    charged = positive + negative
    polar = charged + ['Q', 'N', 'S', 'T', 'C']
    alcohol = ['S', 'T']
    tiny = ['G', 'A', 'S']
    small = tiny + ['V', 'T', 'D', 'N', 'P', 'C']
    big = ['K', 'F', 'I', 'L', 'M', 'Q', 'R', 'W', 'Y', 'E']
    all_aa = ['G', 'A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W', 'H', 'C', 'P', 'K', 'R', 'D', 'E', 'Q', 'N', 'S', 'T']

    aa_groups_colors = {
        'a': [aromatic, '#2C68F3'],
        'l': [alifatic, '#2CF3EA'],
        'h': [hydrophobic, '#F3E42C'],
        '+': [positive, '#2C68F3'],
        '-': [negative, '#F50EF1'],
        'c': [charged, '#38F50E'],
        'p': [polar, '#0EF5A1'],
        'o': [alcohol, '#AE5BF8'],
        'u': [tiny, '#EE9C0C'],
        's': [small, '#DA1477'],
        'b': [big, '#A28694'],
        '.': [all_aa, '#FFFFFF']
    }

    aln_r = _get_df_to_img(seqobj, consensus=consensus, remove_gaps=remove_gaps, annotations=annotations)

    aln_r = aln_r.fillna('    ').astype(str)
    aln_r = aln_r.T.reset_index(drop=True).T
    aln_r = aln_r.reset_index().drop_duplicates('index').set_index('index')
    aln_r.index.name = None

    xx = aln_r.loc['consensus'].to_frame()
    no_color_row = np.where(aln_r.index.str.contains('ss_from', regex=True))[0].tolist()
    conservedAA = xx.iloc[:, 0].str.isupper().where(lambda x: x == True).dropna().index.to_list()
    cell_text = aln_r.values.tolist()

    fig, ax = plt.subplots()
    ax.axis('auto')
    ax.axis('off')

    the_table = ax.table(
        cellText=cell_text,
        cellLoc='center',
        rowLabels=aln_r.index.to_list(),
        colLabels=None,
        loc='center'
    )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(10)
    the_table.auto_set_column_width(col=aln_r.columns.to_list())

    for x in the_table.get_celld():
        the_table.get_celld()[x].set_facecolor('white')
        the_table.get_celld()[x].set_edgecolor('white')

        if x[1] in xx[(xx['consensus'] != '.')].index:
            if str(xx.iloc[x[1], 0]).isupper():
                if the_table.get_celld()[x].get_text().get_text() == xx.iloc[x[1], 0]:
                    the_table.get_celld()[x].set_facecolor('black')
                    the_table.get_celld()[x].set_edgecolor('black')
                    the_table.get_celld()[x].get_text().set_color('white')
            elif xx.iloc[x[1], 0] in aa_groups_colors:
                if the_table.get_celld()[x].get_text().get_text() in aa_groups_colors[xx.iloc[x[1], 0]][0]:
                    color = aa_groups_colors[xx.iloc[x[1], 0]][1]
                    the_table.get_celld()[x].set_facecolor(color)
                    the_table.get_celld()[x].set_edgecolor(color)

        if x[0] in no_color_row:
            val = the_table.get_celld()[x].get_text().get_text()
            if val in ['H', 'h']:
                the_table.get_celld()[x].set_edgecolor('red')
                the_table.get_celld()[x].get_text().set_color('white')
                the_table.get_celld()[x].set_facecolor('red')
            elif val in ['e', 'E']:
                the_table.get_celld()[x].set_edgecolor('blue')
                the_table.get_celld()[x].get_text().set_color('white')
                the_table.get_celld()[x].set_facecolor('blue')
            else:
                the_table.get_celld()[x].set_edgecolor('white')
                the_table.get_celld()[x].get_text().set_color('black')
                the_table.get_celld()[x].set_facecolor('white')

    plt.savefig(outfile, bbox_inches='tight', dpi=300)
    return f'{outfile} saved on the working path'


