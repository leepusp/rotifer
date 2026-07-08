from asyncio import tasks
import os
from io import StringIO
import sys
import numpy as np
import pandas as pd
import pyhmmer as ph
import rotifer
import rotifer.devel.beta.sequence as rdbs
from rotifer.devel.alpha import trsantos as rdat
from rotifer.db import ncbi
from rotifer.taxonomy import utils as rtu
from rotifer.interval import utils as riu
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ete3
from tqdm import tqdm
import subprocess
from pathlib import Path
from collections.abc import Iterable
import tempfile
import math
from joblib import Parallel, delayed
from multiprocessing import Pool, pool
from rotifer.genome import utils as rgu
from rotifer.genome import io as rgio
from rotifer.devel.alpha import igem as rdai

def get_matrix(df, filter_list, rows, columns, n=10, filter_by='pid'):
        filtered_df = df[df[filter_by].isin(filter_list)]
        frequent_rows = filtered_df[rows].value_counts(rows).nlargest(n).index.tolist()
        filtered_df = filtered_df[filtered_df[rows].isin(frequent_rows)]
        matrix = filtered_df.pivot_table(index=rows, columns=columns, aggfunc='size', fill_value=0)
        return matrix

def update_lineage(ndf, preferred="/home/leep/epsoares/projects/databases/data/preferred_taxa.txt"):
    '''
    A function to update the lineage of Gene Neighborhood cursor;
    '''
    tc = ncbi.TaxonomyCursor()
    tax = tc.fetchall(ndf.taxid.dropna().drop_duplicates().tolist())
    tax['taxid'] = tax.taxid.astype(int)
    ndf.classification = ndf.taxid.map(tax.set_index('taxid').classification.to_dict())
    ndf.lineage = rtu.lineage(ndf.classification, preferred_taxa=[ x.strip() for x in open(preferred,"rt")]).tolist()
    return ndf

def extract_envelope(df, seqobj=None, start='estart', end='eend', expand=10, local_database_path='/databases/fadb/nr/nr'):
    '''
    Function to extract the envelope of a model match.
    '''
    if len(df) == 0:
        print(f'Empty input')
        seqobj = rdbs.sequence()

    if seqobj == None:
        seqobj = rdbs.sequence(df.sequence.drop_duplicates().to_list(), local_database_path=local_database_path)

    seqobj.df = seqobj.df.merge(df.rename({'sequence':'id'},axis=1), how='left')
    seqobj.df['end'] = seqobj.df[end].fillna(seqobj.df.length).astype(int)
    seqobj.df['start'] = seqobj.df[start].fillna(1).astype(int)

    if expand:
        seqobj.df['start'] = np.where(seqobj.df.start<=expand,1,seqobj.df.start-expand+1)
        seqobj.df['end'] = np.where(seqobj.df.end+expand > seqobj.df.length, seqobj.df.length, seqobj.df.end+expand)

    seqobj.df.sequence = seqobj.df.apply(lambda row: row['sequence'][row['start'] - 1:row['end']], axis=1)
    seqobj.df['pid'] = seqobj.df['id']
    seqobj.df['id'] = seqobj.df['id'] + '/' + seqobj.df[start].astype(str) + '-' + seqobj.df[end].astype(str)
    seqobj._reset()

    return seqobj

def to_network(df, target=['pfam'], ftype=['CDS'], interaction=True, ignore = [], strand = True):
    if isinstance(ftype, str):
        ftype = [ftype]
    if isinstance(target, str):
        itarget = [target]

    w = df.copy()
    if strand:
        w = df.neighbors(strand='same')
        w['rid'] = list(range(1,len(w)+1))
        w.rid = w.rid * w.strand
        w.sort_values(['rid'], inplace=True)

    # Building the source column
    w = w.query(f'type == "{ftype}"').block_id.reset_index().drop('index', axis=1)
    w['source'] = df[target[0]]
    for col in target[1:]:
        w['source'] = np.where(w['source'].isna(), df[df.type == ftype][col], w['source'])
    w['source'] = w['source'].fillna("?").str.split('+')
    w = w.explode(column='source')
    if ignore:
        w = w[~w.source.isin(ignore)].copy()

    # Building target data
    w['tblock_id'] = w['block_id'].shift(-1)
    w['target'] = w['source'].shift(-1)

    # Selecting same block rows
    w = w[w.block_id == w.tblock_id.shift(1)].copy()

    # Fix source target order when not restricted to the same strand
    sameprotein = (w.index.to_series() == w.index.to_series().shift(1))
    if not strand:
        reverse = (~sameprotein) & (w.source > w.target)
        w.loc[reverse,['source','target']] = w.loc[reverse, ['target', 'source']].values

    if interaction:
        w['interaction'] = np.where(sameprotein, 'fusion', 'neighbor')
        w = w.groupby(['source', 'target', 'interaction'])
    else:
        w = w.groupby(['source', 'target'])

    w = w.agg(weight=('block_id', 'count'), blocks=('block_id', 'nunique')).reset_index()
    return w

def make_palette(categories):
    """
    Generate a color-blind safe palette dictionary 
    for the given list of category names (max 10 entries).
    """
    safe_colors = [
        "#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#999999",
        "#117733", "#882255"
    ]

    # Keep only as many as the safe palette allows
    categories = categories[:len(safe_colors)]

    palette = {cat: safe_colors[i] for i, cat in enumerate(categories)}

    return palette

def attribute_table(
    ndf,
    columns_list=['pid', 'pfam', 'compact_total', 'compact_same_strand', 'kingdom', 'phylum', 'class', 'classification'],
    tax_parser= True,
    color_by_tax='phylum',
    number_of_taxa=5,
    save=None):

    '''
    Make a attribute table using the Neighborhood dataframe with a architecure column, 
    select and reorder columns to match those required by TreeViewer (leaf identifier 
    as first column). Also permits automatic coloring by taxonomy.
    '''

    ndfc = ndf.compact()
    ndfcs = ndf.select_neighbors(strand = True).compact()
    att = ndf[ndf.pid.isin(ndf.query('query == 1').pid)].drop_duplicates('pid')
    att['compact_total'] = att.block_id.map(ndfc['compact'].to_dict())
    att['compact_same_strand'] = att.block_id.map(ndfcs['compact'].to_dict())

    if tax_parser:
        tax = rdat.taxon_summary(att).set_index('taxid')
        for col in tax.columns:
            att[col] = att['taxid'].map(tax[col].to_dict())

    if columns_list:
        att = att[columns_list]

    if color_by_tax:
        palette = make_palette(att[color_by_tax].value_counts().head(number_of_taxa).index.tolist())
        att['Color'] = att[color_by_tax].map(palette).fillna('#000000')

    if save:
        att.to_csv(save, sep = '\t', index = False)
        print(f'Table saved to {save}')

    else:
        return att

def make_heatmap(
        ndf,
        name = None,
        domain_list = [],
        format_table = True,
        color_list = ['white','blue','green','red'],
        cbar = False,
        fmt = 'd',
        annot = False,
        linewidths = 1,
        tree_file = None,
        tree = False,
        figsize=None
    ):
    '''
    Creates a heatmap from a dataframe output of full annotate
    to represent main occurrence. The colors are correspondent 
    to the values in the dataframe, by exemple, four colors 
    corresponds to the values: 0, 1, 2 and 3, respectively. 
    The columns in the dataframe will be the columns in the 
    heatmap, and each line represents a occurrence from the 
    respective domain and the colors the type of correspondence, 
    that is determined by the numbers utilized. 0 equals to non 
    correlation, 1 correlation by neighborhood, 2 correlation by 
    fusions and 3 by both. If a tree file is delivered by the user
    the heatmap will sort the occurrences by the order of the tree.
    '''
 
    if format_table == True:
        hm = ndf.set_index('pid')[['block_id']]
        hm['neighbors'] = hm.block_id.map(ndf.query('query == 0').groupby('block_id').agg(neighbors = ('pfam','sum')).neighbors.to_dict())
        hm['pfam'] = hm.block_id.map(ndf.query('query == 1').set_index('block_id').pfam.to_dict())
        for x in domain_list:
            inpfam = hm['pfam'].str.contains(x, na = False)
            inneighbors = hm['neighbors'].str.contains(x, na = False)
            hm[f'{x}'] = ((inpfam & inneighbors).astype(int)*3 + (inpfam & ~inneighbors).astype(int)*2 + (~inpfam & inneighbors).astype(int)*1 + (~inpfam & ~inneighbors).astype(int)*0)
        df = hm[domain_list]
                                                                   
    if tree == True:
        t = ete3.Tree(f'{tree_file}')
        df = df.reindex(df.reindex(t.get_leaf_names()).index.str.replace('\'','')).fillna(0).astype(int)
    
    if figsize:
        figsize=figsize
    else:
        figsize=(len(df.columns),len(df)/5)
  
    cmap = colors.ListedColormap(color_list)
    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(data = df, cmap = cmap, cbar = cbar, fmt = fmt, annot = annot, linewidths = linewidths, ax=ax)
    plt.savefig(f'{name}', bbox_inches='tight')
    plt.close(fig)

def digitalize_seqobj(seqobj, alignment=None, msa_name='alignment'):

    """
    Convert sequences from a sequence object into digital format for use with pyhmmer.

    This function takes a sequence object (expected to have a `.df` attribute with
    columns "id" and "sequence"), digitizes the sequences using a aminoacid alphabet,
    and optionally returns a pyhmmer DigitalMSA object.

    Args:
        seqobj: Object
        alignment: bool, optional (default: None)
            If True, returns a DigitalMSA object containing all digitized sequences.
            If False or None, returns a list of pyhmmer DigitalSequence objects.
        msa_name: str, optional (default: 'alignment')
            Name to assign to the DigitalMSA object if `alignment=True`.

    Returns:
        pyhmmer.easel.DigitalMSA or list
            - If `alignment=True`: returns a `DigitalMSA` object containing all
              digitized sequences, suitable for building an HMM.
            - Otherwise: returns a list of `DigitalSequence` objects, one per sequence.

    Example:
        # Digitize sequences without creating a DigitalMSA
        digital_seqs = digitalize_seqobj(seqobj)

        # Create a DigitalMSA for HMM building
        msa = digitalize_seqobj(seqobj, alignment=True, msa_name='my_alignment')
    """

    abc = ph.easel.Alphabet.amino()
    digital_sequences = []
    
    for _, row in seqobj.df.iterrows():
        name = str(row["id"]).encode()
        seq = str(row["sequence"])
        text_seq = ph.easel.TextSequence(name=name, sequence=seq)
        digital_sequences.append(text_seq.digitize(abc))
           
    if alignment:
        for dseq in digital_sequences:
            msa = ph.easel.DigitalMSA(abc, sequences=digital_sequences, name = (msa_name.encode()))    
        return msa
    
    else:
        return digital_sequences
        
def hmmbuild(seqobj, hmm_name='alignment', save='alignment.hmm'):
      
    '''
    Build a Hidden Markov Model (HMM) from a sequence object using pyhmmer.
    This function digitizes sequences from a aligned sequence object builds an 
    HMM using pyhmmer's Plan7 Builder, and optionally saves the HMM to disk.

    Args:
        seqobj: Object
        hmm_name: str, optional (default: 'alignment')
            Name assigned to the MSA/HMM. Used internally in the digitalization step.
        save: str or None, optional (default: 'alignment.hmm')
            File path to save the resulting HMM. If None or False, the HMM is not saved.

    Returns:
        tuple of pyhmmer.plan7.HMM
            The HMM object(s) generated from the MSA.

    Example:
        # Build and save an HMM from a sequence object
        hmm = make_hmm(seqobj, hmm_name='my_model', save='my_model.hmm')

        # Build an HMM without saving to disk
        hmm = make_hmm(seqobj, save=None)
    '''
    
    abc = ph.easel.Alphabet.amino()
    msa = digitalize_seqobj(seqobj, alignment=True, msa_name=hmm_name)
    builder = ph.plan7.Builder(abc)
    background = ph.plan7.Background(abc)
    hmm = builder.build_msa(msa, background)
     
    if save:
        with open (save, 'wb') as x:
            hmm[0].write(x)
        print(f'HMM saved in {save}')
        
    return hmm
    
def pyhmmer_to_df(output, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):

    # Creation of the list to store the results
    r = []

    # Loop to process each hit
    for tophits in output:
        target_name = tophits.query.name if tophits.query.name else None
        for hits in tophits:
            for domain in hits.domains:
                r.append({
                    # pyhmmer.plan7.Domain attributes
                    "hit":                   domain.hit,
                    "bias":                  domain.bias,
                    "c_evalue":              domain.c_evalue,
                    "correction":            domain.correction,
                    "env_from":              domain.env_from,
                    "env_to":                domain.env_to,
                    "env_score":             domain.envelope_score,
                    "i_evalue":              domain.i_evalue,
                    "pvalue":                domain.pvalue,
                    "score":                 domain.score,

                    # pyhmmer.plan7.Alignment attributes
                    "aln_domain":            domain.alignment.domain,
                    "aln_hmm_accession":     domain.alignment.hmm_accession,
                    "aln_hmm_from":          domain.alignment.hmm_from,
                    "aln_hmm_name":          domain.alignment.hmm_name,
                    "aln_hmm_sequence":      domain.alignment.hmm_sequence,
                    "aln_hmm_to":            domain.alignment.hmm_to,
                    "aln_hmm_length":        domain.alignment.hmm_length,
                    "aln_identity_sequence": domain.alignment.identity_sequence,
                    "aln_target_from":       domain.alignment.target_from,
                    "aln_target_name":       domain.alignment.target_name,
                    "aln_target_sequence":   domain.alignment.target_sequence,
                    "aln_target_to":         domain.alignment.target_to,
                    'aln_target_length':     domain.alignment.target_length
                    })

    df = pd.DataFrame(r)

    if df.empty:
        print("No results found in HMMER output.")
        return pd.DataFrame(columns=columns)

    if columns:
        df = df[columns]

    if rename:
        df.rename({'aln_target_name': 'sequence', 'aln_hmm_name': 'model', 'i_evalue': 'evalue', 'env_from': 'estart', 'env_to': 'eend'}, axis=1, inplace=True)

    return df

def hmmscan_linear(sequences, file=None, models_path=['/databases/pfam/Pfam-A.hmm'], cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):
    
    '''
    Perform an hmmscan of protein sequences against a Pfam HMM database.
    This function accepts an sequence object or a FASTA file path, and 
    returns a pandas DataFrame summarizing domain hits and alignment attributes.
    ----------
    Parameters
    ----------
    sequences: sequence object
    file: str, optional
        Path to a FASTA file containing the query sequences. If provided, `sequences`
        is ignored and sequences are read from this file directly.
    pfam_database_path : str, optional
        Filesystem path to the Pfam-A HMM database file.
    cpus: int, optional
        Number of CPU threads to allocate for the hmmscan search. A value of 0
        lets HMMER autodetect available cores.
    columns: list of str, optional
        Subset of result column names to include in the output DataFrame.
        Default columns include basic domain and alignment metrics.
    '''
    #Progress bar callback
    def callback(hmm, hits):
        pbar.update(1)
    
    if isinstance(models_path, str):
        models_path = [models_path]

    results = []
    #HMM load
    for model in models_path:
        with ph.plan7.HMMFile(model) as hmm_file:
           if hmm_file.is_pressed:
               hmms = list(hmm_file.optimized_profiles())
           else:
               hmms = list(hmm_file)
	    
        #Sequences load
        abc = ph.easel.Alphabet.amino()
        if file:
           seqs = ph.easel.SequenceFile(file, digital = True, alphabet=abc)
        
        else:
           if type(sequences) == list:
               seqobj = rdbs.sequence(sequences)
               seqs = digitalize_seqobj(seqobj)

           elif type(sequences) == rotifer.devel.beta.sequence.sequence:
               seqs = digitalize_seqobj(sequences)
           
        pbar = tqdm(total=len(seqs), desc='hmmscan')
        
        #Hmmscan run and file processment
        h = list(ph.hmmer.hmmscan(seqs, hmms, cpus=cpus, callback=callback))
        df = pyhmmer_to_df(h, columns=columns, rename=rename)
        df['source'] = model 
        results.append(df)
        
    dfs = pd.concat(results)
    
    return dfs

def filter_models_overlaps(df, overlap_filter=0.1):
    """
    Resolve overlapping domains across models using a vectorized greedy approach.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain: ['estart', 'eend', 'score']
    overlap_threshold : float
        Max allowed overlap fraction (default=0.1 → 10%)

    Returns
    -------
    pd.DataFrame
        Filtered non-overlapping domains
    """

    if df.empty:
        return df

    # Sort by score (descending)
    df = df.sort_values('score', ascending=False).reset_index(drop=True)

    starts = df['estart'].to_numpy()
    ends   = df['eend'].to_numpy()
    lengths = ends - starts

    n = len(df)
    keep = np.ones(n, dtype=bool)

    for i in range(n):
        if not keep[i]:
            continue

        # Compute overlap with ALL remaining intervals
        overlap_start = np.maximum(starts[i], starts)
        overlap_end   = np.minimum(ends[i], ends)
        overlap_len   = np.maximum(0, overlap_end - overlap_start)

        # Normalize by the smaller interval (robust choice)
        min_len = np.minimum(lengths[i], lengths)
        frac_overlap = np.zeros_like(overlap_len, dtype=float)
        valid = min_len > 0
        frac_overlap[valid] = overlap_len[valid] / min_len[valid]

        # Suppress overlapping domains BELOW in ranking
        mask = (frac_overlap > overlap_filter)

        # Only remove those with lower priority (j > i)
        mask[:i+1] = False

        keep[mask] = False

    return df[keep]

def add_arch_to_df(df, column='pid', file=None, column_arch_name='arch', evalue_filter=0.1, score_filter=0, overlap_filter = 0.1, 
                   models_path=['/databases/pfam/Pfam-A.hmm'], inplace=False, run_hmmscan=True, workers=4, cpus_per_worker=8,
                   build_arch_by_source=False):

    '''
    Add a column pfam with the domain architecture for the input accessions.
    '''

    if inplace == False:
        df = df.copy()

    if run_hmmscan:
        h = hmmscan(df[column].dropna().tolist(), workers=workers, cpus_per_worker=cpus_per_worker, file=file, models_path=models_path)

    else:
        h = df

    h.rename({'aln_target_name':'sequence','aln_hmm_name':'model','i_evalue':'evalue','env_from':'estart', 'env_to':'eend'}, axis=1, inplace=True)
    h = h.loc[:, ~h.columns.duplicated()]
    h = h.drop_duplicates().reset_index(drop=True)
    h = h[h['evalue'] <= evalue_filter]    
    h = h[h['score'] >= score_filter]
    
    if build_arch_by_source == True:
        sources = h.source.drop_duplicates().str.split('/').str[-1].tolist()
        for x in sources:
            h_source = h[h.source.str.contains(x)]
            h_source = riu.filter_nonoverlapping_regions(h_source, **riu.config['hmmer'])
            h_source = h_source.groupby('sequence', group_keys=False).apply(filter_models_overlaps, overlap_filter=overlap_filter)
            h_source = h_source.sort_values(['sequence', 'estart'])
            arch = h_source.groupby('sequence').agg(pfam = ('model',lambda x: '+'.join(x.astype(str)))).reset_index()
            arch.rename({'sequence':column}, axis = 1, inplace = True)
            arch = arch.set_index(column).pfam.to_dict()
            df[f'{column_arch_name}_{x}'] = df[column].map(arch)
        return df

    else:            
        h = riu.filter_nonoverlapping_regions(h, **riu.config['hmmer'])
        h = h.groupby('sequence', group_keys=False).apply(filter_models_overlaps, overlap_filter=overlap_filter)
        h = h.sort_values(['sequence', 'estart'])
        arch = h.groupby('sequence').agg(pfam = ('model',lambda x: '+'.join(x.astype(str)))).reset_index()
        arch.rename({'sequence':column}, axis = 1, inplace = True)
        arch = arch.set_index(column).pfam.to_dict()
        df[column_arch_name] = df[column].map(arch)
        return df
      
def hmmsearch(models_path, query_db, cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):
    
    if isinstance(models_path, str):
        models_path = [models_path]

    results = []
    
    for model in models_path:
        with ph.plan7.HMMFile(model) as hmm_file:
           if hmm_file.is_pressed:
               hmms = list(hmm_file.optimized_profiles())
           else:
               hmms = list(hmm_file)

        db = ph.easel.SequenceFile(query_db, digital=True, alphabet=ph.easel.Alphabet.amino())
        out = list(ph.hmmer.hmmsearch(hmms, db, cpus=cpus))
        df = pyhmmer_to_df(out, columns=columns, rename=rename)
        df['source'] = model 
        results.append(df)
        
    dfs = pd.concat(results)

    return dfs

def _compute_chunk_size(total_sequences, workers, task_factor=4):
    """
    Compute an adaptive chunk size for balanced parallel workloads.

    Parameters
    ----------
    total_sequences : int
        Number of sequences to process.

    workers : int
        Number of parallel worker processes.

    task_factor : int, optional
        Multiplier controlling how many tasks are created per worker.
        Higher values improve load balancing but increase scheduling overhead.

    Returns
    -------
    int
        Recommended chunk size.
    """

    total_tasks = workers * task_factor

    chunk_size = math.ceil(total_sequences / total_tasks)

    return max(chunk_size, 1)

def _load_models(models_path):
    """
    Load HMM models from the provided database paths.

    This function loads all HMM profiles once in the parent process.
    Worker processes created afterward inherit these objects through
    copy-on-write memory sharing.

    Parameters
    ----------
    models_path : list[str]
        Paths to HMM database files.

    Returns
    -------
    None
        Models are stored in the global variable `HMM_MODELS`.
    """

    global HMM_MODELS

    HMM_MODELS = {}

    for model in models_path:

        with ph.plan7.HMMFile(model) as hmm_file:

            if hmm_file.is_pressed:
                hmms = list(hmm_file.optimized_profiles())
            else:
                hmms = list(hmm_file)

        HMM_MODELS[model] = hmms

def _load_sequences(sequences, file):
    """
    Load and digitalize protein sequences.

    Parameters
    ----------
    sequences : sequence object or list
        Sequence container used in the rotifer ecosystem.
    file : str or None
        Path to a FASTA file.

    Returns
    -------
    list
        List of digital pyhmmer sequences.
    """

    abc = ph.easel.Alphabet.amino()

    if file:
        seqs = list(
            ph.easel.SequenceFile(
                file,
                digital=True,
                alphabet=abc
            )
        )

    else:

        if type(sequences) == list:
            seqobj = rdbs.sequence(sequences)
            seqs = list(digitalize_seqobj(seqobj))

        elif type(sequences) == rotifer.devel.beta.sequence.sequence:
            seqs = list(digitalize_seqobj(sequences))

        else:
            raise ValueError("Unsupported sequence input type")

    return seqs

def _chunk_sequences(seqs, chunk_size):
    """
    Split sequence list into chunks.

    Parameters
    ----------
    seqs : list
        List of sequences.
    chunk_size : int
        Number of sequences per chunk.

    Returns
    -------
    list
        List of sequence chunks.
    """

    return [
        seqs[i:i + chunk_size]
        for i in range(0, len(seqs), chunk_size)
    ]

def _hmmscan_worker(args):

    model, seq_chunk, cpus_per_worker, columns, rename = args

    hmms = HMM_MODELS[model]

    hits = list(
        ph.hmmer.hmmscan(
            seq_chunk,
            hmms,
            cpus=cpus_per_worker
        )
    )

    df = pyhmmer_to_df(
        hits,
        columns=columns,
        rename=rename
    )

    df["source"] = model

    # return both dataframe and number of sequences processed
    return df, len(seq_chunk)

def hmmscan(
    sequences=None,
    file=None,
    models_path=['/databases/pfam/Pfam-A.hmm'],
    workers=4,
    cpus_per_worker=8,
    chunk_size=None,
    columns=[
        'aln_target_name',
        'aln_hmm_name',
        'i_evalue',
        'c_evalue',
        'score',
        'env_score',
        'aln_target_from',
        'aln_target_to',
        'aln_target_length',
        'aln_hmm_length',
        'env_from',
        'env_to'
    ],
    rename=True
):
    """
    Perform hmmscan searches against one or more HMM databases.

    This implementation supports high-performance parallel execution
    suitable for very large datasets (hundreds of thousands to millions
    of sequences).

    Key features
    ------------
    • Shared HMM memory across workers
    • Parallel sequence chunk processing
    • Configurable CPU usage per worker
    • Compatible with both sequence objects and FASTA input
    • Efficient scaling on multi-core systems

    Parameters
    ----------
    sequences : sequence object or list, optional
        Sequence container to scan. Ignored if `file` is provided.

    file : str, optional
        Path to FASTA file containing protein sequences.

    models_path : str or list[str], optional
        Path(s) to HMM database files.

    workers : int, optional
        Number of parallel worker processes.
        Default: 4.

    cpus_per_worker : int, optional
        Number of CPU threads used internally by hmmscan
        within each worker process.
        Default = 8
        Total CPU usage ≈ workers × cpus_per_worker

    chunk_size : int, optional
        Number of sequences processed per job.
        Larger chunks:
        • lower scheduling overhead
        • higher memory usage

    columns : list[str], optional
        Columns retained from hmmer output.

    rename : bool, optional
        Whether to apply column renaming in the parser.

    Returns
    -------
    pandas.DataFrame
        Combined hmmscan results for all sequences and models.
    """

    # --------------------------------------------------------------
    # MODEL PATH NORMALIZATION
    # --------------------------------------------------------------

    if isinstance(models_path, str):
        models_path = [models_path]

    # --------------------------------------------------------------
    # LOAD HMM DATABASES
    # --------------------------------------------------------------

    _load_models(models_path)

    # --------------------------------------------------------------
    # LOAD SEQUENCES
    # --------------------------------------------------------------

    seqs = _load_sequences(sequences, file)

    # --------------------------------------------------------------
    # CHUNK SEQUENCES
    # --------------------------------------------------------------
    if chunk_size is None:
        chunk_size = _compute_chunk_size(
            total_sequences=len(seqs),
            workers=workers)

    seq_chunks = _chunk_sequences(seqs, chunk_size)
    # --------------------------------------------------------------
    # BUILD TASK LIST
    # --------------------------------------------------------------
    tasks = []

    for model in models_path:
        for chunk in seq_chunks:
            tasks.append(
                (model, chunk, cpus_per_worker, columns, rename)
            )

    # --------------------------------------------------------------
    # PARALLEL EXECUTION
    # --------------------------------------------------------------

    pool = Pool(workers)

    results = []
    total_sequences = len(seqs) * len(models_path)

    pbar = tqdm(
        total=total_sequences,
        desc="hmmscan",
        unit="seq")

    for df, processed in pool.imap_unordered(_hmmscan_worker, tasks):
        results.append(df)

    # update by number of sequences processed in that task
        pbar.update(processed)
    pbar.close()
    
    pool.close()
    pool.join()

    # --------------------------------------------------------------
    # MERGE RESULTS
    # --------------------------------------------------------------

    dfs = pd.concat(results, ignore_index=True)

    return dfs

def run_fimo_single(meme_file, genome, out_dir=None, extra_args=None):
    """
    Execute FIMO for a single genome.

    Parameters
    ----------
    meme_file : str | Path
        MEME motif file.
    genome : str | Path
        FASTA file.
    out_dir : str | Path | None
        Output directory. If None, uses a temporary directory.
    extra_args : list[str] | None
        Additional CLI arguments for FIMO.

    Returns
    -------
    pd.DataFrame
        Parsed FIMO output with an additional 'genome' column.
    """
    meme_file = Path(meme_file)
    genome = Path(genome)

    if out_dir is None:
        out_dir = Path(tempfile.mkdtemp())
    else:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    cmd = ["fimo", "--oc", str(out_dir)]
    if extra_args:
        cmd.extend(extra_args)
    cmd.extend([str(meme_file), str(genome)])

    subprocess.run(cmd, check=True)

    df = pd.read_csv(out_dir / "fimo.tsv", sep="\t", comment="#")
    df["genome"] = genome.name

    return df


def run_fimo_batch(meme_file, genomes, extra_args=None, n_jobs=1):
    """
    Execute FIMO across multiple genomes.

    Parameters
    ----------
    meme_file : str | Path
    genomes : iterable[str | Path]
    extra_args : list[str] | None
    n_jobs : int
        Parallel jobs (uses joblib if >1)

    Returns
    -------
    pd.DataFrame
        Concatenated FIMO results.
    """
    if isinstance(genomes, (str, Path)):
        genomes = [genomes]
        genomes = list(genomes)
    else:
        genomes = list(genomes)
        
    if n_jobs == 1:
        dfs = [run_fimo_single(meme_file, g, extra_args=extra_args) for g in genomes]
    else:
        dfs = Parallel(n_jobs=n_jobs)(
            delayed(run_fimo_single)(meme_file, g, extra_args=extra_args)
            for g in genomes
        )

    return pd.concat(dfs, ignore_index=True)


def build_gff_index(gffs):
    """
    Build an index of CDS features from one or more GFF/GFF3 files.

    Parameters
    ----------
    gffs : str | Path | iterable[str | Path]
        Single GFF file or collection of GFF files.

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping:
            seqid -> dataframe containing CDS features sorted by coordinates.
    """

    # Accept single file or iterable of files
    if isinstance(gffs, (str, Path)):
        gffs = [gffs]
    elif not isinstance(gffs, Iterable):
        raise TypeError(
            "gffs must be a path or an iterable of paths"
        )

    gff_dict = {}

    columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    dtypes = {
        "seqid": "string",
        "source": "string",
        "type": "string",
        "score": "string",
        "strand": "string",
        "phase": "string",
        "attributes": "string",
    }

    for gff in map(Path, gffs):

        # Find FASTA section if present
        fasta_line = None

        with gff.open() as fh:
            for i, line in enumerate(fh):
                if line.startswith("##FASTA"):
                    fasta_line = i
                    break

        read_kwargs = dict(
            sep="\t",
            comment="#",
            header=None,
            names=columns,
            dtype=dtypes,
            low_memory=False,
        )

        # Read only annotation section when FASTA exists
        if fasta_line is not None:
            read_kwargs["nrows"] = fasta_line

        gffdf = pd.read_csv(gff, **read_kwargs)

        # Keep only valid CDS rows
        cds = gffdf.loc[gffdf["type"] == "CDS"].copy()

        if cds.empty:
            continue

        cds["start"] = pd.to_numeric(cds["start"], errors="coerce")
        cds["end"] = pd.to_numeric(cds["end"], errors="coerce")

        cds = cds.dropna(subset=["start", "end"])

        cds["start"] = cds["start"].astype("int64")
        cds["end"] = cds["end"].astype("int64")

        cds.sort_values(
            ["seqid", "start", "end"],
            inplace=True,
            kind="mergesort",
        )

        for seqid, subdf in cds.groupby("seqid", sort=False):

            subdf = subdf.reset_index(drop=True)

            if seqid in gff_dict:
                gff_dict[seqid] = pd.concat(
                    [gff_dict[seqid], subdf],
                    ignore_index=True,
                )
            else:
                gff_dict[seqid] = subdf

    # Final sort in case same seqid appeared in multiple files
    for seqid in gff_dict:

        gff_dict[seqid] = (
            gff_dict[seqid]
            .sort_values(
                ["start", "end"],
                kind="mergesort"
            )
            .reset_index(drop=True)
        )

    return gff_dict


def get_next_protein(df, gff_dict):
    """
    Annotate each FIMO hit with the nearest downstream/upstream CDS.

    Also extracts a normalized protein ID (pid) from GFF attributes.

    Parameters
    ----------
    df : pd.DataFrame
        FIMO output. Must contain:
        ['sequence_name', 'start', 'stop', 'strand']
    gff_dict : dict[str, pd.DataFrame]
        Output of build_gff_index().

    Returns
    -------
    pd.DataFrame
        Original dataframe with:
        - next_protein : raw GFF attributes
        - pid : extracted protein ID
    """

    def _get(row):
        seq = row["sequence_name"]
        if seq not in gff_dict:
            return None

        cds = gff_dict[seq]

        if row["strand"] == "+":
            hits = cds[cds["start"].values > row["stop"]]
            if hits.empty:
                return None
            return hits.iloc[0]["attributes"]

        else:
            hits = cds[cds["end"].values < row["start"]]
            if hits.empty:
                return None
            return hits.iloc[-1]["attributes"]

    df = df.copy()
    df["next_protein"] = df.apply(_get, axis=1)

    # parse ID (vectorized)
    df["pid"] = df["next_protein"].str.split(";", expand=True)[0].str.replace("ID=cds-", "", regex=False)
    #df.drop('next_protein', axis = 1, inplace = True)

    return df

def get_distances_repeats(df, inplace=True, filter=False, length=20):
    '''
    Get the distance from fimo hit for the next one, and can filter the occurrences
    due a specific value. Default is 20. 
    '''
    if inplace == False:
        df = df.copy()

    df.sort_values(['sequence_name','strand', 'start'], inplace=True)
    df['distance'] = df['start'] - df.groupby(['sequence_name','strand'])['stop'].shift(1)

    if filter == True:
        df = df[df['distance'] <= length]

    return df

def fimo_pipeline(meme_file, genomes, gffs, n_jobs=1, filter=True, length=20):
    """
    End-to-end execution:
    FIMO → annotate next protein → cluster hits.

    Returns
    -------
    pd.DataFrame
    """
    df = run_fimo_batch(meme_file, genomes, n_jobs=n_jobs)
    gff_dict = build_gff_index(gffs)
    df = get_next_protein(df, gff_dict)
    df = get_distances_repeats(df, filter=filter, length=length)

    return df

def igem_pipeline(genome_annotation, genome_format, genome_protein_fasta, genome_nucleotide_fasta, models_path=['/databases/pfam/Pfam-A.hmm'], return_hmmscan=False,
    after=10, before=10, run_fimo=True, meme_file='/home/leep/epsoares/projects/igem/2026/data/heptarepeats2.meme', return_fimo=False, make_figure=True,
    output_report='neighborhood_report.html', color_dict=None, domain_dict=None):
    ''' 
    Doc
    '''

    fimo = fimo_pipeline(meme_file, genome_nucleotide_fasta, genome_annotation).query('2 <= distance <= 15')
    gen = rgu.seqrecords_to_dataframe(rgio.parse(genome_annotation, informat=genome_format), exclude_type=['source', 'gene', 'region'])
    
    if genome_format == 'genbank':
        gen['pid'] = gen.locus

    hscan = hmmscan(file=genome_protein_fasta, models_path=models_path)
    add_arch_to_df(hscan, run_hmmscan=False, inplace=True, column='sequence')
    gen['pfam'] = gen.pid.map(hscan.set_index('sequence').pfam.to_dict())
    ndf = gen.neighbors(gen.pid.isin(fimo.sequence_name), after=after, before=before)
    
    if make_figure:
        rdai.build_html_report(ndf, output_file=output_report, custom_colors=color_dict, rename_map=domain_dict)
        print(f'figure saved in {output_report}')

    if return_fimo and return_hmmscan:
        return ndf, fimo, hscan
    elif return_fimo:
        return ndf, fimo
    elif return_hmmscan:
        return ndf, hscan

    return ndf
