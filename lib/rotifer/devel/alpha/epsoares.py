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
from rotifer.devel.alpha import malu as rdam

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

def _arch_with_coordinates(h, seq_col='sequence', model_col='model', start_col='estart', end_col='eend'):
    '''
    Map each sequence to its domain architecture annotated with the
    envelope coordinates of every domain, e.g. ``HTH[1-120]+TPR[240-300]``.

    Domains are ordered by start coordinate. Returns a ``{sequence: str}``
    dict, ready to feed ``df[column].map(...)``.
    '''
    h = h.sort_values([seq_col, start_col])

    def _fmt(g):
        return '+'.join(
            f'{model}[{int(start)}-{int(end)}]'
            for model, start, end in zip(g[model_col], g[start_col], g[end_col])
        )

    return h.groupby(seq_col).apply(_fmt).to_dict()

def add_arch_to_df(df, column='pid', file=None, column_arch_name='pfam', evalue_filter=0.1, score_filter=0, overlap_filter = 0.1,
                   models_path=['/databases/pfam/Pfam-A.hmm'], inplace=False, run_hmmscan=True, workers=4, cpus_per_worker=8,
                   build_arch_by_source=False, add_coordinates=True, column_coord_name=None):

    '''
    Add a column pfam with the domain architecture for the input accessions.

    Parameters
    ----------
    add_coordinates : bool, default False
        Also add a column holding the architecture annotated with each
        domain's envelope coordinates, e.g. ``HTH[1-120]+TPR[240-300]``.
    column_coord_name : str or None
        Name of that extra column. Defaults to ``f'{column_arch_name}_coord'``.
        When ``build_arch_by_source`` is True the source tag is appended,
        matching the architecture columns (``f'{column_coord_name}_{source}'``).
    '''

    if inplace == False:
        df = df.copy()

    if column_coord_name is None:
        column_coord_name = f'{column_arch_name}_coord'

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
            if add_coordinates:
                df[f'{column_coord_name}_{x}'] = df[column].map(_arch_with_coordinates(h_source))
        return None if inplace else df

    else:            
        h = riu.filter_nonoverlapping_regions(h, **riu.config['hmmer'])
        h = h.groupby('sequence', group_keys=False).apply(filter_models_overlaps, overlap_filter=overlap_filter)
        h = h.sort_values(['sequence', 'estart'])
        arch = h.groupby('sequence').agg(pfam = ('model',lambda x: '+'.join(x.astype(str)))).reset_index()
        arch.rename({'sequence':column}, axis = 1, inplace = True)
        arch = arch.set_index(column).pfam.to_dict()
        df[column_arch_name] = df[column].map(arch)
        if add_coordinates:
            df[column_coord_name] = df[column].map(_arch_with_coordinates(h))
        return None if inplace else df
      
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

# ------------------------------------------------------------------
# GENOME INPUT HANDLING
# ------------------------------------------------------------------
# The FIMO pipeline needs two things out of a genome: the nucleotide
# sequences FIMO scans and the CDS coordinates used to assign each hit to
# the gene it may regulate. A GFF3 with an embedded ##FASTA section carries
# both, and so does a GenBank flatfile (.gbff), so neither input should be
# required to arrive pre-split into a FASTA plus an annotation file.
#
# The helpers below normalize any of these -- including gzipped files --
# into the two objects the pipeline consumes. Parsing is delegated to
# rotifer.genome (rgio.parse/rgu.seqrecords_to_dataframe), which already
# supports GFF and every Bio.SeqIO format, instead of adding a second
# GenBank reader to this module.

_GENOME_FORMAT_BY_SUFFIX = {
    ".gff": "gff",
    ".gff3": "gff",
    ".gb": "genbank",
    ".gbk": "genbank",
    ".gbf": "genbank",
    ".gbff": "genbank",
    ".genbank": "genbank",
    ".embl": "embl",
    ".fa": "fasta",
    ".fas": "fasta",
    ".fna": "fasta",
    ".ffn": "fasta",
    ".faa": "fasta",
    ".fasta": "fasta",
}


def guess_genome_format(path):
    """
    Guess the format of a genome file from its extension.

    Compression extensions ('.gz', '.bz2') are ignored, so
    'GCF_000005845.2_genomic.gbff.gz' is recognized as GenBank.

    Parameters
    ----------
    path : str | Path

    Returns
    -------
    str | None
        'gff', 'genbank', 'embl', 'fasta' or None when unrecognized.
    """
    suffixes = [s.lower() for s in Path(str(path)).suffixes]

    while suffixes and suffixes[-1] in (".gz", ".bz2"):
        suffixes.pop()

    if not suffixes:
        return None

    return _GENOME_FORMAT_BY_SUFFIX.get(suffixes[-1])


def parse_genome(genome, informat=None, **kwargs):
    """
    Iterate over the Bio.SeqRecord objects of an annotated genome file.

    Thin wrapper around rotifer.genome.io.parse that guesses the format from
    the file name and transparently opens compressed files.

    Parameters
    ----------
    genome : str | Path
        GFF (with a ##FASTA section) or any Bio.SeqIO format, such as the
        GenBank flatfiles (.gbff) distributed by NCBI. May be gzipped.
    informat : str | None
        Format name. Guessed from the extension when None.

    Yields
    ------
    Bio.SeqRecord.SeqRecord
    """
    import shutil
    import rotifer.core.functions as rcf

    genome = str(genome)
    informat = informat or guess_genome_format(genome)

    if informat is None:
        raise ValueError(
            f"Cannot guess the format of {genome}: pass informat explicitly"
        )
    if informat == "fasta":
        raise ValueError(
            f"{genome} is a plain FASTA file and carries no annotation"
        )

    compressed = Path(genome).suffix.lower() in (".gz", ".bz2")

    # rgio.gff() takes a path and opens the file itself, while Bio.SeqIO
    # formats are fed an open stream so that compressed files are read
    # without a detour through the disk.
    handles = []
    try:
        if informat == "gff":
            if compressed:
                plain = tempfile.NamedTemporaryFile(
                    mode="wt", suffix=".gff", delete=False
                )
                handles.append(plain)
                with rcf.open_compressed(genome, mode="rt") as fh:
                    shutil.copyfileobj(fh, plain)
                plain.close()
                source = plain.name
            else:
                source = genome
        else:
            source = rcf.open_compressed(genome, mode="rt")
            handles.append(source)

        try:
            for seqrecord in rgio.parse(source, informat=informat, **kwargs):
                yield seqrecord
        except IndexError:
            if informat == "gff":
                raise ValueError(
                    f"{genome} has no ##FASTA section: its sequences must be "
                    "provided separately"
                )
            raise

    finally:
        for handle in handles:
            handle.close()
            if informat == "gff" and os.path.exists(handle.name):
                os.unlink(handle.name)


def load_genome_annotation(genome, informat=None, exclude_type=["source", "gene", "region"], **kwargs):
    """
    Parse an annotated genome into a rotifer genome dataframe.

    Equivalent to ``rgu.seqrecords_to_dataframe(rgio.parse(...))`` but with
    format guessing and transparent decompression (see parse_genome).

    Parameters
    ----------
    genome : str | Path | iterable[str | Path]
        One or more GFF/GenBank/EMBL files.
    informat : str | None
        Format name, guessed per file when None.
    exclude_type : list
        Feature types dropped while building the table.

    Returns
    -------
    rotifer.genome.data.NeighborhoodDF
    """
    if isinstance(genome, (str, Path)):
        genome = [genome]

    dfs = []
    for path in genome:
        dfs.append(
            rgu.seqrecords_to_dataframe(
                parse_genome(path, informat=informat),
                exclude_type=exclude_type,
                **kwargs,
            )
        )

    if len(dfs) == 1:
        return dfs[0]

    return pd.concat(dfs, ignore_index=True)


def genome_to_fasta(genome, informat=None, out=None):
    """
    Return the path of a nucleotide FASTA for a genome, extracting it when
    the input is an annotated file (GFF with ##FASTA, GenBank, EMBL...).

    FASTA inputs are returned untouched, so passing a plain .fna costs
    nothing; compressed FASTA is decompressed because FIMO cannot read it.

    Parameters
    ----------
    genome : str | Path
    informat : str | None
        Format name, guessed from the extension when None.
    out : str | Path | None
        Destination file. A temporary file is created when None.

    Returns
    -------
    str
        Path to an uncompressed FASTA file. It is the input path itself when
        no conversion was needed, which is how callers know whether the file
        is theirs to delete.
    """
    import shutil
    from Bio import SeqIO
    import rotifer.core.functions as rcf

    genome = str(genome)
    informat = informat or guess_genome_format(genome)
    compressed = Path(genome).suffix.lower() in (".gz", ".bz2")

    # Nothing to extract: hand FIMO the file it was given.
    if informat in (None, "fasta") and not compressed:
        return genome

    if out is None:
        handle = tempfile.NamedTemporaryFile(
            mode="wt", suffix=".fna", delete=False
        )
        out = handle.name
    else:
        out = str(out)
        handle = open(out, "wt")

    try:
        if informat in (None, "fasta"):
            with rcf.open_compressed(genome, mode="rt") as fh:
                shutil.copyfileobj(fh, handle)
        else:
            SeqIO.write(parse_genome(genome, informat=informat), handle, "fasta")
    finally:
        handle.close()

    return out


def genome_to_protein_fasta(genome, informat=None, out=None, codontable="Bacterial"):
    """
    Write the translation of every CDS of an annotated genome to a FASTA file.

    Protein identifiers follow the same rule used by
    rgu.seqrecords_to_dataframe (protein_id, else the GFF 'ID' attribute), so
    the headers match the 'pid' column of the genome dataframe and of the
    FIMO table produced by this pipeline.

    Translations come from the /translation qualifier when present
    (GenBank) and are computed from the genomic sequence otherwise, which
    means a GFF only works when it carries a ##FASTA section.

    Parameters
    ----------
    genome : str | Path
    informat : str | None
    out : str | Path | None
        Destination file. A temporary file is created when None.
    codontable : int | str
        Genetic code used when a CDS has no /translation qualifier.

    Returns
    -------
    str
        Path of the FASTA file written.
    """
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    if out is None:
        out = tempfile.NamedTemporaryFile(suffix=".faa", delete=False).name
    else:
        out = str(out)

    records = []
    for seqrecord in parse_genome(genome, informat=informat):
        for ft in seqrecord.features:
            if ft.type != "CDS" or "pseudo" in ft.qualifiers:
                continue

            qualifiers = ft.qualifiers
            pid = None
            for tag in ("protein_id", "ID"):
                if tag in qualifiers:
                    pid = qualifiers[tag][0]
                    break
            if pid is None:
                continue

            if "translation" in qualifiers:
                protein = qualifiers["translation"][0]
            else:
                table = codontable
                if "transl_table" in qualifiers:
                    table = int(qualifiers["transl_table"][0])
                try:
                    # Pass the Seq, not the SeqRecord: SeqFeature.translate
                    # mirrors the type it is given.
                    protein = str(
                        ft.translate(
                            seqrecord.seq, table=table, cds=False, to_stop=True
                        )
                    )
                except Exception:
                    continue

            records.append(
                SeqRecord(
                    Seq(protein),
                    id=pid,
                    description=qualifiers.get("product", [""])[0],
                )
            )

    if not records:
        raise ValueError(
            f"No protein translation could be extracted from {genome}"
        )

    SeqIO.write(records, out, "fasta")

    return out


def run_fimo_single(meme_file, genome, out_dir=None, extra_args=None, informat=None):
    """
    Execute FIMO for a single genome.

    Parameters
    ----------
    meme_file : str | Path
        MEME motif file.
    genome : str | Path
        Nucleotide FASTA, or any annotated genome carrying its sequences
        (GenBank/.gbff, GFF with a ##FASTA section), optionally compressed.
        Annotated inputs are converted to a temporary FASTA for FIMO.
    out_dir : str | Path | None
        Output directory. If None, uses a temporary directory.
    extra_args : list[str] | None
        Additional CLI arguments for FIMO.
    informat : str | None
        Format of `genome`, guessed from the extension when None.

    Returns
    -------
    pd.DataFrame
        Parsed FIMO output with an additional 'genome' column, which always
        names the file given by the caller rather than the temporary FASTA.
    """
    meme_file = Path(meme_file)
    genome = Path(genome)

    if out_dir is None:
        out_dir = Path(tempfile.mkdtemp())
    else:
        out_dir = Path(out_dir)
        out_dir.mkdir(parents=True, exist_ok=True)

    fasta = genome_to_fasta(genome, informat=informat)

    try:
        cmd = ["fimo", "--oc", str(out_dir)]
        if extra_args:
            cmd.extend(extra_args)
        cmd.extend([str(meme_file), str(fasta)])

        subprocess.run(cmd, check=True)

        df = pd.read_csv(out_dir / "fimo.tsv", sep="\t", comment="#")

    finally:
        # Only remove what this call created.
        if str(fasta) != str(genome):
            os.unlink(fasta)

    df["genome"] = genome.name

    return df


def run_fimo_batch(meme_file, genomes, extra_args=None, n_jobs=1, informat=None):
    """
    Execute FIMO across multiple genomes.

    Parameters
    ----------
    meme_file : str | Path
    genomes : iterable[str | Path]
        FASTA and/or annotated genomes (see run_fimo_single).
    extra_args : list[str] | None
    n_jobs : int
        Parallel jobs (uses joblib if >1)
    informat : str | None
        Format of the inputs, guessed per file when None.

    Returns
    -------
    pd.DataFrame
        Concatenated FIMO results.
    """
    if isinstance(genomes, (str, Path)):
        genomes = [genomes]
    else:
        genomes = list(genomes)

    if n_jobs == 1:
        dfs = [
            run_fimo_single(meme_file, g, extra_args=extra_args, informat=informat)
            for g in genomes
        ]
    else:
        dfs = Parallel(n_jobs=n_jobs)(
            delayed(run_fimo_single)(
                meme_file, g, extra_args=extra_args, informat=informat
            )
            for g in genomes
        )

    return pd.concat(dfs, ignore_index=True)


def build_gff_index(gffs):
    """
    Build an index of CDS features from one or more GFF/GFF3 files.

    This is the fast path for GFF: the annotation section is read straight
    into pandas, so plain GFF files -- with or without a ##FASTA section --
    are supported. Other formats go through build_annotation_index().

    Parameters
    ----------
    gffs : str | Path | iterable[str | Path]
        Single GFF file or collection of GFF files.

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping:
            seqid -> dataframe containing CDS features sorted by coordinates.

        The 'pid' column follows the same rule as
        rgu.seqrecords_to_dataframe (the protein_id attribute, else ID), so
        identifiers match those of a genome dataframe parsed from the same
        file.
    """

    import rotifer.core.functions as rcf

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

        # open_compressed, so gzipped GFFs are scanned like plain ones;
        # pandas infers the compression from the extension on its own.
        with rcf.open_compressed(str(gff), mode="rt") as fh:
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

        # Protein ID, using the same precedence as rgu.seqrecords_to_dataframe
        attributes = cds["attributes"].fillna("").astype(str)
        cds["pid"] = (
            attributes.str.extract(r"(?:^|;)protein_id=([^;]*)", expand=False)
            .fillna(attributes.str.extract(r"(?:^|;)ID=([^;]*)", expand=False))
        )

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


def cds_index_from_dataframe(gendf):
    """
    Build a CDS index from a rotifer genome dataframe.

    Accepts anything shaped like the output of rgu.seqrecords_to_dataframe
    (columns 'nucleotide', 'start', 'end', 'strand', 'pid', 'type'), which is
    how GenBank/EMBL annotation reaches this pipeline -- and also lets a
    caller that already parsed its genome reuse that table instead of
    reading the annotation a second time.

    Parameters
    ----------
    gendf : pd.DataFrame

    Returns
    -------
    dict[str, pd.DataFrame]
        Same structure as build_gff_index().
    """
    missing = [
        c for c in ("nucleotide", "start", "end", "strand", "type")
        if c not in gendf.columns
    ]
    if missing:
        raise ValueError(
            f"Not a genome dataframe, missing column(s): {', '.join(missing)}"
        )

    cds = gendf.loc[gendf["type"] == "CDS"].copy()

    index = {}
    if cds.empty:
        return index

    cds["start"] = pd.to_numeric(cds["start"], errors="coerce")
    cds["end"] = pd.to_numeric(cds["end"], errors="coerce")
    cds = cds.dropna(subset=["start", "end"])

    pid = cds["pid"] if "pid" in cds.columns else pd.Series(np.nan, index=cds.index)

    # Rebuild a GFF-like attribute string so 'next_protein' stays as
    # informative as it is for GFF input.
    attributes = "ID=" + pid.fillna("").astype(str)
    for column, tag in (("locus", "locus_tag"), ("gene", "gene"), ("product", "product")):
        if column in cds.columns:
            extra = cds[column].fillna("").astype(str)
            attributes = attributes.mask(
                extra.ne(""), attributes + ";" + tag + "=" + extra
            )

    cds = pd.DataFrame(
        {
            "seqid": cds["nucleotide"].astype(str),
            "start": cds["start"].astype("int64"),
            "end": cds["end"].astype("int64"),
            "strand": cds["strand"].map(normalize_strand),
            "pid": pid,
            "attributes": attributes,
        }
    )

    for seqid, subdf in cds.groupby("seqid", sort=False):
        index[seqid] = (
            subdf.sort_values(["start", "end"], kind="mergesort")
            .reset_index(drop=True)
        )

    return index


def _merge_cds_index(index, other):
    """
    Merge CDS index `other` into `index`, concatenating shared seqids.
    """
    for seqid, subdf in other.items():
        if seqid in index:
            index[seqid] = (
                pd.concat([index[seqid], subdf], ignore_index=True)
                .sort_values(["start", "end"], kind="mergesort")
                .reset_index(drop=True)
            )
        else:
            index[seqid] = subdf

    return index


def build_annotation_index(annotation, informat=None):
    """
    Build a CDS index from any supported annotation input.

    This is the format-agnostic entry point used by the FIMO pipeline. It
    accepts:

    * GFF/GFF3 files, read by build_gff_index()
    * GenBank (.gbff, .gbk), EMBL and any other Bio.SeqIO format, parsed
      through rotifer.genome and converted by cds_index_from_dataframe()
    * an already parsed genome dataframe (rgu.seqrecords_to_dataframe)
    * an index already built by one of the functions above, returned as is

    Compressed (.gz/.bz2) GenBank files are read directly.

    Parameters
    ----------
    annotation : str | Path | iterable | pd.DataFrame | dict
    informat : str | None
        Format of the file(s), guessed from each extension when None.

    Returns
    -------
    dict[str, pd.DataFrame]
        seqid -> CDS features sorted by coordinates.
    """
    if isinstance(annotation, dict):
        return annotation

    if isinstance(annotation, pd.DataFrame):
        return cds_index_from_dataframe(annotation)

    if isinstance(annotation, (str, Path)):
        paths = [annotation]
    elif isinstance(annotation, Iterable):
        paths = list(annotation)
    else:
        raise TypeError(
            "annotation must be a path, an iterable of paths, a genome "
            "dataframe or a CDS index"
        )

    index = {}
    for path in paths:
        fmt = informat or guess_genome_format(path)

        if fmt == "gff":
            part = build_gff_index(path)
        elif fmt in (None, "fasta"):
            raise ValueError(
                f"{path} carries no annotation: pass a GFF or GenBank file, "
                "or set informat explicitly"
            )
        else:
            part = cds_index_from_dataframe(
                load_genome_annotation(path, informat=fmt)
            )

        _merge_cds_index(index, part)

    return index


_STRAND_ALIASES = {"+": "+", "1": "+", 1: "+", "-": "-", "-1": "-", -1: "-"}


def normalize_strand(value):
    """
    Map the several strand conventions found in FIMO/GFF/rotifer tables
    ('+'/'-', 1/-1, '1'/'-1') to '+' or '-'. Anything else, missing values
    included, becomes None so it can be skipped instead of silently matching.
    """
    if isinstance(value, str):
        value = value.strip()
    elif value is None or value is pd.NA:
        return None
    elif pd.isna(value):
        return None

    return _STRAND_ALIASES.get(value)


def get_next_protein(df, annotation, max_distance=50, informat=None):
    """
    Annotate each FIMO hit with the nearest downstream CDS on the same
    strand as the repeat (i.e. the gene lying after the last repeat in the
    direction of transcription).

    Two constraints are enforced:

    * strand: only CDSs transcribed in the same direction as the repeat are
      eligible, and "downstream" is read in that direction (higher coordinates
      for '+', lower coordinates for '-'). A convergent or divergent gene that
      happens to sit closer is never picked.
    * distance: the CDS must start within ``max_distance`` bp of the end of the
      repeat. Without this cutoff any hit would be assigned the next CDS on the
      replicon, no matter how far away it is; the cutoff also means that in a
      tandem array only the last repeat -- the one actually abutting the gene --
      gets annotated, while the upstream copies are left empty.

    Parameters
    ----------
    df : pd.DataFrame
        FIMO output. Must contain:
        ['sequence_name', 'start', 'stop', 'strand']
    annotation : dict | str | Path | pd.DataFrame
        Anything build_annotation_index() accepts: a CDS index built by
        build_gff_index()/build_annotation_index(), a GFF or GenBank file, or
        a parsed genome dataframe.
    max_distance : int, default 50
        Maximum number of base pairs allowed between the end of the repeat and
        the start of the downstream CDS. Use None to disable the cutoff.
    informat : str | None
        Annotation format, guessed from the file name when None.

    Returns
    -------
    pd.DataFrame
        Original dataframe with:
        - next_protein : annotation attributes of the CDS
        - next_protein_distance : gap in bp between repeat and CDS
        - next_protein_strand : strand of the CDS, always equal to the repeat's
        - pid : protein ID of the CDS
    """
    gff_dict = build_annotation_index(annotation, informat=informat)

    # Index the CDSs by (seqid, strand) once, so a repeat can only ever look at
    # genes co-oriented with it. Rows whose strand is missing or unrecognized
    # are dropped here instead of leaking into the comparisons below.
    by_strand = {}
    for seqid, cds in gff_dict.items():
        strands = cds["strand"].map(normalize_strand).values
        for strand in ("+", "-"):
            sub = cds[strands == strand]
            if not sub.empty:
                by_strand[(seqid, strand)] = (
                    sub.sort_values(["start", "end"], kind="mergesort")
                    .reset_index(drop=True)
                )

    # Indexes built from a genome dataframe, and those built by
    # build_gff_index(), already carry the protein ID; only a hand-made index
    # would leave it out, and there it is parsed from the attributes below.
    has_pid = all("pid" in cds.columns for cds in gff_dict.values())

    def _get(row):
        strand = normalize_strand(row["strand"])
        if strand is None:
            return None, None, None, None

        cds = by_strand.get((row["sequence_name"], strand))
        if cds is None:
            return None, None, None, None

        if strand == "+":
            # Transcription runs left to right: take the first CDS starting
            # after the end of the repeat.
            hits = cds[cds["start"].values > row["stop"]]
            if hits.empty:
                return None, None, None, None
            hit = hits.iloc[0]
            distance = hit["start"] - row["stop"]

        else:
            # Transcription runs right to left: "downstream" means the last CDS
            # ending before the repeat, and the repeat's own start is its 3' end
            # in FIMO's plus-strand coordinates.
            hits = cds[cds["end"].values < row["start"]]
            if hits.empty:
                return None, None, None, None
            hit = hits.iloc[-1]
            distance = row["start"] - hit["end"]

        if max_distance is not None and distance > max_distance:
            return None, None, None, None

        pid = hit["pid"] if has_pid else None

        return hit["attributes"], distance, strand, pid

    df = df.copy()
    df[["next_protein", "next_protein_distance", "next_protein_strand", "pid"]] = (
        df.apply(_get, axis=1, result_type="expand")
    )

    if not has_pid:
        df["pid"] = (
            df["next_protein"].str.split(";", expand=True)[0]
            .str.replace("ID=cds-", "", regex=False)
        )

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

def filter_repeat_arrays(df, min_distance=2, max_distance=15, min_repeats=2, inplace=False):
    """
    Keep only repeats that belong to a valid tandem array.

    A heptarepeat is only meaningful as a regulatory element when it occurs as
    an array: consecutive copies must be spaced within
    ``min_distance < gap <= max_distance`` bp, and the array must contain at
    least ``min_repeats`` copies. Isolated hits, and hits whose spacing breaks
    the array, are discarded.

    Note this is a different condition from the repeat-to-CDS cutoff applied by
    get_next_protein(): that one says how far the gene may be from the last
    repeat, this one says how the repeats must be arranged among themselves.

    Parameters
    ----------
    df : pd.DataFrame
        FIMO output, with a 'distance' column as produced by
        get_distances_repeats(). It is computed on the fly when absent.
    min_distance : int, default 2
        Exclusive lower bound for the gap between consecutive repeats.
    max_distance : int, default 15
        Inclusive upper bound for the gap between consecutive repeats.
    min_repeats : int, default 2
        Minimum number of copies for an array to be valid.
    inplace : bool, default False
        Operate on the input dataframe instead of a copy.

    Returns
    -------
    pd.DataFrame
        Rows of the valid arrays only, with:
        - repeat_array : array identifier
        - repeat_count : number of copies in that array

        'distance' is back-filled for the first copy of each array with the gap
        to its neighbour, so every member of a valid array carries an in-window
        spacing instead of the NaN that opens the array. Otherwise the copy
        adjacent to the gene on the minus strand -- which is the leftmost, hence
        the one with no predecessor by coordinate -- would be silently dropped
        by downstream distance filters such as malu.filter_fimo().
    """
    if not inplace:
        df = df.copy()

    if "distance" not in df.columns:
        df = get_distances_repeats(df, inplace=True)

    df.sort_values(["sequence_name", "strand", "start"], inplace=True)

    # A repeat continues the previous array only when the gap to its predecessor
    # falls inside the allowed window. The first hit of each sequence/strand has
    # a NaN distance, so it always opens a new array.
    linked = (df["distance"] > min_distance) & (df["distance"] <= max_distance)
    df["repeat_array"] = (~linked).cumsum()
    df["repeat_count"] = df.groupby("repeat_array")["repeat_array"].transform("size")

    # The copy opening an array has no predecessor by coordinate; give it the
    # gap to the next copy so its spacing is defined in both orientations.
    df["distance"] = df["distance"].fillna(
        df.groupby("repeat_array")["distance"].shift(-1)
    )

    return df[df["repeat_count"] >= min_repeats]


def fimo_pipeline(meme_file, genomes=None, annotation=None, informat=None,
                  n_jobs=1, filter=True, length=20, max_distance=50, gffs=None):
    """
    End-to-end execution:
    FIMO → annotate next protein → cluster hits.

    Sequences and annotation may come from the same file. A GenBank flatfile
    (.gbff) or a GFF carrying a ##FASTA section is enough on its own:

        fimo_pipeline(meme, 'GCF_000005845.2_genomic.gbff.gz')

    while a plain GFF, which has no sequences, is passed alongside the
    nucleotide FASTA:

        fimo_pipeline(meme, 'genome.fna', 'genome.gff')

    Parameters
    ----------
    meme_file : str | Path
        MEME motif file.
    genomes : str | Path | iterable | None
        Sequences FIMO scans: FASTA, or annotated genomes carrying their
        sequences. Defaults to `annotation` when None.
    annotation : str | Path | iterable | pd.DataFrame | dict | None
        Annotation used to find the gene downstream of each hit. Accepts GFF,
        GenBank/EMBL, an already parsed genome dataframe or a CDS index (see
        build_annotation_index). Defaults to `genomes` when None.
    informat : str | None
        Format of the input files, guessed from each file name when None.
        It applies to `genomes` and `annotation` alike, so set it only when
        one file plays both roles or both share the same format.
    n_jobs : int
        Parallel FIMO jobs.
    filter : bool
        Drop hits farther than `length` bp from the previous hit.
    length : int
        Cutoff used when filter is True.
    max_distance : int, default 50
        Maximum distance, in bp, between the end of the repeat and the start of
        the downstream CDS (see get_next_protein). Use None to disable.
    gffs : deprecated
        Former name of `annotation`.

    Returns
    -------
    pd.DataFrame
    """
    if annotation is None:
        annotation = gffs

    if genomes is None and annotation is None:
        raise ValueError("Either genomes or annotation must be given")

    # A single annotated file plays both roles.
    if genomes is None:
        if isinstance(annotation, (pd.DataFrame, dict)):
            raise ValueError(
                "genomes is required: a parsed annotation has no sequences "
                "for FIMO to scan"
            )
        genomes = annotation
    if annotation is None:
        annotation = genomes

    df = run_fimo_batch(meme_file, genomes, n_jobs=n_jobs, informat=informat)
    df = get_next_protein(df, annotation, max_distance=max_distance, informat=informat)
    df = get_distances_repeats(df, filter=filter, length=length)

    return df

def igem_pipeline(genome_annotation, genome_format=None, genome_protein_fasta=None, genome_nucleotide_fasta=None, models_path=['/databases/pfam/Pfam-A.hmm', '/home/leep/epsoares/projects/igem/2026/data/all_models.hmm'],
    search_models='/home/leep/epsoares/projects/igem/2026/data/search_models.hmm', hmmsearch_score_filter=30, hmmsearch_evalue_filter=1e-4, return_hmmscan=False, after=10, before=10, run_fimo=True,
    meme_file='/home/leep/epsoares/projects/igem/2026/data/heptarepeats2.meme', return_fimo=False, make_figure=True, output_report='neighborhood_report.html', 
    repeat_max_distance=50, repeat_min_spacing=2, repeat_max_spacing=15, min_repeats=2, 
    color_dict=None, domain_dict=None, seed=4, patience=4, max_distance=50, max_extend=30, 
    domains_filter='/home/leep/epsoares/projects/igem/2026/data/hmm_modelnames.tsv', organism=None, add_sequences=True, normalize_orientation=False,
    filter_columns=['seq_type', 'assembly', 'gene', 'origin', 'topology', 'taxid', 'lineage', 'classification', 'feature_order', 'internal_id', 'is_fragment']):
    '''
    Run the whole iGEM analysis on one genome.

    The genome annotation may be a GFF or a GenBank flatfile (.gbff, also
    gzipped). When the annotation carries its sequences -- always for
    GenBank, and for a GFF with a ##FASTA section -- genome_nucleotide_fasta
    and genome_protein_fasta are optional and extracted from it.

    Parameters
    ----------
    genome_annotation : str | Path
        GFF or GenBank annotation of the genome.
    genome_format : str | None
        'gff', 'genbank', ... Guessed from the file name when None.
    genome_protein_fasta : str | Path | None
        Protein FASTA used by hmmscan/hmmsearch. Extracted from
        genome_annotation when None.
    genome_nucleotide_fasta : str | Path | None
        Nucleotide FASTA scanned by FIMO. Extracted from genome_annotation
        when None.
    '''
    genome_format = genome_format or guess_genome_format(genome_annotation)

    # Parse the annotation once and reuse the table: it drives both the
    # repeat-to-gene assignment inside fimo_pipeline() and every neighborhood
    # operation below, so their protein IDs cannot disagree, whatever the
    # annotation format.
    gen = load_genome_annotation(genome_annotation, informat=genome_format)

    derived = []
    if genome_nucleotide_fasta is None:
        genome_nucleotide_fasta = genome_to_fasta(genome_annotation, informat=genome_format)
        derived.append(genome_nucleotide_fasta)
    if genome_protein_fasta is None:
        genome_protein_fasta = genome_to_protein_fasta(genome_annotation, informat=genome_format)
        derived.append(genome_protein_fasta)

    try:
        # filter=False: filter_repeat_arrays() needs every hit, including the first
        # copy of each array (whose distance is NaN), to count array sizes.
        fimo = fimo_pipeline(meme_file, genome_nucleotide_fasta, gen, max_distance=repeat_max_distance, filter=False)
        fimo = filter_repeat_arrays(fimo, min_distance=repeat_min_spacing, max_distance=repeat_max_spacing, min_repeats=min_repeats)

        fimo = rdam.filter_fimo(fimo, gen).query('intragenic == False')
        hscan = hmmscan(file=genome_protein_fasta, models_path=models_path)
        hsearch = hmmsearch(search_models, genome_protein_fasta)
        hsearch_hits = riu.filter_nonoverlapping_regions(hsearch, **riu.config['hmmer']).query(f'score >= {hmmsearch_score_filter} and evalue <= {hmmsearch_evalue_filter}')
        l = hsearch_hits.sequence.tolist()
        pids_list = fimo.pid.dropna().tolist() + l
        add_arch_to_df(hscan, run_hmmscan=False, inplace=True, column='sequence')
        gen['pfam'] = gen.pid.map(hscan.set_index('sequence').pfam.to_dict())
        # ndf = gen.neighbors(gen.pid.isin(pids), after=after, before=before)
        ndf = rdam.filter_neighbors_plus(gen, pids=pids_list, mode='strict', annotate=False, after=after, before=before, max_distance=max_distance,
                                     max_extend=max_extend, seed=seed, patience=patience, reqdom=domains_filter)
        ndf['repeat_start'] = ndf.pid.map(fimo.set_index('pid').start.to_dict())
        ndf['repeat_end'] = ndf.pid.map(fimo.set_index('pid').stop.to_dict())
        ndf['repeat_strand'] = ndf.pid.map(fimo.set_index('pid').strand.to_dict())
        ndf['pfam_coord'] = ndf.pid.map(hscan.set_index('sequence').pfam_coord.to_dict())

        if add_sequences:
            seqs = rdbs.sequence(genome_protein_fasta)
            ndf['sequence'] = ndf.pid.map(seqs.df.set_index('id').sequence.to_dict())
            matched = int(ndf.sequence.notna().sum())
            total = int(ndf.pid.notna().sum())
            print(f'Sequences attached to {matched} of {total} proteins')
            if total and not matched:
                print(f'  WARNING: no pid matched a header in {genome_protein_fasta} -- '
                      'the report will have no sequences')

        # Tag every neighborhood with the search that recovered its query: the
        # heptarepeat MEME/FIMO scan ('Heptarepeat') or the HMM that matched it in
        # hmmsearch (the model's own name). A query found by both searches, or by
        # more than one model, carries every tag joined by '+'.
        hepta_pids = set(fimo.pid.dropna())
        model_by_pid = (hsearch_hits.astype({'model': str}).groupby('sequence')['model'].agg(lambda names: '+'.join(dict.fromkeys(names))).to_dict())

        def _query_source(pid):
            tags = ['Heptarepeat'] if pid in hepta_pids else []
            if pid in model_by_pid:
                tags.append(model_by_pid[pid])
            return '+'.join(tags) if tags else np.nan

        query_source = {pid: _query_source(pid) for pid in dict.fromkeys(pids_list)}
        ndf['query_source'] = ndf.pid.map(query_source)

        # neighbors inherit the tag(s) of their block's query row(s)
        block_source = (
            ndf.loc[ndf['query_source'].notna(), ['block_id', 'query_source']]
            .groupby('block_id')['query_source']
            .agg(lambda tags: '+'.join(dict.fromkeys('+'.join(tags).split('+'))))
            .to_dict()
        )
        ndf['query_source'] = ndf['block_id'].map(block_source)

        if filter_columns:
            ndf = ndf.drop(columns=[c for c in filter_columns if c in ndf.columns])

        if organism:
            ndf['organism'] = organism

        if make_figure:
            rdai.build_html_report(ndf, output_file=output_report, custom_colors=color_dict, rename_map=domain_dict,
                                   normalize_orientation=normalize_orientation)
            print(f'figure saved in {output_report}')

        if return_fimo and return_hmmscan:
            return ndf, fimo, hscan
        elif return_fimo:
            return ndf, fimo
        elif return_hmmscan:
            return ndf, hscan

        return ndf
    finally:
        # Remove only the FASTA files this call extracted from the annotation.
        for path in derived:
            if os.path.exists(path):
                os.unlink(path)
