import os
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
import math
from multiprocessing import Pool, pool

logger = rotifer.logging.getLogger(__name__)

# Characters Easel reads as gaps ('-', '.', '_') or missing data ('~'). Left in a
# sequence to be searched, HMMER scores them as residues: hits are lost and the
# coordinates count alignment columns instead of residues (issue #137).
GAP_CHARACTERS = '-._~'
_UNGAP = str.maketrans('', '', GAP_CHARACTERS)

def _is_pressed(hmm_file):
    """
    Tell whether an HMM database was pressed with hmmpress.

    `HMMFile.is_pressed` is a method in recent pyhmmer releases (0.12.3, at
    least) and a property in older ones. Testing the method without calling it
    is always true, which sends unpressed files to `optimized_profiles()` and
    makes it fail, so support both forms.
    """
    pressed = hmm_file.is_pressed
    if callable(pressed):
        pressed = pressed()
    return bool(pressed)

def _sequence_rows(seqobj):
    """
    Rows of a sequence object that hold sequences, leaving out rows such as
    consensus lines added by `add_consensus`.
    """
    df = seqobj.df
    if 'type' in df.columns:
        df = df[df.type == 'sequence']
    return df

def _digitize(records, ungap):
    """
    Digitize (name, sequence) pairs with the amino acid alphabet.

    Returns
    -------
    tuple
        The list of digital sequences and the number of sequences that had
        gaps removed (always 0 when `ungap` is false).
    """
    abc = ph.easel.Alphabet.amino()
    digital_sequences = []
    ungapped = 0
    for name, seq in records:
        seq = str(seq)
        if ungap:
            residues = seq.translate(_UNGAP)
            if len(residues) != len(seq):
                ungapped += 1
            seq = residues
        text_seq = ph.easel.TextSequence(name=str(name).encode(), sequence=seq)
        digital_sequences.append(text_seq.digitize(abc))
    return digital_sequences, ungapped

def digitalize_seqobj(seqobj, alignment=None, msa_name='alignment'):
    """
    Convert sequences from a sequence object into digital format for use with pyhmmer.

    This function takes a sequence object (expected to have a `.df` attribute with
    columns "id" and "sequence"), digitizes the sequences using a aminoacid alphabet,
    and optionally returns a pyhmmer DigitalMSA object.

    Only rows of type "sequence" are digitized: rows such as the consensus
    lines added by `add_consensus` are left out.

    Gaps ('-', '.', '_', '~') are removed unless `alignment` is true, so that
    an aligned sequence object can be searched like unaligned sequences. With
    `alignment=True` the gaps are kept, as building an HMM needs them.

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
            - Otherwise: returns a list of `DigitalSequence` objects, one per
              sequence, without gaps.

    Example:
        # Digitize sequences without creating a DigitalMSA
        digital_seqs = digitalize_seqobj(seqobj)

        # Create a DigitalMSA for HMM building
        msa = digitalize_seqobj(seqobj, alignment=True, msa_name='my_alignment')
    """
    rows = _sequence_rows(seqobj)
    digital_sequences, _ = _digitize(zip(rows.id, rows.sequence), ungap=not alignment)

    if alignment:
        abc = ph.easel.Alphabet.amino()
        return ph.easel.DigitalMSA(abc, sequences=digital_sequences, name=msa_name.encode())
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

def pyhmmer_to_df(output, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to','aln_target_length','aln_hmm_length','env_from','env_to']):
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

    return df

def hmmscan_linear(sequences, file=None, models_path=['/databases/pfam/Pfam-A.hmm'], cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):
    
    '''
    Perform an hmmscan of protein sequences against a Pfam HMM database.
    This function accepts an sequence object or a FASTA file path, and 
    returns a pandas DataFrame summarizing domain hits and alignment attributes.

    Gaps are removed from the input before the search, so an alignment (a
    sequence object or an aligned FASTA file) gives the same results as the
    unaligned sequences, and coordinates refer to the ungapped sequences.
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

    #Sequences load, once for all models
    seqs = _load_sequences(sequences, file)

    results = []
    #HMM load
    for model in models_path:
        with ph.plan7.HMMFile(model) as hmm_file:
           if _is_pressed(hmm_file):
               hmms = list(hmm_file.optimized_profiles())
           else:
               hmms = list(hmm_file)

        #Hmmscan run and file processment
        pbar = tqdm(total=len(seqs), desc='hmmscan')
        h = list(ph.hmmer.hmmscan(seqs, hmms, cpus=cpus, callback=callback))
        df = pyhmmer_to_df(h, columns=columns)
        df['source'] = model 
        results.append(df)

    df = pd.concat(results)
    if rename:
        df.rename({'aln_target_name': 'sequence', 'aln_hmm_name': 'model', 'i_evalue': 'evalue', 'env_from': 'estart', 'env_to': 'eend'}, axis=1, inplace=True)
    return df

def add_arch_to_df(df, column='pid', file=None, name='archtecture', evalue_filter=0.1, score_filter=0, overlap_filter = 0.1, 
                   models_path=['/databases/pfam/Pfam-A.hmm'], inplace=False, workers=4, cpus_per_worker=8):
    '''
    Add a column pfam with the domain architecture for the input accessions.
    '''
    from copy import deepcopy
    if inplace == False:
        df = df.copy()

    h = hmmscan(df[column].dropna().tolist(), workers=workers, cpus_per_worker=cpus_per_worker, file=file, models_path=models_path)
    h = h[h['evalue'] <= evalue_filter]
    h = h[h['score'] >= score_filter]

    arch = riu.architecture(h, **{ **riu.config['hmmer'], 'maximum_overlap':overlap_filter })
    arch = arch.set_index('sequence').architecture.to_dict()

    df[name] = df[column].map(arch)
    return None if inplace else df

def hmmsearch(models_path, query_db, cpus=0, columns=['aln_target_name', 'aln_hmm_name','i_evalue','c_evalue','score','env_score','aln_target_from','aln_target_to', 'aln_target_length', 'aln_hmm_length', 'env_from', 'env_to'], rename=True):
    if isinstance(models_path, str):
        models_path = [models_path]

    results = []
    for model in models_path:
        with ph.plan7.HMMFile(model) as hmm_file:
           if _is_pressed(hmm_file):
               hmms = list(hmm_file.optimized_profiles())
           else:
               hmms = list(hmm_file)
        db = ph.easel.SequenceFile(query_db, digital=True, alphabet=ph.easel.Alphabet.amino())
        out = list(ph.hmmer.hmmsearch(hmms, db, cpus=cpus))
        df = pyhmmer_to_df(out, columns=columns)
        df['source'] = model 
        results.append(df)

    df = pd.concat(results)
    if rename:
        df.rename({'aln_target_name': 'sequence', 'aln_hmm_name': 'model', 'i_evalue': 'evalue', 'env_from': 'estart', 'env_to': 'eend'}, axis=1, inplace=True)
    return df

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
            if _is_pressed(hmm_file):
                hmms = list(hmm_file.optimized_profiles())
            else:
                hmms = list(hmm_file)
        HMM_MODELS[model] = hmms

def _load_sequences(sequences, file):
    """
    Load and digitalize protein sequences, without gaps.

    Gaps are removed from every sequence, so an alignment is searched like
    the unaligned sequences. A warning is logged when any were found.

    Parameters
    ----------
    sequences : sequence object or list
        Sequence container used in the rotifer ecosystem.
    file : str or None
        Path to a FASTA file, aligned or not.

    Returns
    -------
    list
        List of digital pyhmmer sequences.
    """
    if file:
        # Easel's FASTA parser rejects gap characters, digital or not, so
        # read the file with Biopython and remove the gaps here.
        from Bio import SeqIO
        records = ((r.id, str(r.seq)) for r in SeqIO.parse(file, 'fasta'))
    else:
        if isinstance(sequences, list):
            seqobj = rdbs.sequence(sequences)
        elif isinstance(sequences, rdbs.sequence):
            seqobj = sequences
        else:
            raise ValueError("Unsupported sequence input type")
        rows = _sequence_rows(seqobj)
        records = zip(rows.id, rows.sequence)

    seqs, ungapped = _digitize(records, ungap=True)
    if ungapped:
        logger.warning(
            f'{ungapped} of {len(seqs)} input sequences contained gaps: gaps were '
            'removed before scanning and coordinates refer to the ungapped sequences.'
        )
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
    df = pyhmmer_to_df(hits, columns=columns)
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
    • Aligned input is searched without its gaps

    Gaps are removed from the input before the search, so an alignment (a
    sequence object or an aligned FASTA file) gives the same results as the
    unaligned sequences, and coordinates refer to the ungapped sequences. A
    warning is logged when gaps were removed.

    Parameters
    ----------
    sequences : sequence object or list, optional
        Sequence container to scan. Ignored if `file` is provided.
        Only rows of type "sequence" are scanned (consensus rows are not).

    file : str, optional
        Path to FASTA file containing protein sequences, aligned or not.

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
    pbar = tqdm(total=total_sequences, desc="hmmscan", unit="seq")
    for df, processed in pool.imap_unordered(_hmmscan_worker, tasks):
        results.append(df)
        pbar.update(processed) # update by number of sequences processed in that task
    pbar.close()
    pool.close()
    pool.join()

    # --------------------------------------------------------------
    # MERGE RESULTS
    # --------------------------------------------------------------
    df = pd.concat(results, ignore_index=True)
    if rename:
        df.rename({'aln_target_name': 'sequence', 'aln_hmm_name': 'model', 'i_evalue': 'evalue', 'env_from': 'estart', 'env_to': 'eend'}, axis=1, inplace=True)
    return df
