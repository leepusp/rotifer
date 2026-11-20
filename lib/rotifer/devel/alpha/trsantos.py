import numpy as np
import pandas as pd
from collections import Counter
from rotifer.db import ncbi
from rotifer.db.ncbi import entrez
from rotifer.devel.alpha import epsoares as rdae
from rotifer.devel.beta import sequence as rdbs
from rotifer.interval import utils as riu
ic = ncbi.IPGCursor()
gnc = ncbi.GeneNeighborhoodCursor()
tc = ncbi.TaxonomyCursor()
import os
os.environ["MKL_THREADING_LAYER"] = "GNU"
import subprocess
import time
import psutil
import tempfile
import pathlib
import traceback
from glob import glob
from typing import Optional, Union, List, Dict
from io import StringIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo
import seaborn as sns
import pygraphviz as pgv
import html

# --- GPU support via pynvml ---
try:
    import pynvml
    pynvml.nvmlInit()
    gpu_available = True
except Exception:
    gpu_available = False

def get_gpu_stats():
    """
    Collect current GPU utilization and memory usage for all available GPUs.

    Returns
    -------
    list[dict]
        A list of dictionaries, one per GPU, each containing:
        - "gpu_id" (int): The GPU index.
        - "gpu_util" (float): Current GPU utilization percentage.
        - "gpu_mem" (float): GPU memory usage in megabytes (MB).

    Notes
    -----
    - If `pynvml` is not available or fails to initialize, returns an empty list.
    - This function assumes NVIDIA GPUs accessible via NVML.
    """
    stats = []
    if not gpu_available:
        return stats
    device_count = pynvml.nvmlDeviceGetCount()
    for i in range(device_count):
        handle = pynvml.nvmlDeviceGetHandleByIndex(i)
        util = pynvml.nvmlDeviceGetUtilizationRates(handle)
        mem = pynvml.nvmlDeviceGetMemoryInfo(handle)   # ✅ fixed
        stats.append({
            "gpu_id": i,
            "gpu_util": util.gpu,            # percent
            "gpu_mem": mem.used / (1024**2)  # MB
        })
    return stats


def profile_command(command, poll_interval=0.1, output_file=None, shell=False):
    """
    Run an external command while profiling its CPU, memory, and (if available) GPU usage.

    Parameters
    ----------
    command : list[str] | str
        The command to execute, either as a list of arguments (recommended) or as a shell string.
    poll_interval : float, optional
        Time interval (in seconds) between sampling events. Default is 0.1 seconds.
    output_file : str | None, optional
        File path to write command stdout/stderr. If None, a temporary log file is created.
    shell : bool, optional
        Whether to execute the command through the shell. Default is False.

    Returns
    -------
    dict
        A dictionary summarizing the execution and profiling results:
        - "runtime_sec" (float): Total runtime of the command in seconds.
        - "avg_cpu_percent" (float | None): Average CPU utilization percentage across all processes.
        - "peak_mem_mb" (float | None): Peak memory usage in megabytes.
        - "avg_gpu_util_percent" (float | None): Average GPU utilization percentage (if available).
        - "peak_gpu_mem_mb" (float | None): Peak GPU memory usage in megabytes (if available).
        - "return_code" (int | None): Process return code.
        - "error" (str | None): Error message if profiling failed.
        - "traceback" (str | None): Full traceback if an error occurred.
        - "output_file" (str): Path to the log file with stdout/stderr output.

    Notes
    -----
    - Uses `psutil` for CPU and memory sampling, and `pynvml` for GPU statistics.
    - Captures resource usage across all child processes spawned by the command.
    - The function blocks until the process finishes.
    
    Example
    ------
all_stats = []
for pdb_file in sorted(glob("WeakConsensus/*.pdb")):
    pdb_stem = pathlib.Path(pdb_file).stem
    out_path = f"eguchicnn_weakconsensus_output/{pdb_stem}.txt"
    log_out = f"eguchicnn_weakconsensus_logs/{pdb_stem}.log"
    pathlib.Path("eguchicnn_weakconsensus_output").mkdir(exist_ok=True)
    pathlib.Path("eguchicnn_weakconsensus_logs").mkdir(exist_ok=True)

    result = profile_command(
        [
            "python", "../../data/software/EguchiCNN/eguchicnn_wrapper.py",
            "-i", str(pdb_file),
            "-o", out_path
        ],
        poll_interval=0.05,
        output_file=log_out,
        shell=False
    )

    result.update({
        "pdb": pdb_stem,
        "eguchicnn_out": out_path,
        "log_out": log_out
    })
    all_stats.append(result)
    """
    # ensure we always have a file to write stdout/stderr to
    temp_file = None
    if output_file is None:
        tf = tempfile.NamedTemporaryFile(prefix="profile_cmd_", suffix=".log", delete=False)
        output_file = tf.name
        tf.close()
        temp_file = output_file

    start_time = time.time()

    try:
        f = open(output_file, "w")
    except Exception as e:
        return {
            "runtime_sec": None,
            "avg_cpu_percent": None,
            "peak_mem_mb": None,
            "avg_gpu_util_percent": None,
            "peak_gpu_mem_mb": None,
            "return_code": None,
            "error": f"Could not open output file {output_file}: {e}",
            "traceback": traceback.format_exc(),
            "output_file": output_file
        }

    try:
        process = subprocess.Popen(command, stdout=f, stderr=f, shell=shell)
    except Exception as e:
        f.close()
        return {
            "runtime_sec": None,
            "avg_cpu_percent": None,
            "peak_mem_mb": None,
            "avg_gpu_util_percent": None,
            "peak_gpu_mem_mb": None,
            "return_code": None,
            "error": f"Popen failed: {e}",
            "traceback": traceback.format_exc(),
            "output_file": output_file
        }

    # attach psutil for CPU/mem tracking
    try:
        p = psutil.Process(process.pid)
        # prime cpu_percent so that the first nonzero measurement works
        p.cpu_percent(interval=None)
    except Exception:
        p = None

    cpu_samples, mem_samples = [], []
    gpu_samples, gpu_mem_samples = [], []

    try:
        while True:
            if process.poll() is not None:
                break

            try:
                if p is not None:
                    procs = [p] + p.children(recursive=True)
                    total_cpu = sum(proc.cpu_percent(interval=None) for proc in procs if proc.is_running())
                    total_mem = sum(proc.memory_info().rss for proc in procs if proc.is_running())
                else:
                    total_cpu, total_mem = None, None

                if total_cpu is not None:
                    cpu_samples.append(total_cpu)
                if total_mem is not None:
                    mem_samples.append(total_mem / (1024**2))  # MB

                if gpu_available:
                    gstats = get_gpu_stats()
                    if gstats:
                        avg_gpu = np.mean([g["gpu_util"] for g in gstats])
                        sum_gmem = np.sum([g["gpu_mem"] for g in gstats])
                        gpu_samples.append(avg_gpu)
                        gpu_mem_samples.append(sum_gmem)

            except psutil.NoSuchProcess:
                break

            time.sleep(poll_interval)

        return_code = process.wait()

    except Exception as e:
        try:
            process.kill()
        except Exception:
            pass
        return_code = getattr(process, "returncode", None)
        f.close()
        return {
            "runtime_sec": time.time() - start_time,
            "avg_cpu_percent": (np.nanmean(cpu_samples) if cpu_samples else None),
            "peak_mem_mb": (np.nanmax(mem_samples) if mem_samples else None),
            "avg_gpu_util_percent": (np.nanmean(gpu_samples) if gpu_samples else None),
            "peak_gpu_mem_mb": (np.nanmax(gpu_mem_samples) if gpu_mem_samples else None),
            "return_code": return_code,
            "error": f"Sampling failed: {e}",
            "traceback": traceback.format_exc(),
            "output_file": output_file
        }

    f.close()

    return {
        "runtime_sec": time.time() - start_time,
        "avg_cpu_percent": (np.nanmean(cpu_samples) if cpu_samples else None),
        "peak_mem_mb": (np.nanmax(mem_samples) if mem_samples else None),
        "avg_gpu_util_percent": (np.nanmean(gpu_samples) if gpu_samples else None),
        "peak_gpu_mem_mb": (np.nanmax(gpu_mem_samples) if gpu_mem_samples else None),
        "return_code": return_code,
        "error": None,
        "traceback": None,
        "output_file": output_file
    }

def taxon_summary(
    df,
    rank=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
    update=False,
    new_ndf=False
):
    """
    Generate a taxonomic summary from a DataFrame with NCBI TaxIDs, including query and assembly counts
    grouped by the specified taxonomic rank.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame (usually a ndf) containing at least a 'taxid' column (NCBI Taxonomy ID),
        a 'query' column (1 or 0 indicating query status), and an 'assembly' column (assembly accession).
        
    rank : list of str, optional
        List of NCBI taxonomic ranks to extract. The final rank in the list (e.g. 'species') will be used
        to group query and assembly counts. Defaults to the standard 7-rank Linnaean hierarchy:
        ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'].
        
    update : bool, optional
	If True, updates the local NCBI taxonomy database via ETE3 before fetching lineages.
	This is recommended the first time you run the function or if your taxonomy data may be outdated.

    new_ndf : bool, optional
        If True, returns the original DataFrame with taxonomic columns (from 'rank') merged in.
        Not recommended to run this function on the new_ndf created, usually results in error.

    Returns
    -------
    a : pandas.DataFrame
        The original DataFrame with taxonomic columns (from `rank`) merged in.

    taxon_df : pandas.DataFrame
        A summary DataFrame, one row per unique taxon at the final rank level,
        with columns:
        - 'taxid' and taxonomic ranks (as defined in `rank`)
        - 'query_tax_count': Number of query sequences mapped to each taxon
        - 'query_tax_pct': Percentage of queries mapped to each taxon
        - 'tax_assembly_count': Number of unique assemblies mapped to each taxon
        - 'tax_assembly_pct': Percentage of assemblies mapped to each taxon
    """
    
    tc = ncbi.TaxonomyCursor() 
    
    # Update database so we can use the most up-to-date taxonomy
    if update:
        tc.cursors['ete3'].update_database()
    
    a = df.copy() 

    # Fetching all taxids
    ndftax = tc.fetchall(a.taxid.drop_duplicates().tolist()) 
    
    # Creating a table to display the requested taxonomies
    z = pd.Series(tc.cursors['ete3'].ete3.get_lineage_translator(ndftax.taxid.drop_duplicates().tolist()), name='lineage').reset_index().rename({'index':'taxid'}, axis=1).explode('lineage')
    z['rank'] = z.lineage.map(tc.cursors['ete3'].ete3.get_rank(z.lineage.drop_duplicates().tolist()))
    z['taxon'] = z.lineage.map(tc.cursors['ete3'].ete3.get_taxid_translator(z.lineage.drop_duplicates().tolist()))
    taxon_df = z[z['rank'].isin(rank)].pivot(index='taxid', columns='rank', values='taxon').reset_index()
    
    # Merging the summary to ndf
    a = a.merge(taxon_df, on='taxid', how='left')
    
    # Measuring the amount of queries that belong to each taxa
    query_df = a[a['query'] == 1]
    query_counts = (query_df.groupby(rank[-1]).size().rename('query_tax_count').reset_index())
    total_queries = query_counts['query_tax_count'].sum()
    query_counts['query_tax_pct'] = (query_counts['query_tax_count'] / total_queries) * 100

    # Measuring the amount of assemblies that belong to each taxa
    assembly_df = a[['assembly', rank[-1]]].drop_duplicates()
    assembly_counts = (assembly_df.groupby(rank[-1]).size().rename('tax_assembly_count').reset_index())
    total_assemblies = assembly_counts['tax_assembly_count'].sum()
    assembly_counts['tax_assembly_pct'] = (assembly_counts['tax_assembly_count'] / total_assemblies) * 100

    taxon_df = taxon_df.merge(query_counts, on=rank[-1], how='left')
    taxon_df = taxon_df.merge(assembly_counts, on=rank[-1], how='left')

    taxon_df[['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']] = (taxon_df[['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']].fillna(0))
    taxon_df['query_tax_count'] = taxon_df['query_tax_count'].astype(int)
    taxon_df['tax_assembly_count'] = taxon_df['tax_assembly_count'].astype(int)    
    
    taxon_df = taxon_df.drop_duplicates(subset=[rank[-1]])
    columns = ['taxid'] + rank + ['query_tax_count', 'query_tax_pct', 'tax_assembly_count', 'tax_assembly_pct']
    taxon_df = taxon_df[columns]
    taxon_df = taxon_df.sort_values(by=['tax_assembly_pct', 'query_tax_pct'], ascending=[False, False])
    
    if new_ndf:
        return a, taxon_df
    else:
        return taxon_df

def shannon(self, ignore_gaps=True):
    """
    Computes Shannon entropy for each column in the alignment.

    Parameters
    ----------
    ignore_gaps : bool
        Whether to exclude gaps ('-' or '.') from entropy calculation.

    Notes
    -----
    You need to run the following commands in order for this method to work properly:

    from rotifer.devel.beta import sequence as rdbs
    rdbs.sequence.compute_shannon_entropy = compute_shannon_entropy

    Returns
    -------
    entropy : pd.Series
        A pandas Series with entropy values indexed by column number.
    """

    # Make sure we only work with sequences
    sequences = self.df[self.df['type'] == 'sequence']['sequence'].tolist()

    if not sequences:
        raise ValueError("No sequences found to compute entropy.")

    # Transpose the alignment into columns
    alignment_length = len(sequences[0])
    columns = zip(*sequences)  # Converts list of sequences to columns

    entropy_values = []
    for col in columns:
        if ignore_gaps:
            col = [res for res in col if res not in ['-', '.']]
        freqs = Counter(col)
        total = sum(freqs.values())
        if total == 0:
            entropy = 0
        else:
            probs = [count / total for count in freqs.values()]
            entropy = -sum(p * np.log2(p) for p in probs if p > 0)
        entropy_values.append(entropy)

    return pd.Series(entropy_values, name='shannon_entropy')

def flag_best_id(group):
    """
    Flag one "best" id per group.

    Priority: pick a random row with id_type 'RefSeq', else 'EMBL-CDS'.
    Returns the same group DataFrame with a new integer 'flag' column (0/1).
    Example:
        radsamorg = pd.concat([flag_best_id(g) for _, g in radsamorg.groupby('qid')])
    """
    n = len(group)
    flags = np.zeros(n, dtype=int)

    id_types = group['id_type'].to_numpy()
    for t in ("RefSeq", "EMBL-CDS"):
        pos = np.flatnonzero(id_types == t)   # positions within the group (0..n-1)
        if pos.size:
            choice = np.random.choice(pos)
            flags[choice] = 1
            break

    group = group.copy()     # avoid SettingWithCopyWarning
    group['flag'] = flags
    return group

def create_attable(ids, update=False, ndf=None, ipg=None):
    """
    Create an enriched sequence object with attribute table, neighborhoods, and Pfam architecture.

    This function initializes a sequence object, fetches genome neighborhoods, IPG information,
    and taxonomy, and compiles them into a unified attribute table (.df). Neighborhoods are
    mapped only for query proteins, and Pfam domains are integrated from HMM scans.

    Parameters
    ----------
    ids : list
        Accession IDs to build the attribute table from.
    update : bool, optional
        If True, updates the ETE3 taxonomy database before fetching.
    ndf : pd.DataFrame, optional
        Precomputed genome neighborhood DataFrame. If given, skips fetching.
    ipg : pd.DataFrame, optional
        Precomputed IPG information. If given, skips fetching.

    Returns
    -------
    a : rotifer.devel.beta.sequence
        Sequence object with:
        - a.df : main attribute table with assembly, taxonomy ranks, neighborhood info, and Pfam architecture
        - a.ndf : genome neighborhood DataFrame
        - a.ipg : IPG information
        - a.tax : taxonomy table (pivoted)
        - a.hsdf : Pfam HMM scan results
    """

    # Step 1: Initialize sequence object
    print("Step 1: Initializing sequence object...")
    a = rdbs.sequence(ids)
    print("Step 1 done.")

    # Step 2: Fetch IPG info if not provided
    if ipg is None:
        print("Step 2: Fetching IPG info...")
        a.ipg = ic.fetchall(ids)
        print("Step 2 done.")
    else:
        a.ipg = ipg
        print("Step 2 skipped (using provided IPG).")

    # Step 3: Fetch genome neighborhood if not provided
    if ndf is None:
        print("Step 3: Fetching genome neighborhood...")
        a.ndf = rdae.add_arch_to_df(gnc.fetchall(ids, ipgs=a.ipg))
        print("Step 3 done.")
    else:
        a.ndf = ndf
        print("Step 3 skipped (using provided NDF).")

    # Step 4: Fetch taxonomy info
    print("Step 4: Fetching taxonomy info...")
    if update:
        print("Updating ETE3 taxonomy database...")
        tc.cursors['ete3'].update_database()
    unique_taxids = a.ndf['taxid'].drop_duplicates()
    tax = tc.fetchall(unique_taxids)
    print("Step 4 done.")

    # Step 5: Expand lineages
    print("Step 5: Expanding lineages...")
    lineage_dict = tc.cursors['ete3'].ete3.get_lineage_translator(unique_taxids.tolist())
    lineage_series = pd.Series(lineage_dict).explode()
    z = lineage_series.reset_index().rename(columns={'index': 'taxid', 0: 'lineage'})
    rank_dict = tc.cursors['ete3'].ete3.get_rank(z['lineage'].tolist())
    taxon_dict = tc.cursors['ete3'].ete3.get_taxid_translator(z['lineage'].tolist())
    z['rank'] = z['lineage'].map(rank_dict)
    z['taxon'] = z['lineage'].map(taxon_dict)
    print("Step 5 done.")

    # Step 6: Pivot taxonomy table
    print("Step 6: Pivoting taxonomy table...")
    ranks = ['domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    taxon_df = (
        z[z['rank'].isin(ranks)]
        .pivot(index='taxid', columns='rank', values='taxon')
        .reset_index()
    )
    print("Step 6 done.")

    # Step 7: Merge taxonomy to main df
    print("Step 7: Mapping taxonomy and assembly to .df...")
    tmp_tax = a.ndf[['pid', 'taxid', 'assembly']].drop_duplicates()
    tmp_tax = tmp_tax.merge(taxon_df, on='taxid', how='left')
    a.df = a.df.merge(tmp_tax, left_on='id', right_on='pid', how='left', suffixes=('', '_y'))
    a.df.drop(columns=['pid'], inplace=True, errors='ignore')
    a.tax = taxon_df
    print("Step 7 done.")

    # Step 8: Map neighborhoods only for query proteins
    print("Step 8: Mapping neighborhood and ss_neighborhood...")
    cndf = a.ndf.compact().reset_index()          # default compact
    cndf_ss = a.ndf.compact(strand=True).reset_index()  # strand-specific

    # Map block_id to .df only for query proteins
    a.df['block_id'] = a.df.id.map(dict(zip(a.ndf[a.ndf['query'] == 1].pid, a.ndf[a.ndf['query'] == 1].block_id)))

    a.df['neighborhood'] = a.df.block_id.map(dict(zip(cndf.block_id, cndf.compact)))
    a.df['ss_neighborhood'] = a.df.block_id.map(dict(zip(cndf_ss.block_id, cndf_ss.compact)))

    # Step 9: Fill missing values
    a.df.fillna('Missing', inplace=True)
    a.tax.fillna('Missing', inplace=True)

    print("Step 8-9 done. Attribute table ready.")

    # Step 10: Pfam integration
    hsdf = rdae.hmmscan(a.df.id.to_list())
    a.hsdf = hsdf
    hsdff = riu.filter_nonoverlapping_regions(hsdf, **riu.config['hmmer'])
    idx = hsdf.groupby('sequence')['evalue'].idxmin()
    best_hits = hsdf.loc[idx, ['sequence', 'model']].set_index('sequence')['model']
    a.df['arch'] = a.df['id'].map(best_hits).fillna('Missing')

    print("Step 10 done. Pfam column added.")

    return a

class PhyloBuilder:
    """
    Build and hold a simple phylogeny (NJ or UPGMA) from a rotifer.sequence object
    (or from a Biopython MultipleSeqAlignment / list of SeqRecords).

    Key features:
    - Automatically builds the tree when an aligned rdbs.sequence is provided.
    - If a rdbs.sequence is not aligned, it will call `.align()` (non-inplace) to align.
    - Stores the alignment (msa), distance matrix (Bio object and pd.DataFrame),
      the tree (Bio.Phylo tree), and metadata.
    - Visualization via Bio.Phylo.draw and export to Newick.

    Parameters
    ----------
    seqs : optional
        A rotifer.devel.beta.sequence object (preferred), or a Bio.Align.MultipleSeqAlignment,
        or a list of Bio.SeqRecord, or a list/Series of sequence strings (ids must be provided).
    method : int, optional
        1 => UPGMA, 2 => Neighbor-Joining (default 2).
    matrix : str, optional
        Substitution matrix name accepted by Bio.Phylo.TreeConstruction.DistanceCalculator.
        Default 'blosum62'.
    metadata : dict, optional
        Arbitrary metadata to store with the object.
    build_on_init : bool, optional
        If True and seqs is provided, automatically run build() on init (default True).
    """
    def __init__(
        self,
        seqs: Optional[Union['rdbs.sequence', MultipleSeqAlignment, List[SeqRecord], List[str]]] = None,
        method: int = 2,
        matrix: str = 'blosum62',
        metadata: Optional[Dict] = None,
        build_on_init: bool = True,
    ):
        self.seqs_input = seqs
        self.method = int(method)
        self.matrix = matrix
        self.metadata = dict(metadata) if metadata is not None else {}
        self.msa: Optional[MultipleSeqAlignment] = None
        self.distance_matrix = None  # Bio DistanceMatrix object
        self.dm_df: Optional[pd.DataFrame] = None
        self.tree = None  # Bio.Phylo tree
        # If user passed a rotifer sequence object (rdbs.sequence), keep reference
        self._rseq = None

        if seqs is not None and build_on_init:
            self.build()

    def _to_msa_from_rdbs(self, rseq) -> MultipleSeqAlignment:
        """
        Convert a rotifer.sequence object to a Biopython MultipleSeqAlignment.
        If the rseq is not aligned, call its .align() (non-inplace) to obtain an aligned object.
        """
        # Prefer using rseq.align() if sequences are not same length (i.e. not aligned)
        seq_df = rseq.df.query('type == "sequence"').copy()
        # check if aligned
        if seq_df.sequence.str.len().nunique() != 1:
            # call align non-inplace to avoid mutating input object
            try:
                aligned = rseq.align(inplace=False)
            except Exception as e:
                raise RuntimeError(f"Failed to align sequence object automatically: {e}")
            seq_df = aligned.df.query('type == "sequence"').copy()
            rseq = aligned

        # Build SeqRecords
        seq_records = []
        for _, row in seq_df.iterrows():
            seq_records.append(SeqRecord(Seq(row['sequence']), id=str(row['id']), description=''))

        msa = MultipleSeqAlignment(seq_records)
        self._rseq = rseq  # store reference to possibly aligned sequence object
        return msa

    def _to_msa_from_generic(self, seqs) -> MultipleSeqAlignment:
        """
        Convert a generic input to MultipleSeqAlignment:
        - If already MultipleSeqAlignment: return as is
        - If list of SeqRecord: wrap
        - If list/Series of strings: require that it's paired with a list of ids (not handled here)
        """
        if isinstance(seqs, MultipleSeqAlignment):
            return seqs
        # list of SeqRecord
        if isinstance(seqs, list) and all(isinstance(x, SeqRecord) for x in seqs):
            return MultipleSeqAlignment(seqs)
        # list/Series of strings -> require ids embedded? try to use indices as ids
        if isinstance(seqs, (list, pd.Series)) and all(isinstance(x, str) for x in seqs):
            seq_records = [SeqRecord(Seq(s), id=f"seq{i}") for i, s in enumerate(seqs)]
            return MultipleSeqAlignment(seq_records)

        raise TypeError("Unsupported seqs type for conversion to MultipleSeqAlignment.")

    def build(self, overwrite: bool = True):
        """
        Build the distance matrix and tree from the provided sequences.

        overwrite: if False, will not rebuild if a tree already exists.
        """
        if (self.tree is not None) and (not overwrite):
            return

        # Prepare MSA
        seqs = self.seqs_input
        if seqs is None:
            raise ValueError("No sequences provided to build the phylogeny.")

        if hasattr(seqs, '__class__') and seqs.__class__.__name__ == 'sequence':
            # assume rotifer.devel.beta.sequence class instance
            msa = self._to_msa_from_rdbs(seqs)
        else:
            msa = self._to_msa_from_generic(seqs)

        self.msa = msa

        # Calculate distances
        try:
            calculator = DistanceCalculator(self.matrix)
        except Exception as e:
            raise ValueError(f"DistanceCalculator could not be initialized with model '{self.matrix}': {e}")

        # get_distance expects a MultipleSeqAlignment
        dm = calculator.get_distance(msa)
        self.distance_matrix = dm

        # keep a pandas copy for convenience
        names = dm.names
        mat = dm.matrix
        # dm.matrix is a lower-triangular list of lists (Bio's DistanceMatrix), convert to full DataFrame
        df = pd.DataFrame(0.0, index=names, columns=names, dtype=float)
        for i, n in enumerate(names):
            row = mat[i]
            # row length equals i+1
            for j, val in enumerate(row):
                df.iat[i, j] = val
                df.iat[j, i] = val
        self.dm_df = df

        # Build tree
        constructor = DistanceTreeConstructor()
        if self.method == 1:
            self.tree = constructor.upgma(dm)
        else:
            self.tree = constructor.nj(dm)

        # store some metadata
        self.metadata.setdefault('method', 'upgma' if self.method == 1 else 'nj')
        self.metadata.setdefault('matrix', self.matrix)
        self.metadata.setdefault('nseq', len(self.msa))
        return self

    def draw(self, show: bool = True, **phylo_draw_kwargs):
        """
        Draw the current tree using Bio.Phylo.draw.
        Returns the matplotlib Axes if available.
        """
        if self.tree is None:
            raise RuntimeError("No tree built yet. Call build() first.")
        ax = Phylo.draw(self.tree, **phylo_draw_kwargs)
        if show:
            try:
                import matplotlib.pyplot as plt
                plt.show()
            except Exception:
                pass
        return ax

    def to_newick(self) -> str:
        """Return tree as Newick string (or empty string if no tree)."""
        if self.tree is None:
            raise RuntimeError("No tree built yet. Call build() first.")
        buf = StringIO()
        Phylo.write(self.tree, buf, "newick")
        return buf.getvalue().strip()

    def export(self, path: str, fmt: str = "newick"):
        """
        Export tree to file in a supported format (newick, nexus, phyloxml, etc.).
        """
        if self.tree is None:
            raise RuntimeError("No tree built yet. Call build() first.")
        Phylo.write(self.tree, path, fmt)

    def get_distance_dataframe(self) -> pd.DataFrame:
        """
        Return the distance matrix as a pandas DataFrame.
        """
        if self.dm_df is None:
            raise RuntimeError("Distance matrix not available; call build() first.")
        return self.dm_df.copy()

    def attach_metadata(self, d: Dict):
        """Merge provided metadata dict into builder metadata."""
        if not isinstance(d, dict):
            raise TypeError("metadata must be a dict")
        self.metadata.update(d)

    def __repr__(self):
        info = f"<PhyloBuilder nseq={len(self.msa) if self.msa is not None else 0} method={'UPGMA' if self.method==1 else 'NJ'} matrix={self.matrix}>"
        return info

def evalue_hist(df, step=20):
    """
    Automatically creates logarithmic bins for E-values.
    'step' defines how many orders of magnitude to group (e.g., 20 = 1e-80 to 1e-60).
    """
    # We find the min/max of non-zero values to define our scale
    non_zero = df[df['evalue'] > 0]['evalue']
    
    if non_zero.empty:
        # Fallback if all values are 0
        return pd.DataFrame({'evalue range': ['0'], 'unique_proteins': [df['sequence'].nunique()]})

    # 2. Find the powers of 10 for the scale
    min_exp = int(np.floor(np.log10(non_zero.min())))
    max_exp = int(np.ceil(np.log10(non_zero.max())))
    
    # 3. Create the bin edges (powers of 10)
    # We create a range of exponents from min to max
    bin_exps = np.arange(min_exp, max_exp + step, step)
    bins = [0] + list(np.power(10.0, bin_exps.astype(float)))
    
    # 4. Generate the labels automatically
    labels = [f"0 to 1e{bin_exps[0]}"]
    for i in range(len(bin_exps) - 1):
        labels.append(f"1e{bin_exps[i]} to 1e{bin_exps[i+1]}")
    
    # 5. Bin the data and count unique sequences
    df_temp = df.copy()
    df_temp['range'] = pd.cut(df_temp['evalue'], bins=bins, labels=labels, include_lowest=True)
    
    summary = df_temp.groupby('range')['sequence'].nunique().reset_index()
    summary.columns = ['evalue range', 'unique_proteins']
    
    # Return only rows with data
    return print(summary[summary['unique_proteins'] > 0])

def operon_fig(df, group_col='block_id', label_col='pfam', org_col='organism',
                       output_file='operon_fig_out.svg', max_colors=5, 
                       highlight_query=True, query_same_direction=False, 
                       font_size=10, ignore_domains=None):
    
    # Define generic/structural domains to ignore when coloring
    if ignore_domains is None:
        ignore_domains = ['TM', 'SP', 'LP', 'LIPO', 'SIG']

    a = df.copy()
    
    # --- Safe Column Mapping ---
    a['ID'] = a.get(group_col, 'Unknown_Block')
    a['org_name'] = a.get(org_col, 'Unknown Organism')
    
    # Safely handle missing fallback columns to prevent .fillna(None) crashes
    if label_col == 'pfam':
        a['domain'] = a['pfam'] if 'pfam' in a.columns else pd.Series([np.nan]*len(a))
        if 'aravind' in a.columns:
            a['domain'] = a['domain'].fillna(a['aravind'])
        if 'product' in a.columns:
            a['domain'] = a['domain'].fillna(a['product'])
        a['domain'] = a['domain'].fillna('unk')
    else:
        a['domain'] = a[label_col].fillna('unk') if label_col in a.columns else 'unk'

    # Parse the binary query flag
    if 'query' in a.columns:
        a['is_query'] = (a['query'] == 1) | (a['query'] == '1') | (a['query'] == True)
    else:
        a['is_query'] = False

    a['pid_order'] = pd.factorize(a['ID'])[0]
    
    # --- Filter by Query Direction (Optional) ---
    if query_same_direction:
        filtered_dfs = []
        for pid, group in a.groupby('pid_order'):
            query_rows = group[group['is_query']]
            if not query_rows.empty:
                # Get the strand of the query
                q_strand = query_rows.iloc[0].get('strand')
                if pd.notna(q_strand):
                    # Keep only rows matching the query strand
                    group = group[group['strand'] == q_strand]
            filtered_dfs.append(group)
        a = pd.concat(filtered_dfs)

    a = a.reset_index(drop=True).reset_index()

    # --- Priority Coloring Logic ---
    ignore_lower = [x.lower() for x in ignore_domains] + ['unk', ' ', '-', '?']
    
    # Identify domains that should NOT get a color
    is_ignorable = (
        a['domain'].astype(str).str.lower().isin(ignore_lower) | 
        a['domain'].astype(str).str.lower().str.contains('hypothetical', na=False)
    )
    
    # 1. ALWAYS grab the unique query domains first
    query_domains = a[a['is_query'] & ~is_ignorable]['domain'].unique().tolist()
    
    # 2. Grab the most frequent non-query domains
    remaining_slots = max(0, max_colors - len(query_domains))
    freq_domains = a[~a['domain'].isin(query_domains) & ~is_ignorable]['domain'].value_counts().head(remaining_slots).index.tolist()
    
    # Combine and map to a pastel palette
    final_color_list = query_domains + freq_domains
    palette = sns.color_palette("pastel", len(final_color_list)).as_hex()
    color_map = {domain: color for domain, color in zip(final_color_list, palette)}

    # --- Alignment Preparation ---
    max_id_len = a['ID'].astype(str).str.len().max()
    max_org_len = a['org_name'].astype(str).str.len().max()
    query_pids = a[a['is_query']]['pid'].astype(str)
    max_q_len = query_pids.str.len().max() if not query_pids.empty else 8
    max_width = max(max_id_len, max_org_len, max_q_len)

    def pad_and_escape(text):
        padded_text = str(text).ljust(max_width)
        return html.escape(padded_text).replace(" ", "&nbsp;")

    # --- Graph Initialization ---
    A = pgv.AGraph(directed=True)
    A.graph_attr.update(nodesep=0.05, ranksep=0.15)
    
    first_nodes_per_row = []

    for pid, group in a.groupby('pid_order'):
        block_id_str = pad_and_escape(group['ID'].iloc[0])
        org_str = pad_and_escape(group['org_name'].iloc[0])
        
        query_rows = group[group['is_query']]
        query_pid_raw = query_rows['pid'].iloc[0] if not query_rows.empty else "No Query"
        query_pid_str = pad_and_escape(query_pid_raw)

        # --- Draw HTML Label Node ---
        label_node_id = f"label_{pid}"
        html_label = (
            f'<<TABLE BORDER="0" CELLBORDER="0" CELLPADDING="0" CELLSPACING="0">'
            f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}"><B>{query_pid_str}</B></FONT></TD></TR>'
            f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas" POINT-SIZE="{font_size}">{block_id_str}</FONT></TD></TR>'
            f'<TR><TD ALIGN="LEFT"><FONT FACE="Consolas italic" POINT-SIZE="{font_size}">{org_str}</FONT></TD></TR>'
            f'</TABLE>>'
        )
        
        A.add_node(label_node_id, label=html_label, shape='none', margin=0.1)
        node_indices = [label_node_id]
        
        # --- Draw Genes ---
        for orig_idx, row in group.iterrows():
            is_target = highlight_query and row['is_query']
            border_color = 'red' if is_target else 'black'
            border_width = '3' if is_target else '1'
            
            strand_val = row.get('strand', 1)
            node_shape = 'rarrow' if strand_val == 1 else 'larrow' if strand_val == -1 else 'box'
            
            bg_color = color_map.get(row['domain'], '#ffffff')
            node_id = f"node_{orig_idx}"
                
            A.add_node(
                node_id, 
                label=str(row['domain']), 
                shape=node_shape,
                style='filled', 
                fixedsize='false', 
                margin='0.1,0.05',  
                height='0.4', 
                color=border_color, 
                penwidth=border_width,
                fillcolor=bg_color,
                fontsize=font_size,
                fontname="Consolas" 
            )
            node_indices.append(node_id)
            
        A.add_subgraph(node_indices, rank="same")
        
        for i in range(len(node_indices) - 1):
            A.add_edge(node_indices[i], node_indices[i+1], style='invis', penwidth=0)
            
        first_nodes_per_row.append(label_node_id)

    # --- Enforce Strict Left Alignment ---
    for i in range(len(first_nodes_per_row) - 1):
        A.add_edge(first_nodes_per_row[i], first_nodes_per_row[i+1], style='invis', penwidth=0, weight=10000)

    A.draw(output_file, prog="dot")
    return a
