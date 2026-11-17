from gc import enable
import rotifer
import rotifer.pipeline
import rotifer.db.ncbi as ncbi
from rotifer.pandas import functions as rpf
from rotifer.core.functions import loadConfig

from copy import deepcopy
from io import StringIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rotifer.core import functions as rcf
import pandas as pd
import numpy as np
import os
import re
logger = rotifer.logging.getLogger(__name__)

# Defaults
config = loadConfig(__name__, defaults = {
    'fetch': 'pfetch',
    'pdb_dir': os.path.join(rotifer.config['data'],"pdb"),
    'databases': ['pdb70','pfam'],
    'databases_path': os.path.join(rotifer.config['data'],"hhsuite"),
    'local_database_path': [ os.path.join(rotifer.config['data'],"fadb","nr","nr") ],
    'html_colors': rcf.findDataFiles(":colors/aminoacids_HEX.yml"), 
})

class sequence(rotifer.pipeline.Annotatable):
    """
    This class represents a multiple sequence alignment (MSA)
    or a collection of unaligned sequences.

    It provides methods for MSA slicing, visualization, alignment
    rebuilding and annotation.

    Alignment annotations
    ---------------------
    Four types of annotations are supported by this class:

    - Row annotations
      Refers to properties of alignment rows, such as cluster
      assignments, functional annotations for entire sequences,
      such as EC numbers, or organism identifiers.

      Sequence annotations are stored as columns in the internal
      Pandas DataFrame and the user may manipulate such columns
      directly via the ``df`` attribute, as in the example below,
      where we remove the ``c80i30`` added by ``add_cluster``:

      >>> aln = aln.add_cluster(identity = 0.3)
      >>> aln.df.drop('c80i30', inplace=True, axis=1)

    - Character-based residue annotations
      Refers to a single column and are specified as strings of the
      same length as the alignment itself. Examples include
      **consensus sequences** and per-residue **secondary structure**
      codes.

      This type of annotation can be added to the internal Pandas
      DataFrame directly, i.e. using Pandas's ``append`` or
      ``insert`` methods but the user must be **careful** to set the
      value or the ``type`` column to something other than
      ``sequence``, which is researved for the aligned sequences.

    - Numerical values assigned to individual residues
      This type of annotation also assigns annotations to single
      columns of the alignment but are best described by numerical
      values, instead of single characters. Examples include the
      frequency of a given residue at each alignment column, the
      number of gaps per column or the column entropies.

    - Alignment or sequence regions
      Describe properties for subsets of columns in the alignment
      or for regions within a specific sequence. The former uses
      column indices to describe alignment features, while the later
      specifies a region within a given sequence using residue 
      coordinates that ignore gaps.

    Specialized methods are provided for the most common annotations
    (see, for example, ``add_cluster``, ``add_consensus``).

    Parameters
    ----------
    input_data : str or file
        Input file name, string or open file handler
    input_format : str
        Input alignment format. See Bio.SeqIO and/or
        Bio.AlignIO for a list of supported formats
    name:
        An (optional) name for the multiple sequence alignment

    See also
    --------
    Bio.SeqIO   : BioPython parser for aligned/unaligned sequences
    Bio.AlignIO : BioPython parser for aligned sequences
    pandas      : library used internally to store sequences

    Examples
    --------
    Simplest usage: visualize a (multi)FASTA file.

    >>> from rotifer.devel.beta.sequence import sequence
    >>> aln = sequence("alignment.aln")
    >>> aln.view()

    Load alignment from another file format, such as
    StockHolm and show the consensus above the alignment:

    >>> from rotifer.devel.beta.sequence import sequence
    >>> sto = sequence("alignment.sto", format="stockholm")
    >>> sto.add_consensus().view()

    Parsing a FASTA alignment from an explicit string:

    >>> from rotifer.devel.beta.sequence import sequence
    >>> aln = sequence(">Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n")
    >>> aln.add_consensus().view()
    """
    def __init__(self,
            input_data=None,
            input_format='fasta',
            name=None,
            local_database_path=config['local_database_path'],
            *args, **kwargs):
        from io import IOBase
        self._reserved_columns = ['id','sequence','length','type','description']
        self.input_format = input_format
        self.name = name
        if not isinstance(local_database_path,list):
            local_database_path = [ local_database_path ]

        # Cursors: use method arguments or defaults
        kwargs['local_database_path'] = local_database_path
        kwargs.update({ k:v for k,v in deepcopy(config).items() if k not in kwargs })
        self.cursors = {
            'neighborhood': ncbi.GeneNeighborhoodCursor(**kwargs),
            'fasta': ncbi.FastaCursor(**kwargs),
        }

        # Generate empty object
        if input_data is None:
            self.df = pd.DataFrame({}, columns=self._reserved_columns)
            self.input_format = None
            self.file_path = None

        # Initialize for strings
        elif isinstance(input_data,str):
            if os.path.exists(input_data):
                self.file_path = input_data
                if not name:
                    self.name = os.path.basename(input_data)
            else:
                self.file_path = 'StringIO'
                input_data = StringIO(input_data)
                if not name:
                    self.name = str(input_data)
            input_data = SeqIO.parse(input_data, input_format)
            self.df = self._seqrecords_to_dataframe(input_data)

        # Initialize for IO streams
        elif isinstance(input_data, IOBase):
            input_data = SeqIO.parse(input_data, input_format)
            self.df = self._seqrecords_to_dataframe(input_data)

        # Initialize for list of Bio.SeqRecords
        elif isinstance(input_data, list):
            if isinstance(input_data[0],SeqRecord):
                self.df = self._seqrecords_to_dataframe(input_data)
            elif isinstance(input_data[0],str):
                seqrecords = self.cursors['fasta'].fetchall(input_data)
                self.df = self._seqrecords_to_dataframe(seqrecords)

        # Initialize for Pandas DataFrame
        elif isinstance(input_data, pd.DataFrame):
            other = [ x for x in input_data.columns if x not in self._reserved_columns ]
            self.df = input_data[['id','sequence']]
            self.df['type'] = 'sequence'
            self.df[other] = input_data[other].copy()

        # Initialize for Pandas Series
        elif isinstance(input_data, pd.Series):
            self.df = input_data.dropna().reset_index()
            self.df.columns = ['id','sequence']
            self.df['type'] = 'sequence'

        # For creating an empty object
        else:
            self.df = pd.DataFrame(columns=self._reserved_columns)

        # Still unnamed? Try using the input filename
        if not name:
            if hasattr(input_data,"stream"):
                if hasattr(input_data.stream,"filename") and callable(input_data.stream.filename):
                    self.name = os.path.basename(input_data.stream.filename())
                elif hasattr(input_data.stream,"name"):
                    self.name = os.path.basename(input_data.stream.name)
            elif not self.df.empty:
                self.name = self.df.id.iloc[0]

        # Make sure the new object is clean!
        self._reset()

    # Functionis to alow item assignment to the sequence class.
    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getitem__(self, key):
        return self.filter(keep=[key])

    def _seqrecords_to_dataframe(self, data):
        cols = self._reserved_columns
        parsed = []
        for s in data:
            current = [ s.id, str(s.seq) ]
            current.append(len(current[1].replace("-","")))
            current.append('sequence')
            current.append(s.description)
            parsed.append(current)
        df = pd.DataFrame(parsed, columns=cols)
        return df

    def _reset(self):
        # Recalculate each sequence's length
        maxlen = self.df.sequence.str.len().max()
        self.df['length'] = np.where(
                self.df.type == "sequence",
                self.df.sequence.str.replace('-', '').str.len(),
                maxlen)
        if not self.name and not self.df.empty:
            self.name = self.df.id.iloc[0]
        self.df.reset_index(drop=True, inplace=True)

    def __len__(self):
        return len(self.df.query('type == "sequence"'))

    def __eq__(self, other):
        if type(self) != type(other):
            logger.error(f'Cannot compare sequence object and {type(other)}')
            return None
        return self.checksum == other.checksum

    @property
    def checksum(self):
        import hashlib
        cksum = self.df.sequence.sort_values().to_list()
        cksum = "\n".join(cksum).encode()
        cksum = hashlib.md5(cksum).hexdigest()
        return cksum

    def copy(self, reset=False):
        '''
        Copy alignment to new object.
        '''
        from copy import deepcopy
        result = deepcopy(self)
        if reset:
            result._reset()
        return result

    def filter(self, query=None, minlen=0, maxlen=0, keep=[], remove=[], regex=None):
        '''
        Function to filter / select rows from the MSA.

        Parameters
        ----------
        query: str
          String to use as input for the ``query`` method of the
          internal Pandas Dataframe
        minlen : int
          Length of the shortest sequence, without gaps
        maxlen : int
          Length of the longest sequence, without gaps
        keep : list of strings
          List of sequence IDs to preserve

          The ``keep`` parameter has the highest precedence over
          any other parameters, i.e. any sequence listed in
          ``keep``, is **never** removed.

        remove : list of strings
          List of sequences IDs to be removed.

          Sequences in this list are removed even if failing the
          criteria set by other parameters, except ``keep``.

        regex : regular expression
          Regular expression to identify rows based on
          residue patterns matching their aminoacid/nucleotide
          sequences.

        Returns
        -------
        New MSA

        Examples
        --------
        Filter sequences with at least 50 and at most 100 residues

        >>> aln.filter(minlen=50, maxlen=100)

        Filter sequences based on a cluster but make sure NPE27555.1
        is also kept

        >>> aln = aln.add_cluster(identity=0.6)
        >>> aln.filter('c80i80 == "WP_091936315.1"', keep="NPE27555.1")
        '''
        result = self.copy(reset=True)

        # Build query statement
        querystr = []
        if query:
            querystr.append(query)
        if minlen > 0:
            querystr.append('(length >= @minlen)')
        if maxlen > 0:
            querystr.append('(length <= @maxlen)')
        if regex:
            querystr.append('(sequence.str.replace("-", "").str.contains(@regex, regex=True))')
        if remove:
            querystr.append('(id not in @remove)')
        querystr = " and ".join(querystr)
        if len(keep) > 0:
            if querystr:
                querystr = f'(id in @keep) or ({querystr})'
            else:
                querystr = f'id in @keep'
        if querystr:
            result.df = result.df.query(querystr)
        result.df.reset_index(drop=True, inplace=True)
        return result

    def get_alignment_length(self):
        '''
        Retrieve the number of columns in the alignment.

        Returns
        -------
        Integer

        Examples
        --------
        >>> al = aln.get_alignment_length()
        '''
        if self.df.empty:
            return 0
        else:
            return int(self.df.query('type == "sequence"').sequence.str.len().max())

    def get_group(self, group: dict):
        gname = tuple(group.values())
        if len(gname) == 1:
            gname = gname[0]
        aln = self.copy()
        aln.df = self.df.groupby(list(group.keys())).get_group(gname)
        aln._reset()

    def sample(self, n: 'int | None' = None, frac: 'float | None' = None, replace = False):
        other = self.filter('type != "sequence"').copy()
        seqs =  self.filter('type == "sequence"').copy()
        seqs.df = seqs.df.sample(n=n, frac=frac, replace=replace)
        return concat([ other, seqs ])

    def position_to_column(self, position, reference):
        '''
        Map sequence coordinates to alignment columns

        Parameters
        ----------
        position: Pandas Series, list, tuple or integer
          The user may provide one or several positions.

        reference: string
          Identifier of the sequence that must be used as reference
          for the sequence-based coordinate system.

        Returns
        -------
          List of integers

        Coordinate systems
        ------------------
        This method maps from a reference sequence-based coordinate
        system to aligment-based coordinates, i.e. column positions.

        * Alignment-based coordinates
          Refer to the column indices in the original alignment

        * Sequence-based coordinates
          Refer to residues in a reference sequence, without gaps

        All types of coordinates systems above are **one-based**
        (residue-based) and closed on both ends, i.e. the first
        column or residue position is 1 and the last column equals
        the length of the alignment or sequence.

        Examples
        --------

        - Locate the columns corresponding to residues 10 to 80 of
          the sequence WP_003247817.1

          >>> aln.position_to_column((10,80),"WP_003247817.1"))

        '''
        if isinstance(position, int):
            position = [position]
        refseq = self.df.query(f'''id == "{reference}"''')
        if refseq.empty:
            logger.error(f"Reference sequence {reference} could not befound in this alignment.")
            return None
        refseq = pd.Series(list(refseq.sequence.values[0]))
        refseq.index = refseq.index + 1 # Adjust alignment coordinates to interval [1,length]
        refseq = refseq.where(lambda x: x != '-').dropna().reset_index().rename({'index':'mapped_position'}, axis=1)
        refseq.index = refseq.index + 1 # Adjust reference sequence coordinates
        missing = set(position) - set(refseq.index)
        if missing:
            logger.error(f"The following positions could not befound in {reference}: {missing}")
            return None
        else:
            return refseq.loc[position].mapped_position.tolist()

    def slice(self, position):
        '''
        Select and concatenate one or more sets of columns.

        Parameters
        ----------
        position: (list of) tuple of two integers

        Coordinate systems
        ------------------
        This method cuts the input aligment based on two types of
        coodinate systems:

        * Alignment-based coordinates
          Refer to the column indices in the original alignment

        * Sequence-based coordinates
          Refer to residues in a reference sequence, without gaps

        All types of coordinates systems above are **one-based**
        (residue-based) and closed on both ends, i.e. the first
        column or residue position is 1 and the last column equals
        the length of the alignment or sequence.

        Notice that coordinates must be provided as tuples (see below).

        Examples
        --------

        - Alignment-based coordinates
          Fetch columns 20 to 130:

          >>> aln.slice((20,130))

        - Sequence-based coordinates
          Fetch columns corresponding to residues 10 to 80 of
          the sequence WP_003247817.1

          >>> aln.slice((10,80,"WP_003247817.1"))

        - You can also fetch several non-overlapping regions
          and these will be concatenated after removal of the
          interviening regions:

          >>> aln.slice([ (10,80),(110,160) ])

        '''
        result = self.copy()
        if isinstance(position[0],int):
            position = [position]

        # Cut slices
        sequence  = []
        for pos in position:
            pos = [*pos]
            if len(pos) == 3:
                pos[0:2] = self.position_to_column(pos[0:2],reference=pos[2])
            sequence.append(result.df.sequence.str.slice(pos[0]-1, pos[1]))

        # Concatenate slices per row
        sequence = pd.concat(sequence, axis=1, ignore_index=True)
        sequence = sequence.sum(axis=1).tolist()
        result.df['sequence'] = sequence

        # Return new sequence object
        result._reset()
        return result

    def sort(self, by=['length'], ascending=True, inplace=False, id_list=None, tree_file=None, tree_format='newick'):
        """
        Sort alignment rows.

        Parameters
        ----------
        by : list of strings, default ['length']
            List of criteria for sorting.

            The following criteria are always supported:
            - id : sort by sequence identifiers
            - length : sequence length
            - list : order follows a list of sequence IDs
                See parameter ``id_list``.
            - name : sort by sequence identifiers
            - tree : sort by phylogeny
                See parameter ``tree_file``.
                Leaf names must match sequence identifiers.
            - residues: alignment columns

            Additionally, user-supplied sequence annotations may
            also be used to sort alignment rows. See the section
            ```Alignment annotations```.

        ascending : bool or list of bools, default True
            If True, sort in ascending order, otherwise descending.

        inplace : bool, default False
            Sort in place, i.e. without copying of the object

        id_list : list of strings or file name
            A list of sequence identifiers to use in sorting.
            If a file name is given as a string, the list will be
            automatically loaded from that file.

        tree_file : path to phylogeny
            For this parameter to be effective, at least a subset
            of the leaves in the phylogeny must match sequence
            identifiers in the alignment.

            Leaf names that do not match sequence identifiers are
            ignored. Sequence identifiers from the alignment that
            are absent in the phylogeny are listed after

        tree_format: phylogeny file format
            All formats supported by Ete3 or BioPython are valid

        Returns
        -------
        MSA with rows in a new order

        See also
        --------
        sequence.annotations : 

        Examples
        --------
        Default: sort by sequence length, ignoring gaps,
        from shortest to longest:

        >>> aln2 = aln.sort()

        Sort by length and then by the sequence cluster
        (see ``add_cluster``).

        >>> aln = aln.add_cluster(coverage=0.8, identity=0.3)
        >>> aln = aln.sort(['length','c80i30'])

        Make sure sequence ``PUA33204.1`` shows up first.
        Sort the others by length.

        >>> aln.sort(['list','length'], id_list=['PUA33204.1']).view()

        Sort in-place by columns 32 and 27 and by length:

        >>> aln.sort([32,27,'length'], inplace=True)

        """
        import sys

        # Inplace mode
        if inplace:
            result = self
        else:
            result = self.copy()
        result.df.reset_index(inplace=True, drop=True)

        # Prepare local variables
        cols = result.df.columns
        isnumeric = [ x for x in by if x not in cols and isinstance(x,int) ]
        if isnumeric:
            matrix = self.residues

        fields = []
        identifiers = result.df.id
        for item in by:
            if item == 'tree':
                from ete3 import Tree
                tree = Tree(tree_file)
                leaves = pd.Series([ x.name.strip("'") for x in tree.get_leaves() ], name='leaves')
                leaves = leaves.reset_index().rename({'index':'order'}, axis=1).set_index('leaves').order.to_dict()
                fields.append(identifiers.where(lambda x: x.isin(leaves)).replace(leaves).rename("tree"))

            elif item == 'list':
                if isinstance(id_list,str) and os.path.exists(id_list):
                    ids = pd.read_csv(id_list, names=['id']).id
                elif isinstance(id_list, list):
                    ids = pd.Series(id_list, name='id')
                ids = ids.reset_index().rename({'index':'order'}, axis=1).set_index('id').order.to_dict()
                fields.append(identifiers.where(lambda x: x.isin(ids)).replace(ids).rename("list"))

            elif item == 'name':
                fields.append(identifiers.rename("name"))

            elif item in cols:
                fields.append(result.df[item])

            elif isinstance(item,int):
                fields.append(matrix.loc[:,item])

            else:
                logger.error(f'sequence.sort: Unsupported criteria or missing annotation ({item})')

        # Concatenate sorting fields
        fields = pd.concat(fields, axis=1).sort_values(by=by, ascending=ascending)
        #result.df = result.df.reindex(fields.index)
        result.df = result.df.loc[fields.index]
        if not inplace:
            return result

    def add_cluster(self, coverage=0.8, identity=0.7, cascade=None, inplace=False):
        '''
        Add or update MMseqs2 clustering data for all sequences.

        Parameters
        ----------
        coverage : float
            Minimum coverage required for both sequences in each
            pairwise alignment.
        identity : float
            Minimum percentage identity required for both sequences
            in each pairwise alignment.
        cascade: (list of) tuples, default: None
            Generate multiple hierarchically related clusterings
            through cascading. The arguments must be pairs of 
            coverage and identity values.

            Note:
             When using cascade, coverage and identity are ignored
        inplace: bool
            Add column inplace, without creating a copy of the MSA
            If set to True, this method return nothing.

        Returns
        -------
        A new MSA object

        Examples
        --------
        Add an annotation column listing representatives of clusters
        defined by 60% alignment coverage and 30% identity for both
        sequences

        >>> aln.add_cluster(coverage=0.6, identity=0.3)
        '''
        import tempfile
        from subprocess import Popen, PIPE, STDOUT

        # Adjust parameters
        if not cascade:
            cascade = [(coverage, identity)]
        elif not isinstance(cascade,list):
            cascade = [cascade]

        if inplace:
            result = self
        else:
            result = self.copy()

        path = os.getcwd()
        lastname = None
        for coverage, identity in cascade:
            if identity > 1:
                identity /= 100
            if coverage > 1:
                coverage /= 100
            name = f'c{int(float(coverage)*100)}i{int(float(identity)*100)}'

            with tempfile.TemporaryDirectory() as tmpdirname:
                os.chdir(tmpdirname)
                if lastname:
                    result.filter(f'id == {lastname}').to_file(f'{tmpdirname}/seqaln', remove_gaps=True)
                else:
                    result.to_file(f'{tmpdirname}/seqaln', remove_gaps=True)
                Popen(f'mmseqs easy-cluster {tmpdirname}/seqaln nr --min-seq-id {identity} -c {coverage} tmp', stdout=PIPE,shell=True).communicate()
                d = pd.read_csv(f'{tmpdirname}/nr_cluster.tsv', sep="\t", names=['cluster', 'pid'])
                d = d.set_index('pid').cluster.to_dict()
                if lastname:
                    result.df[name] = result.df[lastname].replace(d)
                else:
                    result.df[name] = result.df.id.replace(d)
            os.chdir(path)
            lastname = name

        if inplace:
            return None
        else:
            return result

    def add_consensus(self, cutoffs=(50, 60, 70, 80, 90), separator='=', consensus_gap="."):
        """
        Add consensus rows to the MSA object.

        Parameters
        ----------
            cutoffs : list of integers
                Minimum frequency, in each column, required for
                the conserved residue or residue category
            separator : character
                Add a row with filled with the user-defined
                character to separate the alignment and the
                consensus lines

        Returns
        -------
        A copy of the MSA with consensus sequences
        """
        result = self.copy()
        cx = []
        frequencies = result.residue_frequencies
        for cutoff in reversed(sorted(cutoffs)):
            consensus = self.consensus(cutoff, frequencies=frequencies, consensus_gap=consensus_gap)
            cx.append([ f'Consensus {cutoff}%', consensus, len(consensus.replace('.','')), "consensus", f'Consensus {cutoff}%' ])
        if separator:
            cx.append([ 'separator', "".join([ separator for x in range(0,self.get_alignment_length()) ]), self.get_alignment_length(), "view", "separator" ])
        result.df = pd.concat([pd.DataFrame(cx, columns=self._reserved_columns), result.df])
        return result

    def add_hhpred(self, hhpred_data, hhpred_id, hhsearch=False):
        '''
        Merge the MSA and a pairwise alignment generated by HH-suite.

        The resulting MSA will include both the sequence and 
        secondary structure of the selected HH-suite result.

        Parameters
        ----------
        hhpred_data : string or file path
            HH-suite output file or string
        hhpred_id : string
            HH-suite hit identifier
        '''
        import re

        # Parser
        try:
            with open(hhpred_data, 'r') as f:
                texto = f.read()
        except:
            texto = hhpred_data

        # Parse
        match = re.findall(f'No {hhpred_id}(.+?)No ', texto, re.DOTALL)[0].strip()
        target = re.findall('T Consensus.*?\nT\s(.*?)\s',match, re.MULTILINE)[0]
        if hhsearch:
            query = re.findall('Query\s+(.*?)\s',texto, re.MULTILINE)
        else:
            query = re.findall('Q ss_pred.*?\nQ\s(.*?)\s',match, re.MULTILINE)
        
        if query:
            query = query[0]
        T_con = ''.join(re.findall('T Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_target = ''.join(re.findall(f'T\s+{target}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        sequence_query = ''.join(re.findall(f'Q\s+{query}.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        T_ss_pred = ''.join(re.findall('T ss_pred\s+(.*?)$',match,  re.MULTILINE))
        T_ss_dssp = ''.join(re.findall('T ss_dssp\s+(.*?)$',match,  re.MULTILINE))
        Q_ss_pred = ''.join(re.findall('Q ss_pred\s+(.*?)$',match,  re.MULTILINE))
        Q_con = ''.join(re.findall('Q Consensus.*?\d (.*?)\s*\d',match,  re.MULTILINE))
        if hhsearch:
            structure_df = pd.DataFrame({'query':list(sequence_query), 'target':list(sequence_target), 'ss':list(T_ss_dssp)})
        else:
            structure_df = pd.DataFrame({'query':list(sequence_query), 'query_pred': list(Q_ss_pred), 'target':list(sequence_target), 'ss':list(T_ss_dssp)})

        # join aln to hhpred results
        result = self.copy()
        aln = pd.Series(list(result.df.query('id == @query').sequence.iloc[0])).where(lambda x: x!='-').dropna().reset_index().rename({'index':'position', 0:'sequence'}, axis=1)
        s_aln = ''.join(aln.sequence).find(''.join(structure_df['query'].where(lambda x : x !='-').dropna()).upper())
        e_aln = s_aln + len(''.join(structure_df['query'].where(lambda x : x !='-').dropna())) -1
        aln.loc[s_aln : e_aln , 'dssp'] = structure_df[structure_df['query'] != '-'].ss.to_list()
        if not hhsearch:
            aln.loc[s_aln : e_aln , 'sspred'] = structure_df[structure_df['query'] != '-'].query_pred.to_list()
        aln.loc[s_aln : e_aln , target] = structure_df[structure_df['query'] != '-'].target.to_list()
        aln = pd.Series(list(result.df.iloc[0].sequence)).to_frame().join(aln.set_index('position', drop=True)).fillna('-')
        aln = aln.T.apply(lambda x: ''.join(list(x)), axis=1).reset_index().rename({'index': 'id', 0: 'sequence'}, axis=1) 
        result.df =  pd.concat([aln,result.df]).iloc[2:].reset_index(drop=True)

        return result

    def add_pdb(self, pdb_id, chain_id='A', pdb_file=None, pdb_dir=config['pdb_dir']):
        """
        Add structure annotation from PDB.

        Parameters
        ----------
        pdb_id : string
            4-letter PDB code or arbitrary string
        chain_id : string, defaul 'A'
            Target chain identifier.
        pdb_file : string
            PDB file or URL of remote PDB file.
            If pdb_file set to 'esm', sequence selected by pdb_id will be
            forwarded to esmfold api, only applies to sequences with less than 400aa in length.
        pdb_dir : string
            Path to a local PDB mirror or a directory
            where PDB files are stored

        See also
        --------
        Rotifer's configuration
        """
        import tempfile
        from Bio import pairwise2
        from Bio.PDB import PDBParser
        from Bio.PDB.DSSP import DSSP
        import Bio.PDB.PDBList as PDBList
        result = self.copy()

        if not isinstance(pdb_id,list):
            pdb_ids = [ pdb_id ]
        else:
            pdb_ids = pdb_id

        for i in range(0, len(pdb_ids)):
            pdb_id = pdb_ids[i]
            # Find local file or download and then open it!
            if pdb_file:
                if not os.path.exists(pdb_file):
                    if os.path.exists(os.path.join(pdb_dir,pdb_file)):
                        pdb_file = os.path.exists(os.path.join(pdb_dir,pdb_file))
                    else:
                        if pdb_file == 'esm':
                            # Sends first sequence to ESM-Fold API 
                            import urllib
                            data = self.filter(keep=pdb_id).to_string(output_format='fasta-2line', remove_gaps=True).split('\n')[1].encode('utf-8')
                            req = urllib.request.Request(url="https://api.esmatlas.com/foldSequence/v1/pdb/", data=data, method='POST')
                            pdb_data = urllib.request.urlopen(req).read()
                            pdb_file = open(rotifer.config['cache']+"/"+pdb_id[0]+".pdb", "wb")
                            pdb_file.write(pdb_data)
                            pdb_file = open(rotifer.config['cache']+"/"+pdb_id[0]+".pdb", "r")
                            pdb_file.flush()
                            pdb_file.seek(0)
                        else:
                            # Try using pdb_file as URL
                            import urllib
                            pdb_data = urllib.request.urlopen(pdb_file).read().decode()
                            pdb_file = open(rotifer.config['cache']+ "/" + pdb_id + ".pdb", "wt")
                            pdb_file.write(pdb_data)
                            pdb_file = open(rotifer.config['cache']+ "/" + pdb_id + ".pdb", "r")
                            pdb_file.flush()
                            pdb_file.seek(0)

            else:
                # No file!
                if os.path.exists(os.path.join(pdb_dir,pdb_id[1:3].lower(),"pdb"+pdb_id.lower()+".ent.gz")):
                    # Search PDB code in local PDB mirror
                    pdb_file = os.path.join(pdb_dir, pdb_id[1:3].lower(), "pdb"+pdb_id.lower()+".ent.gz")
                else:
                    pdb_file = PDBList(verbose=False).retrieve_pdb_file(pdb_id.lower(), file_format='pdb', pdir=rotifer.config['cache'])

            if isinstance(pdb_file,str) and os.path.exists(pdb_file):
                if pdb_file[-3:] == '.gz':
                    import gzip
                    orig = gzip.open(pdb_file,"rt")
                    pdb_file = tempfile.NamedTemporaryFile(suffix=".pdb", delete=True, mode="r+t")
                    pdb_file.write(orig.read())
                    pdb_file.flush()
                    pdb_file.seek(0)
                    orig.close()
                else:
                    pdb_file = open(pdb_file,"rt")

            # Parse PDB file to a Pandas DataFrame
            dssp_columns  = ['pdb_id','chain','c1','c2','c3','idx','aa','secstr','rASA','phi','psi']
            dssp_columns += ['NH_O_1_relidx','NH_O_1_energy','O_NH_1_relidx','O_NH_1_energy']
            dssp_columns += ['NH_O_2_relidx','NH_O_2_energy','O_NH_2_relidx','O_NH_2_energy']
            dssp = PDBParser().get_structure(pdb_id, pdb_file)
            dssp = DSSP(dssp[0], pdb_file.name, file_type="PDB")
            dssp = pd.DataFrame([ (pdb_id,x[0],*x[1],*dssp[x]) for x in dssp.keys() ], columns=dssp_columns)
            dssp_to_ehc = {"I":"C","S":"C","H":"H","E":"E","G":"C","B":"E","T":"C","C":"C"}
            dssp.secstr = dssp.secstr.map(dssp_to_ehc).fillna('-')
            dssp = dssp.query('chain == @chain_id')
            pdb_file.close()
            if os.path.exists(pdb_file.name) and pdb_file.name[0:len(rotifer.config['cache'])] == rotifer.config['cache']:
                os.remove(pdb_file.name)

            # Search for pdb sequence in the aligment
            pdb_from_pdb = ''.join(dssp.aa.to_list())
            pdb_index = result.df[result.df.id.str.contains(pdb_id, case=False)]
            pdb_index = pdb_index[~pdb_index.id.str.contains('ss_from')] #Avoids annotating itself
            pdb_index = pdb_index.index[0]
            pdb_in_aln = result.df.loc[pdb_index].sequence.replace('-', '')
            ali = pairwise2.align.localxx(pdb_in_aln, pdb_from_pdb)[0]
            ss = pd.Series(list(ali.seqB)).where(lambda x : x != '-').dropna().to_frame()
            ss['structure'] = dssp.secstr.to_list()
            pdb_df = pd.Series(list(ali.seqA)).rename('seq').to_frame().join(ss['structure']).fillna('-').query(' seq != "-"')
            to = pd.Series(list(result.df.loc[pdb_index].sequence)).where(lambda x: x !='-').dropna().rename('ung').to_frame()
            to['structure'] = pdb_df.structure.to_list()
            pdbn = f'ss_from:{pdb_id}_{chain_id}'
            pdbss = ''.join(pd.Series(list(result.df.loc[pdb_index].sequence)).to_frame().join(to).fillna('-').structure.to_list())
            result.df = pd.concat([pd.DataFrame([[pdbn,pdbss,len(pdbss),"residue_annotation",f'Secondary structure from {pdb_id}']], columns=self._reserved_columns),self.df])

        return result

    def add_seq(self, seq_to_add, cpu=12, fast=False, local_database_path=[os.path.join(rotifer.config['data'],'fadb','nr','nr')]):
        import tempfile
        from Bio.SeqRecord import SeqRecord
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.devel.beta import sequence as rdbs

        # Make sure input is a list
        if not isinstance(seq_to_add,list):
            seq_to_add = [ seq_to_add ]

        # Dump input sequences
        with tempfile.TemporaryDirectory() as tmpdirname:
            # Save current sequences to temporary file
            self.to_file(f'{tmpdirname}/seqaln') 

            # Save new sequences to temporary file
            if isinstance(seq_to_add[0],type(self)):
                seq_to_add.to_file(f'{tmpdirname}/acc.fa')
            elif isinstance(seq_to_add[0],SeqRecord):
                SeqIO.write(seq_to_add, f'{tmpdirname}/acc.fa', "fasta")
            else:
                seqobj = rdbs.sequence(seq_to_add, local_database_path=local_database_path)
                seqobj.df = seqobj.df.drop_duplicates() # Avoid duplicates
                seqobj.to_file(f'{tmpdirname}/acc.fa')

            # Run MAFFT
            child = f'mafft --thread {cpu} --add {tmpdirname}/acc.fa'
            if not fast:
                child = f'{child} --maxiterate 1000 --localpair'
            child = Popen(f'{child} {tmpdirname}/seqaln', stdout=PIPE, shell=True).communicate()

        # Load output alignment
        result = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        result.file_path = 'from add_seq function'
        return result

    def consensus(self, cutoff=50, frequencies=None, consensus_gap="."):
        '''
        Generate consensus strings for the alignment object

        Parameters
        ----------
        cutoff : integer or float, default 50
          Minimum frequency for conserved residues or categories
        frequencies: Pandas DataFrame
          Frequencies of residues at each column
          See sequence.residue_frequencies

        Returns
        -------
        A string representing the consensus sequence
        '''

        # Ranking of amino acid categories
        aa_type_dict =  {
            'a':6,
            'l':4,
            'h':8,
            '+':3,
            '-':1,
            'c':7,
            'p':10,
            'o':0,
            'u':5,
            's':9,
            'b':11,
            '_':12,
            consensus_gap: 13,
            'gap':14,
        }

        # Copying frequency table and building consensus
        if isinstance(frequencies,pd.DataFrame) and not frequencies.empty:
            result = frequencies
        else:
            result = self.residue_frequencies
        result.rename({'gap':consensus_gap}, inplace=True)
        result = pd.concat([result, pd.DataFrame([[cutoff+100]*len(result.columns)], columns=result.columns, index=[consensus_gap])])
        result = result.melt(ignore_index=False).reset_index().rename({'index':'aa', 'variable':'position', 'value':'freq'}, axis=1)
        result['ranking'] = result.aa.map(aa_type_dict)
        result = result.query(f'freq >= {cutoff}').sort_values(['position','ranking'], na_position='first').drop_duplicates(subset='position')
        return ''.join(result.aa.to_list())

    @property
    def freq_table(self):
        '''
        DEPRECATED: use residue_frequencies instead!
        Residue and residue categories frequencies.
        '''
        return self.residue_frequencies

    @property
    def residues(self):
        '''
        Access alignment columns as a Pandas dataframe.
        '''
        return self.df.query('type == "sequence"').sequence.str.split("", expand=True).loc[:,1:self.get_alignment_length()]

    @property
    def residue_frequencies(self):
        '''
        Column-wise frequencies of residues and residue categories.
        '''
        freq_df = self.residues
        freq_df = freq_df.apply(pd.Series.value_counts).fillna(0).astype(int)/len(freq_df)*100 
        freq_df.rename({'-':'gap','.':'gap','?':'X'}, inplace=True)

        aromatic = ['F','Y', 'W', 'H']
        alifatic = ['I','V','L']
        hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
        positive = ['H', 'K', 'R']
        negative = [ 'D', 'E']
        charged = positive + negative
        polar = charged + ['Q', 'N', 'S', 'T','C']
        alcohol = ['S','T']
        tiny = ['G', 'A', 'S']
        small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
        big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

        aa_type_names =  {'aromatic':[aromatic,'a'],
                         'alifatic':[alifatic,'l'],
                         'hydrophobic': [hydrophobic,'h'],
                         'positive':[positive,'+'], 
                         'negative':[negative,'-'],
                         'charged':[charged,'c'], 
                         'polar':[polar,'p'],
                         'alcohol':[alcohol,'o'],
                         'tiny':[tiny,'u'],
                         'small':[small,'s'],
                         'big':[big,'b'],
                         'all_aa':[all_aa, '_']
                         }

        missing = set(all_aa) - set(freq_df.index.tolist())
        for x in missing:
            freq_df.loc[x] = 0
        for x in aa_type_names.keys():
            freq_df = pd.concat([freq_df,pd.DataFrame({aa_type_names[x][1]:freq_df.loc[aa_type_names[x][0]].sum()}).T])
        return freq_df

    def to_a3m(self, name=None, remove_gaps=False):
        name = self.name
        if not name:
            name = self.df.id.iloc[0]
        a3m = f'#{name}\n'
        ss = self.df.query('type == "structure prediction"')
        if not ss.empty:
            a3m += ">ss_pred " + ss.query('id != "conf"').id.iloc[0].upper() + "\n"
            a3m += ss.query('id != "conf"').sequence.iloc[0].replace(" ","-") + "\n"
            a3m += ">ss_conf " + ss.query('id != "conf"').id.iloc[0].upper() + "\n"
            a3m += ss.query('id == "conf"').sequence.iloc[0].replace(" ","0") + "\n"
        a3m += self.to_string(remove_gaps=remove_gaps)
        return a3m

    def to_bioalign(self):
        """
        Convert the MSA to BioPython's Bio.Align.MultipleSeqAlignment.

        See also
        --------
        Bio.Align.MultipleSeqAlignment : MSA object in BioPython
        Bio.AlignIO : BioPython I/O modules for alignments
        """
        from Bio.Align import MultipleSeqAlignment
        return MultipleSeqAlignment(self.to_seqrecords())

    def to_file(self, file_path=None, output_format='fasta', annotations=None, remove_gaps=False):
        file_ext_dict = {'pkl':'pickle', 'a3m':'a3m', 'hmm':'hmm'}
        file_ext = file_path.rsplit('.', maxsplit=1)
        if len(file_ext) == 2 and file_ext[1] in file_ext_dict:
            output_format = file_ext_dict[file_ext[1]]

        if output_format == "a3m":
            fh = open(file_path,"wt")
            fh.write(self.to_a3m())
            fh.flush()
            fh.close()
            return len(self)
        elif output_format == "pickle":
            import pickle
            with open(file_path,"wb") as fh:
                pickle.dump(self,fh)
            return f'{file_path} pickle file saved'
        elif output_format =="hmm":
            import subprocess
            import tempfile
            import os
            with tempfile.TemporaryDirectory() as temp_dir:
                fasta_file = os.path.join(temp_dir, "alignment.fa")
                self.to_file(fasta_file)
                hmm_file = os.path.join(temp_dir, 'output.hmm')
                cmd = ["hmmbuild", "-n", self.name, hmm_file, fasta_file]
                subprocess.run(cmd, check=True)
                with open(hmm_file, "r") as f:
                    lines = f.readlines()
                for i, line in enumerate(lines):
                    if line.startswith("NAME"):
                        lines.insert(i + 1, f"DESC  {self.description}\n")
                        break    
                with open(file_path, "w") as f:
                    f.writelines(lines)
            return f'{file_path} hmm file saved'

        else:
            return SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), file_path, output_format)

    def to_hist(self, bins=10):
        """
        This method displays the distribution of sequence
        lengths for all the sequences in the alignment, 
        after removal of all their gaps.

        Parameters
        ----------
        bins : int
            Number of histogram bins to discretize the distribution

        Returns
        -------
        None

        Examples
        --------
        >>> aln.hist(50)
        """
        from ascii_graph import Pyasciigraph
        import collections
        collections.Iterable = collections.abc.Iterable
        print(f'Total proteins: {len(self.df)}')
        a = self.df['length'].value_counts().to_frame().reset_index()
        a.columns = ['length', 'count']
        a = a.sort_values('length')
        a['raw_bin'] = pd.cut(a['length'],bins,precision=0)
        a['bin'] = a.raw_bin.apply(lambda x : '{} - {}'.format(int(x.left),int(x.right)))
        test = a.groupby('bin').agg({'count':'sum'}).reset_index().apply(tuple, axis=1)
        graph = Pyasciigraph()
        for line in  graph.graph('count \t sequence size', test):
            print(line)

    def to_seqrecords(self, annotations=None, remove_gaps=False):
        """
        Convert the MSA to a list of BioPython's SeqRecord objects.
        """
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq

        # Access and filter internal dataframe
        df = self.df
        if annotations:
            df = df.query('type == "sequence" or type in @annotations').reset_index(drop=True)
        else:
            df = df.query('type == "sequence"')

        # Convert each row to a SeqRecord
        result = []
        for row in list(df.T.to_dict().values()):
            if remove_gaps:
                row["sequence"] = row["sequence"].replace("-","")
            if "description" not in row:
                row["description"] = "unknown sequence"
            if not isinstance(row["description"], str):
                row["description"] = "unknown sequence"
            result.append(SeqRecord(id=row["id"], seq=Seq(row["sequence"]), description=row["description"]))
        return result

    def to_string(self, output_format='fasta', annotations=None, remove_gaps=False):
        sio = StringIO("")
        SeqIO.write(self.to_seqrecords(annotations=annotations, remove_gaps=remove_gaps), sio, output_format)
        return sio.getvalue()

    def align(self, method='famsa', cpu=12, region=False, inplace=False):
        """
        Rebuild the alignment using Mafft.

        Parameters
        ----------
        method : string, default is famsa
            famsa ,
            mafft runs maaft with defaul parameter 
            linsi runs mafft  using the Smith-Waterman algorithm with options ```--maxiterate 1000``` and ```--localpair```
            kalign 
        cpu : integer, default is 12
            Number of threads to use, kalign do not have thread option
        region: tuple of integers, default: False
            Alignment columns to realign
        inplace: boolean, default False
            Align sequences inplace, i.e. don't create a new object

        Returns
        -------
        A new MSA object.
        """
        from subprocess import Popen, PIPE, STDOUT

        # Make sure the index is defined by this function
        self.df.reset_index(drop=True) 

        # Prepare temporary file
        seq_string = self.copy()
        seq_string.df = seq_string.df.reset_index().drop(['id'], axis=1).rename({'index':'id'}, axis=1).astype({'id':'str'})
        if region:
            seq_string = seq_string.slice((region[0],region[1]))
        seq_string = seq_string.to_string(remove_gaps=True).encode()

        # Run alignment software and load output
        if method =='mafft':
            child = Popen(f'cat|mafft  --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method == 'linsi':
            child = Popen(f'cat|mafft  --maxiterate 1000 --localpair --thread {cpu} -' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method =='famsa':
            child = Popen(f'cat|famsa -t {cpu} -v /dev/stdin STDOUT' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        elif method =='kalign':
            child = Popen(f'cat|kalign' , stdin=PIPE, stdout=PIPE,shell=True).communicate(input=seq_string)
        aligned = self.from_string(child[0].decode("utf-8"), input_format = 'fasta')
        # Make sure aligned index matches the index from self
        aligned.df = aligned.df.astype({'id':int}).set_index('id')

        if inplace:
            result = self
        else:
            result = self.copy()

        # The following rows depende on the index of r1, r2, 
        if region:
            r1 = self.slice((1, region[0] - 1))
            r2 = self.slice((region[1] + 1, len(self.df.iloc[0,1])))
            result.df.sequence = r1.df.sequence + aligned.df.sequence + r2.df.sequence
        else:
            result.df.sequence = aligned.df.sequence

        if not inplace:
            return result


    def _view_groups(self,
            groupby=None,
            min_group_size=2,
            group_separator=np.nan,
            color=True,
            scale=True,
            pager='less -SR',
            *args, **kwargs):
        df = []

        # Add consensus for the full dataframe
        if 'consensus' in kwargs and kwargs['consensus']:
            df = self._view_sequence(color=False, scale=True, pager=None, *args, **kwargs)
            df = df.query('type != "sequence"').copy()
            df.id = df.id.str.replace("Consensus",f'Full ')
            df[groupby] = 'Full'
            df = [ df ]
            if pager and color:
                header = pd.Series(df[0].columns, index=df[0].columns).to_frame().T.astype(str)
                df = [ header ] + df
            if group_separator:
                empty = [[ group_separator * df[0][x].astype(str).str.len().max() for x in df[0].columns ]]
                empty = pd.DataFrame(empty, columns=df[0].columns)
                df.append(empty)

        # Processing each group
        small = []
        blocks = self.df.groupby(groupby)
        for group in blocks.size().sort_values(ascending=False).reset_index().values.tolist():
            sz = group.pop()
            if len(group) == 1:
                group = group[0]

            # Load sub-alignment
            aln = self.copy()
            aln.df = blocks.get_group(group)

            # Isolate small groups
            if sz < min_group_size:
                small.append(aln.df)
                continue

            # Format and add to stack
            aln = aln._view_sequence(color=False, scale=scale, pager=None, *args, **kwargs)
            aln[groupby] = group

            # Add block to stack
            df.append(aln)
            if group_separator:
                empty = pd.DataFrame([[ group_separator * aln[x].astype(str).str.len().max() for x in aln.columns ]], columns=aln.columns)
                df.append(empty)

        # Merge and process small groups
        if small:
            small = pd.concat(small, ignore_index=True)
            aln = self.copy()
            aln.df = small
            aln = aln._view_sequence(color=False, scale=scale, pager=None, *args, **kwargs)
            aln[groupby] = "small"
            df.append(aln)

        # Merge all sub-alignments
        df = pd.concat(df, ignore_index=True).fillna("")

        # Choose coloring
        if color:
            df.sequence = rpf.to_color(df.sequence, padding='left')
            out = df.to_string(header=False, index=False) + "\n"
        else:
            out = df.to_string(index=False) + "\n"

        if pager:
            from IPython.core.page import page
            page(out, pager_cmd=pager)
        else:
            return df

    def _view_sequence(self,
            color=True,
            scale=True,
            consensus=True,
            consensus_gap='.',
            separator="=",
            interval=10,
            columns=True,
            sample=None,
            pager='less -SR',
            *args, **kwargs):
        df = self.copy()

        # Choose columns to display
        if isinstance(columns,bool):
            if not columns:
                df.df = df.df[self._reserved_columns]
        else:
            if isinstance(columns,list):
                df.df = df.df[self._reserved_columns + columns]
            else :
                logger.error('columns should be either a list or a bool')

        if consensus:
            if isinstance(consensus,bool):
                kwargs = {}
            elif isinstance(consensus,int):
                kwargs = {'cutoffs': (consensus,)}
            else:
                kwargs = {'cutoffs': consensus}
            df = df.add_consensus(separator=separator, consensus_gap=consensus_gap, **kwargs)

        # Shift work to the internal dataframe, erasing the sequence object
        df = df.df.copy().astype(str)

        # Add scale bar rows
        if scale:
            scale = self._scale_bar(self.get_alignment_length(), interval=interval).astype(str)
            df = pd.concat([ scale, df ]) 

        # Select sequences
        if sample != None:
            seqs =  df.query('type == "sequence"')
            if len(seqs) > sample:
                seqs = seqs.sample(sample)
                other = df.query('type != "sequence"')
                df = pd.concat([ other, seqs ])

        # Deal with coloring
        header=True
        if color:
            if pager:
                header = pd.Series(df.columns, index=df.columns).to_frame().T
                df = pd.concat([ header, df ])
                header = False
            df.sequence = rpf.to_color(df.sequence, padding='right')

        if pager:
            from IPython.core.page import page
            df = df.to_string(index=False, header=header) + "\n"
            page(df, pager_cmd=pager)
        else:
            return df

    def view(self,
             groupby=None,
             min_group_size=2,
             group_separator=np.nan,
             color=True,
             scale=True,
             consensus=True,
             consensus_gap=".",
             separator="=",
             interval=10,
             columns=True,
             sample=None,
             pager='less -SR',
             *args, **kwargs):
        from rotifer.view import functions
        jupyter = functions.is_running_in_jupyter()
        """
        See a colored version of the alignment with annotations.

        Parameters
        ----------

        Usage
        -----
        >>> from rotifer.devel.beta import sequence as rdbs
        >>> aln = rdbs.sequence("my.aln","fasta")
        >>> aln.add_cluster(coverage=0, identity=0).view(groupby="c0i0")

        Parameters
        ----------
        groupby: string, default None
            See colored versions of all sub-alignments.
            Choose a column to group sequences and show alignments
            and annotations for each group.
        min_group_size: int, default None
            Size of the smallest groups to be shown as a block.
        group_separator:
            Character to add to the rows separating groups.
            Set it to None to disable.
        groupby: string, default None
            Choose a column for grouping rows. A separate
            alignment view will be created for each group.
        min_group_size: integer, default 2
            Minimum number of sequences in groups that will
            be displayed in their own separate view.
        color: bool
            Color sequence residues
        scale: bool, default True
            If set to False, no scale is shown
        consensus: bool, default True
            Whether to display consensus rows
        consensus_gap: string, default '.'
            Character to represent (mostly) gapped columns in the
            consensus string(s).
        separator: character
            Single character to fill row separating the alignment
            and the consensus sequence
        interval: integer, default 10
            Interval between position marks in the scale
        columns: bool or list of strings, default is True
            List of annotation columns to show
            If set to False, only the default columns are shown
            If set to True, all columns in the internal DataFrame are shown
        sample: int, default None
          Randomly select and display up to this number of 
          sequences per group.

          Note: You can set this number to 0 hide ALL sequences!
        pager: str, default 'less -SR'
          External command to use as viewer.

        Returns
        -------
          Nothing if pager is set.
          When pager is set to None, the formatted dataframe is returned.
        """
        method = "_view_sequence"
        if groupby:
            if sample == 0: # Consensus only
                scale = False
                separator = None
            method = '_view_groups'
        method = self.__getattribute__(method)
        if jupyter:
            functions.display_html_popup_from_file(self.to_html())
            return
        return method(groupby=groupby,
                      min_group_size=min_group_size,
                      group_separator=group_separator,
                      color=color,
                      scale=scale,
                      consensus=consensus,
                      consensus_gap=consensus_gap,
                      separator=separator,
                      interval=interval,
                      columns=columns,
                      sample=sample,
                      pager=pager,
                      *args, **kwargs)

    def hhblits(self, databases=config['databases'], database_path=config['databases_path'], view=True, cpu=18):
        """
        Search the alignment against a HMM databases using hhsearch.

        Parameters
        ----------
        databases : list of strings, default ['pfam','pdb70']
            List of HMM databases to include in the search
        database_path : string, default is ROTIFER_DATA/hhsuite
            Path to the directory where the HMM databases are stored
        view : bool, default True

        Returns
        -------
            A tuple of two elements:
            - The HHblits output as a string
            - HHblits output as a Pandas DataFrame
              See rotifer.io.hhsuite

        See also
        --------
            rotifer environment configuration

        Examples
        --------
        Load alignment in multi-FASTA format and compare it to Pfam

        >>> aln = sequence("myaln.aln")
        >>> (hhout,hhdf) = aln.hhblits(databases=['pfam'])
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.io.hhsuite import read_hhr

        dbs = " ".join([ " -d " + os.path.join(database_path, x) for x in databases ])
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/seqaln')
            child = f'hhblits -i {tmpdirname}/seqaln {dbs} -M 50 -cpu {cpu} -o {tmpdirname}/seqaln.hhr'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            hhtable = read_hhr(f'{tmpdirname}/seqaln.hhr')
            with open(f'{tmpdirname}/seqaln.hhr') as f:
                hhr_result = f.read()
            if view:
                from IPython.core.page import page
                page(hhr_result)
            return (hhr_result, hhtable)

    def hhsearch(self, databases=config['databases'], database_path=config['databases_path'], view=True, cpu=18):
        """
        Search the alignment against a HMM databases using hhsearch.

        Parameters
        ----------
        databases : list of strings, default ['pfam','pdb70']
            List of HMM databases to include in the search
        database_path : string, default is ROTIFER_DATA/hhsuite
            Path to the directory where the HMM databases are stored
        view : bool, default True

        Returns
        -------
            A tuple of two elements:
            - The HHsearch output as a string
            - HHsearch output as a Pandas DataFrame
              See rotifer.io.hhsuite

        See also
        --------
            rotifer environment configuration

        Examples
        --------
        Load alignment in multi-FASTA format and compare it to Pfam

        >>> aln = sequence("myaln.aln")
        >>> (hhout,hhdf) = aln.hhsearch(databases=['pfam'])
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer.io.hhsuite import read_hhr

        dbs = " ".join([ " -d " + os.path.join(database_path, x) for x in databases ])
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/seqaln')
            child = f'hhsearch -i {tmpdirname}/seqaln {dbs} -M 50 -cpu {cpu} -o {tmpdirname}/seqaln.hhr'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            hhtable = read_hhr(f'{tmpdirname}/seqaln.hhr')
            with open(f'{tmpdirname}/seqaln.hhr') as f:
                hhsearch_result = f.read()
            if view:
                from IPython.core.page import page 
                page(hhsearch_result)
            return (hhsearch_result, hhtable)

    def community(self):
        import tempfile
        import community
        import networkx as nx
        import numpy as np

        x = self.copy()
        with tempfile.TemporaryDirectory() as tmpdirname:
            self.to_file(f'{tmpdirname}/tt', remove_gaps=True)
            os.system(f'mmseqs easy-search {tmpdirname}/tt {tmpdirname}/tt {tmpdirname}/tt.m8 {tmpdirname}/tmp')
            df = pd.read_csv(f'{tmpdirname}/tt.m8',sep="\t", names='source target pident length mismatch gapopen qstart qend sstart send evalue bitscore'.split())
            dfc = df
            dfc.source, dfc.target = np.where(dfc.source > dfc.target , [dfc.source, dfc.target], [dfc.target, dfc.source])
            dfc = dfc.sort_values('evalue').drop_duplicates(['source', 'target'])
            dfc['evalue_t'] = - np.log10(dfc.evalue)
            G = nx.from_pandas_edgelist(dfc[['source', 'target', 'evalue_t']].query('evalue_t >= 3'), edge_attr='evalue_t')
            partition = community.best_partition(G,weight='weight')
            c = pd.DataFrame.from_dict(partition,orient='index').reset_index().rename(
                {'index': 'id', 0: 'community'}, axis=1)
            x.df = x.df.merge(c)
        return x 

    def trim(self, max_perc_gaps=80, minimum_length=1):
        '''
        Remove alignment columns based on column statistics.

        Parameters
        ----------
        max_perc_gaps: integer or float
          Maximum relative frequency of gaps

        minimum_length : integer, default 1
          Minimum number of consecutive bad quality columns.
          Columns located in regions shorter than this threshold
          will not removed.

        Examples
        --------
        Remove all columns with more than 70% of gaps.
        >>> aln.trim(70)
        '''
        columns_to_keep = self.residue_frequencies.T
        result = self.copy()
        if 'gap' not in columns_to_keep.columns:
            return result
        else:
            columns_to_keep = columns_to_keep.query('gap <= @max_perc_gaps').T.columns.to_list()
            result.df['sequence'] = result.residues.loc[:, columns_to_keep].sum(axis=1)
            result._reset()
            return result

    def add_jpred(self, email=False, symbol=False):
        ''' 
        Function to add secondary structure from the Jpred server
        ∞ = Alpha helix
        > = Beta strand
        - = Turn
        '''
        import tarfile
        import tempfile
        import re
        from subprocess import Popen, PIPE, STDOUT, check_output
        import os
        import io

        if not email:
            from rotifer.db.ncbi import NcbiConfig
            email = NcbiConfig['email']

        result = self.copy()
        with tempfile.TemporaryDirectory() as tmpdirname:
            cd = os.getcwd()
            os.chdir(tmpdirname)
            self.to_file('seqaln')
            child = f'jpredapi submit mode=msa format=fasta email={email} file=seqaln name=rotifer'
            child = check_output(child, shell=True).decode()
            match = re.findall(f'jobid=(.+?)\s', child, re.DOTALL)[0]
            child = f'jpredapi status jobid={match} getResults=yes checkEvery=10'
            child = Popen(child, stdout=PIPE,shell=True).communicate()
            with tarfile.open(f'{match}/{match}.tar.gz', 'r:gz') as tar:
                ss_file = tar.getmember(f'{match}.jnet')
                query_file = tar.getmember(f'{match}.msf.query')
                jnet = pd.read_csv(io.BytesIO(tar.extractfile(ss_file).read()), sep=":", names = ['a', 'b'])
                query = tar.extractfile(query_file).read()
                query = re.findall('>(.+?)\n', query.decode(), re.DOTALL)[0]

            os.chdir(cd)
            jnet = jnet.iloc[0:2,:]
            jnet.b = jnet.b.str.replace(',', '')
            jnet.iloc[0,1] = jnet.query('a  == "jnetpred"').iloc[0,1].replace(',', '')
            if symbol:
                jnet.iloc[0,1] = jnet.iloc[0,1].replace('E', '>').replace('H', "=")

            a1 = pd.Series(list(result.df.query('id ==@query').iloc[0,1])).where(lambda x: x !='-').dropna().rename('seq').reset_index()
            a2 = pd.Series(list(jnet.loc[0, 'b'])).rename('jnet')
            a3 = pd.Series(list(jnet.loc[1, 'b'])).rename('conf')
            a4 = pd.Series(list(result.df.query('id == @query').iloc[0,1]))
            a5 = a1.join(a2).join(a3).set_index('index')
            a6 = pd.concat([a4,a5], axis = 1).fillna(' ').sum().reset_index()
            a6.columns = ['id', 'sequence']
            a6['type'] = 'structure prediction'
            a6.iloc[2:3,2] = 'Jpred'
            a6.iloc[3:4,2] = 'confidance'

            result.df = pd.concat([a6.iloc[2:4,:],result.df])

            return result
    def compact_residue_df(self,
                           consensus,
                           annotations=False,
                           remove_gaps=False,
                           adjust_coordinates=False):
        """
        Function to return a residue df containing a selected consensus as the last line
        consensus: consensus threshold to bild the consensus line
        annotation: can be a list of string of sequence obejct id to be used as annotation on top of aligment
        and will not be colored
        remove gaps: Should add the id of the protein used to be the model to remove gaps on aligment
        adjust coordinate: will fetch sequences to add start and end coordinates 
        """

        aln = self.copy()
        if remove_gaps:
            gaps,gdf, spos,epos = aln.filter('type =="sequence"').gaps_to_numbers(remove_gaps, adjust_coordinates=adjust_coordinates)
            gdf.index = aln.df.query('type == "sequence"').id

        aln_r = aln.residues
        con = pd.Series(list(aln.consensus(consensus)))
        aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
        con.index +=1
        aln_r = pd.concat([aln_r, con.rename(f'consensus/{consensus}%').to_frame().T], axis=0)
        if annotations:
            if isinstance(annotations, str):
                ann = pd.Series(list(aln.df.query('id ==@annotations').sequence.iloc[0]))
                ann.index +=1
                aln_r = pd.concat([ann.rename(annotations).to_frame().T,aln_r])
            else:
                for x in annotations:
                    ann = pd.Series(list(aln.df.query('id ==@x').sequence.iloc[0]))
                    ann.index +=1
                    aln_r = pd.concat([ann.rename(x).to_frame().T,aln_r])

        if remove_gaps:
            if adjust_coordinates:
                if spos not in aln_r.columns:
                    aln_r.insert(0,spos, ['-'] * len(aln_r))
                    aln_r.columns = [x +1 for x in aln_r.columns] # Checar isso!!!!!
                    gaps = [x +1 for x in gaps ]
                    epos = epos + 1
                    gdf.columns = [x if x == 1 else x + 1 for x in gdf.columns]
                if epos not in aln_r.columns:
                    aln_r[epos] = ['-'] * len(aln_r)
            #return aln_r, gdf, gaps       
            aln_r =  aln_r.drop(gaps,axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0','  ')
        return aln_r




    def _to_df_style(self,
                     consensus,
                     annotations=False,
                     remove_gaps=False,
                     adjust_coordinates = False,
                     font_size=6):
        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment
        :output_file: output file name
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :returns: TODO

        """

        import sys
        import pandas as pd
        from rotifer.devel.beta.sequence import sequence
        from rotifer.core.functions import loadConfig
        from rotifer.core  import config as CoreConfig

        #### Loading the color dictionary
        cd = loadConfig(
                ':colors.html_aa_colors',
                system_path=CoreConfig['baseDataDirectory'])

        import numpy as np
        aln = self.copy()

        aln_r = aln.compact_residue_df(consensus,
                     annotations=annotations,
                     remove_gaps=remove_gaps,
                     adjust_coordinates = adjust_coordinates)


        # Funtions to color the algiment:
        def highlight_aln(s):
            import numpy as np
            ### getting the consensus value to map the colors filling na with "  " to color white 
            d = cd[s.fillna('_').iloc[-1]]
            #d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            return np.where(
                s == '  ',
                'color:white;background-color:white',
                np.where(
                    s == s.iloc[-1],
                    f'color:{d["fcolor"]};background-color:{d["color"]}',
                    np.where(
                        s.isin(d['residues']),
                        f'color:{d["fcolor"]};background-color:{d["color"]}',
                        f'color:black;background-color:white')))

        def highlight_consensus(s):
            import numpy as np
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

        #Making slice index where the functions should be applied:
        # One function should be applien only in the consensus row
        # Other function should be appplied only in seq rows

        #### Geting the Consensus line
        idx1 = pd.IndexSlice
        corte = idx1[idx1[f'consensus/{consensus}%'],idx1[:]]
        ###Getting the sequences from the aligment
        #Getting the firs sequence (fs) row to map the slice:
        idx2 = pd.IndexSlice
        fs = aln.df.query('type == "sequence"').id.iloc[0]    
        corte2 = idx2[idx2[fs:],idx2[:]]

        headers = {
            'selector': 'th:not(.index_name)',
            'props': f'''font-size: {font_size}px;
            text-align: left;
            font-family:"Lucida Console", Monaco, monospace;
            color:black;
            background-color:white'''
        }

        if sys.version_info.minor > 8:
            df_style = aln_r.style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide(axis='columns').apply(
                highlight_consensus, subset=corte
            ).set_table_styles(
                [headers]
            )
        else:
            df_style = aln_r.style.set_properties(**{
                'font-size': f'{font_size}px',
                'font-family':'"Lucida Console", Monaco,monospace',
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide_columns().apply(
                highlight_consensus, subset=corte
            ).set_table_styles(
                [headers]
            )
        #if whant to send to latex, replace set_stick... to:to_latex(environment='longtable', convert_css=True)
        return df_style

    def _to_df_style_simple(self, consensus=[50,60,70,80,90,100], annotations=False, remove_gaps=False,consensus_position='top', adjust_coordinates = False, background='black', colors=config['html_colors']):
        from rotifer.core import functions as rcf
        # Defaults
        import yaml
        colors = yaml.load(open(colors), Loader=yaml.Loader)


        def col_fun(df_c):
            df_c = df_c.copy()
            df_c = df_c.map(colors).fillna('').radd("color:")
            return df_c

        if background == 'black':
            bck_color = 'black'
            fg_color = 'white'
        else:
            bck_color = 'white'
            fg_color = 'black'

        aln = self.copy()
        aln_r = aln.residues
        if remove_gaps:
            gaps,gdf,_,_ = aln.gaps_to_numbers(remove_gaps, adjust_coordinates=adjust_coordinates)
            gdf.index = aln.df.query('type == "sequence"').id

        aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
        for x in consensus:
            con = pd.Series(list(aln.consensus(x)))
            con.index +=1
            if consensus_position =='top':
                aln_r = pd.concat([con.rename(f'consensus_{x}%').to_frame().T, aln_r], axis=0)
            else:
                aln_r = pd.concat([aln_r, con.rename(f'consensus_{x}%').to_frame().T], axis=0)
        if annotations:
            if isinstance(annotations, str):
                ann = pd.Series(list(aln.df.query('id ==@annotations').sequence.iloc[0]))
                ann.index +=1
                aln_r = pd.concat([ann.rename(annotations).to_frame().T,aln_r])
            else:
                for x in annotations:
                    ann = pd.Series(list(aln.df.query('id ==@x').sequence.iloc[0]))
                    ann.index +=1
                    aln_r = pd.concat([ann.rename(x).to_frame().T,aln_r])

        if remove_gaps:
            aln_r =  aln_r.drop(gaps,axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0','  ')
                    

        df_style = aln_r.style.set_properties(**{
            'font-size': '12px',
            'font-family':'"Lucida Console", Monaco,monospace',
            "text-align": "center",
            "background-color": bck_color,
            "color": fg_color}
        ).apply(col_fun).hide(axis='columns')
        return df_style    

    def to_html(self,
                output_file=False,
                consensus=[50, 60,70,80,90,100],
                simple=True,
                background='black',
                consensus_position='top',
                annotations=False,
                remove_gaps=False,
                fixed_index=True,
                adjust_coordinates=False):

        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment, if simple=False the consensus cannot be a list
        :output_file: output file name, if no outputfile it will retunr the HTML and can be save as object
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        :Simple: Create a simple html visualization withou fancy highlights etc, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :fixed_index, It makes the index fixed on the html page 
        :returns: TODO

        """
        if simple:
            if background == 'black':
                bck_color = 'black'
                fg_color = 'white'
            else:
                bck_color = 'white'
                fg_color = 'black'
            headers = {
                'selector': 'th:not(.index_name)',
                'props': f'''font-size: 8px;
                text-align: left;
                color:{fg_color};
                background-color:{bck_color}'''
            }
            html = self._to_df_style_simple(
                consensus=consensus,
                annotations=annotations,
                adjust_coordinates = adjust_coordinates,
                remove_gaps=remove_gaps,
                background=background,
                consensus_position=consensus_position
            ).set_table_styles(
            [headers]
        )
            html_template = f'''
            <!DOCTYPE html>
            <html>
            <head>
                <style>
                    body {{
                        background-color: {bck_color};
                        color: {fg_color};
                    }}
                </style>
            </head>
            <body>
                {html.to_html()}
            </body>
            </html>
            '''

        else:
            html = self._to_df_style(
                consensus=consensus,
                annotations=annotations,
                adjust_coordinates = adjust_coordinates,
                remove_gaps=remove_gaps,
            )

        if fixed_index:
            html = html.set_sticky(axis="index")
        if simple:
            if output_file:
                with open(f'{output_file}', 'w') as f:
                    f.write(html_template)
                    return(f'{output_file} saved on the working path')
            else:
                return html_template
        else:        
            with open(output_file, 'w') as f:
                f.write(html.render(table_attributes='cellspacing=0, cellpadding=0'))
                return(f'{output_file} saved on the working path')


    def gaps_to_numbers(self, smodel, adjust_coordinates=False):
        """: This function will retunr a tuple. The first element of the tuple is 
        a list containing the index of columns that contains gap on the model sequence.
        It should be used to drop the columns from the residues df.
        >>> aln.residues.drop(firs element of the tuple output, axis =1).
        
        The seccond element is a df that should be joined with the residues df after the drop:
        >>> aln.residues.join(second elemnet of tuple).sort_index(axis=1)
        :smodel: Sequence used as model to remove the gaps.

        """
        midx = self.df.query('id == @smodel').index
        gaps = self.residues.loc[midx].iloc[0].where(lambda x: x=='-').dropna()
        grouper_map = gaps.reset_index().drop(midx, axis=1).rename(
            {'index':'region'},
            axis=1
        ).eval(
            'inter = ~region.diff().fillna(1).le(1)'
        ).eval(
            'grouper = inter.cumsum()'
        )
        gap_df = grouper_map.groupby('grouper').agg(
            first = ('region','min'),
            last = ('region', 'max'),
            rsize = ('region', 'nunique'),
            rlist = ('region', list)
        )
        num_gaps = grouper_map.set_index('region')[['grouper']].join(
            (self.residues[gaps.index] != '-').T
        ).groupby('grouper').sum().T
        num_gaps.columns = gap_df['first'].to_list()
        #### Creating a list of where the gaps are
        gaps = gaps.index.tolist()
        nongaps = self.residues.loc[midx].iloc[0].where(lambda x: x !='-').dropna().index
        start_position = nongaps[0] 
        end_position =  nongaps[-1]

        if adjust_coordinates:
            def add_cordinates_to_aln(seqobj):
                from rotifer.devel.beta.sequence import sequence
                c = seqobj.copy()
                cx = sequence(c.df.id.tolist())
                cx. df = cx.df.rename({'sequence':'full_sequence', 'length':'full_length'}, axis=1)
                c.df = c.df.merge(cx.df[['id','full_sequence','full_length']], how='left')
                c.df.full_sequence = c.df.full_sequence.fillna(c.df.sequence)
                c.df.full_length = c.df.full_length.fillna(c.df['length'])
                c.df['start'] = c.df.apply(lambda x : 1 + x.full_sequence.find(x.sequence.replace('-','')), axis=1)
                c.df['end'] = c.df.start + c.df['length']
                c.df['C_term'] = c.df['full_length'] - c.df['end']
                c.df =  c.df.drop(['full_sequence'], axis=1)
                #c.df = c.df.fillna(0)
                c.df = c.df.query('type =="sequence"')
                return c
            cord = add_cordinates_to_aln(self).df[['start', 'end', 'C_term']]
            #num_gaps[nongaps[-1] + 1] = cord.C_term.mask(cord.C_term < 0,0 )
            if end_position not in num_gaps.columns:  # Adding column to have the trailling number with the amminoacids counts.
                num_gaps[end_position +1] = 0 
            if 1 not in num_gaps.columns:
                num_gaps.insert(0,1,0) # Adding a column to be used to add the header number of sequences (if gappy there, the other function would already add it)
            
            num_gaps.iloc[:,-1] = cord.end - num_gaps.iloc[:,-1]  ## Checking if fix
            num_gaps.iloc[:,0] = num_gaps.iloc[:,0] +  cord.start ## Checking if fix
            #num_gaps.insert(0,nongaps[0] -1 , cord.start) removed to fix issues
            start_position -= 1
            end_position += 1
            gaps.insert(0, start_position)
            gaps.append(end_position)


        return (gaps,num_gaps, start_position,end_position)

    def _to_df_styleTEX(self, consensus, annotations=False, remove_gaps=False, adjust_coordinates=False):
        """TODO: Docstring for function.

        :consensus: The consensus threshold that should be used to color the aligment
        :output_file: output file name
        :annotation: List of annotations rows that should be keept in the  html file
        The annotation label should be the same as in the id seq object df columm
        :remove_gaps: Query sequence to use as model to remove the gaps, 
        it will add numbers of aminoacid suppressed in the sequence that contain the insertions.
        :returns: TODO

        """
        import sys
        import pandas as pd
        from rotifer.devel.beta.sequence import sequence
        import numpy as np
        aln = self.copy()
        if remove_gaps:
            gaps,gdf, spo,epo = aln.gaps_to_numbers(remove_gaps, adjust_coordinates=adjust_coordinates)
            gdf.index = aln.df.query('type == "sequence"').id

        aromatic = ['F','Y', 'W', 'H']
        alifatic = ['I','V','L']
        hydrophobic = alifatic + [ 'A', 'C', 'F', 'M', 'W', 'Y']
        positive = ['H', 'K', 'R']
        negative = [ 'D', 'E']
        charged = positive + negative
        polar = charged + ['p','Q', 'N', 'S', 'T','C']
        alcohol = ['S','T']
        tiny = ['G', 'A', 'S']
        small = tiny + [ 'V', 'T', 'D', 'N', 'P', 'C']
        big = ['K', 'F', 'I', 'L','M', 'Q', 'R', 'W', 'Y', 'E']
        all_aa = ['G','A','V','I','L','M','F','Y','W','H','C','P','K','R','D','E','Q','N','S','T']

        aa_groups_colors = {'a':[aromatic,  '#2C68F3'],
                            'l':[alifatic, '#2CF3EA'],
                            'h':[hydrophobic,  '#F3E42C90'],
                            '+':[positive,  '#2C68F3'],
                            '-':[negative,  '#F50EF195'],
                            'c':[charged,  '#38F50E'],
                            'p':[polar,  '#0EF5A150'],
                            'o':[alcohol,  '#AE5BF8'],
                            'u':[tiny,  '#EE9C0C'],
                            's':[small,  '#DA147750'],
                            'b':[big,  '#A28694'],
                            '.':[all_aa,  'white'],
                            'G':[all_aa,'white'],
                            'A':[all_aa,'white'],
                            'V':[all_aa,'white'],
                            'I':[all_aa,'white'],
                            'L':[all_aa,'white'],
                            'M':[all_aa,'white'],
                            'F':[all_aa,'white'],
                            'Y':[all_aa,'white'],
                            'W':[all_aa,'white'],
                            'H':[all_aa,'white'],
                            'C':[all_aa,'white'],
                            'P':[all_aa,'white'],
                            'K':[all_aa,'white'],
                            'R':[all_aa,'white'],
                            'D':[all_aa,'white'],
                            'E':[all_aa,'white'],
                            'Q':[all_aa,'white'],
                            'N':[all_aa,'white'],
                            'S':[all_aa,'white'],
                            'T':[all_aa,'white'],
                            ' ':[all_aa,'white'],
                            '  ':[all_aa,'white'],
                            '_':[all_aa,'white']}

        # Geting the residues tabele:
        aln_r = aln.residues
        con = pd.Series(list(aln.consensus(consensus)))
        aln_r = aln_r.set_index(aln.df.query('type == "sequence"').id)
        con.index +=1
        aln_r = pd.concat([aln_r, con.rename('consensus').to_frame().T], axis=0)
        if annotations:
            if isinstance(annotations, str):
                ann = pd.Series(list(aln.df.query('id ==@annotations').sequence.iloc[0]))
                ann.index +=1
                aln_r = pd.concat([ann.rename(annotations).to_frame().T,aln_r])
            else:
                for x in annotations:
                    ann = pd.Series(list(aln.df.query('id ==@x').sequence.iloc[0]))
                    ann.index +=1
                    aln_r = pd.concat([ann.rename(x).to_frame().T,aln_r])

        if remove_gaps:
            aln_r =  aln_r.drop(gaps,axis=1).join(gdf).sort_index(axis=1).fillna(0).astype(int, errors='ignore').astype(str).replace('0','  ')

        # Funtion that works!!!
        def highlight_aln(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            return np.where(
                s == '  ',
                'color:white;background-color:white',
                np.where(
                    s == s.iloc[-1],
                    'color:white;background-color:black',
                    np.where(
                        s.isin(d[0]),
                        f'color:black;background-color:{d[1]}',
                        'color:black;background-color:white')))

        def highlight_consensus(s):
            import numpy as np
            d = aa_groups_colors[s.fillna('  ').iloc[-1]]
            """TODO: Docstring for highlight_consensus.

            :arg1: TODO
            :returns: TODO

            """
            return np.where(
                s.isin(all_aa),
                'color:white;background-color:black',
                f'color:black;background-color:{d[1]}',
                )

        #Making slice index where the functions should be applied:
        # One function should be applien only in the consensus row
        # Other function should be appplied only in seq rows
        idx1 = pd.IndexSlice
        corte = idx1[idx1['consensus'],idx1[:]]
        #Getting the firs sequence (fs) row to map the slice:
        idx2 = pd.IndexSlice
        fs = aln.df.query('type == "sequence"').id.iloc[0]    
        corte2 = idx2[idx2[fs:],idx2[:]]

        if sys.version_info.minor > 8:
            df_style = aln_r.style.set_properties(**{
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide(axis='columns').apply(
                highlight_consensus, subset=corte
            )
        else:
            df_style = aln_r.style.set_properties(**{
                "text-align": "center"}
            ).apply(highlight_aln, axis=0, subset=corte2).hide_columns().apply(
                highlight_consensus, subset=corte
            )
        #if whant to send to latex, replace set_stick... to:to_latex(environment='longtable', convert_css=True)
        return df_style

    def edit(self, consensus=True, scale=True):
        """
        Edit the alignment using VIM.

        Parameters
        ----------
        consensus: bool, deafult True
            Add the default consensus rows
        scale: bool, default True
            Add a scale bar

        Returns
        -------
        rotifer.devel.beta.sequece

        Examples
        --------
        Load alignment in multi-FASTA format and edit

        >>> aln = sequence("myaln.aln")
        >>> aln2 = aln.edit()
        """
        import tempfile
        from subprocess import Popen, PIPE, STDOUT
        from rotifer import config as gc
        import os
        os.environ["VIMMOVE"] = f'{gc["base"]}/etc/vim/vim-move'
        aln = self.copy()
        result = self.copy()
        if consensus:
            aln = aln.add_consensus(separator='=')

        alndf = aln.df.fillna('X')
        cols = list(alndf.dtypes.where(lambda x: x=='object').dropna().index) 
        if scale:
            scale = aln._scale_bar(self.get_alignment_length(), interval=10)
            alndf = pd.concat([ scale, alndf ]).fillna('X') 

        for x in cols:
            alndf[x] = alndf[x].str.pad(alndf[x].str.len().max(), side="right")
        with tempfile.TemporaryDirectory() as tmpdirname:
            alndf.to_csv(f'{tmpdirname}/seqdf.fa',index=False, sep="\t")
            print ("Open NVIM to edit your file.")
            if os.system(f'nvim -u {gc["base"]}/etc/vim/init.vim {tmpdirname}/seqdf.fa') != 0:
                        raise TryNext()
            tmpdf = pd.read_csv(f'{tmpdirname}/seqdf.fa', sep="\t")
            tmpdf = tmpdf[tmpdf['type'].str.contains("sequence")]
            tmpdf = tmpdf.dropna(subset='sequence')
            result = result.filter(keep=tmpdf.id.to_list())
            result.df.sequence = result.df.id.map(tmpdf.set_index('id').sequence.dropna().to_dict())
        
        return result 

    def fetch_neighbors(self, before=5, after=5):
        '''
        Fetch genome neighborhood from the sequences on the sequence object.

        Parameters
        ----------
        after: integer, defult 5.
          Number of neighbors downstream you query, defaul 5.

        before: integer, default 5.
          Number of neighbors upstream you query, defaul 5.

        Examples
        --------
        Collect 3 genes upstream and downstream your query sequence:
        >>> aln.fetch_neighbors(after=3, before=3)
        '''
        gnc = self.cursors['neighborhood']
        gnc.after = after
        gnc.before = before
        pids = self.df.id
        self.ndf = gnc.fetchall(pids)

        return "Neighborhood DF was created at seqobj.ndf" 

    ## Class methods

    @classmethod
    def _scale_bar(self, length, interval=10, start=1):
        # Numbers for the scale bar
        if length > interval:
            position = list(range(int(((start / interval)+1))*interval,length,interval))
            if not length % interval:
                position.append(length)
        else:
            position = [length]

        # Tick marks (vertical bars)
        scale_bar = [ f'{"|":{interval}}' for x in position ]

        # Format numbers
        scale_number = [ f'{str(x):{interval}}' for x in position ]
        if len(str(start)) < position[0] - start - 1:
            scale_number.insert(0, str(start) + " " * (position[0] - start - len(str(start))))
            scale_bar.insert(0, "|" + " " * (position[0] - start - 1))
        else:
            scale_number.insert(0, " " * (position[0] - start))
            scale_bar.insert(0, " " * (position[0] - start))

        # Prepare Series
        scale_bar = "".join(scale_bar)
        scale_number = "".join(scale_number)
        scale_number = scale_number.rstrip() + " " * (length - start - len(scale_number.rstrip()) + 1)
        scale_bar = scale_bar.rstrip() + " " * (len(scale_number) - len(scale_bar.rstrip()))
        scale_dot = "." * len(scale_number)
        scale = pd.Series([scale_number,scale_bar,scale_dot], index=['position', 'bar', 'dot'])
        scale = {'id':scale.index.tolist(), 'sequence':scale, 'length':scale.str.len(), 'type':'view'}
        scale = pd.DataFrame(scale).astype(str)
        return scale

    @classmethod
    def from_seqrecords(cls, input_data):
        """
        Build a MSA object from a list of BioPython objects.

        Parameters
        ----------
        input_data  : list of Bio.SeqRecord.SeqRecord
        """
        return cls(input_data, input_format=type(input_data[0]))

    @classmethod
    def from_string(cls, input_data, input_format='fasta'):
        '''
        This function uses BioPython to parse MSAs from strings.

        Parameters
        ----------
        input_file   : file path or open file handle
        input_format : Biopython supported file format.
                       See Bio.SeqIO and/or Bio.AlignIO.
        '''
        return cls.from_file(StringIO(input_data), input_format=input_format) 

    @classmethod
    def from_file(cls, input_file, input_format='fasta'):
        '''
        Parse multiple sequence alignment files using BioPython.

        Parameters
        ----------
        input_file   : file path or open file handle
        input_format : Biopython supported file format.
                       See Bio.SeqIO and/or Bio.AlignIO.
        '''
        return cls(input_file, input_format=input_format)

# Class method
def concat(alignments, axis=0, by='id'):
    """
    Concatenate data from one or more alignments.

    Parameters
    ----------
    axis: {0/'rows', 1/'columns'}, default 0
      Merge alignments by appending rows (axis == 1), or by
      concatenating sequences sharing the same id (axis = 0).
    """
    if not isinstance(alignments,list):
        logger.error("Expected a list of rotifer.sequence objects as first argument!")
        return
    aln = sequence()
    aln.name = " | ".join([ x.name for x in alignments[1:] ])
    if axis == 0 or axis == "rows":
        aln.df = pd.concat([ x.df.eval(f'source_object = "{x.name}"') for x in alignments ], ignore_index=True)
    else:
        aln.df = alignments[0].df.copy()
        aln.df.description = alignments[0].name
        for y in alignments[1:]:
            y = y.copy()
            y.df.description = y.name
            gap_x = "-" * aln.get_alignment_length()
            gap_y = "-" * y.get_alignment_length()
            aln.df = aln.df.merge(y.df, on=by, how="outer", suffixes=["",f"_{y.name}"])
            aln.df.sequence = aln.df.sequence.fillna(gap_x) + aln.df[f"sequence_{y.name}"].fillna(gap_y)
            aln.df.description = aln.df.description.fillna("") + "," + aln.df[f'description_{y.name}'].fillna("")
            aln.df.drop([ x + f"_{y.name}" for x in y._reserved_columns if x != "id" ], axis=1, inplace=True)
    aln._reset()
    return aln

__doc__ = """
========================
Rotifer Sequence Objects
========================

This module implements classes and methods for dealing with
multiple sequence alignments, also known as MSAs.

Examples
--------
Simplest usage: visualize a (multi)FASTA file.

>>> from rotifer.devel.beta.sequence import sequence
>>> aln = sequence("alignment.aln")
>>> aln.view()

Load alignment from another file format, such as
StockHolm and show the consensus above the alignment:

>>> from rotifer.devel.beta.sequence import sequence
>>> sto = sequence("alignment.sto", format="stockholm")
>>> sto.add_consensus().view()

Parsing a FASTA alignment from an explicit string:

>>> from rotifer.devel.beta.sequence import sequence
>>> aln = sequence('>Seq1\nACFH--GHT\n>Seq2\nACFW--GHS\n')
>>> aln.add_consensus().view()
"""
