__doc__ = """
Helpers for multiple sequence alignments.

``aln_order_by_tree`` replaces the shell function of the same name that
the ROTIFER Lmod module used to define; it is also a command.
"""

def tree_leaves(tree, format='newick'):
    """
    Leaf names of a tree file, from left to right.
    """
    from Bio import Phylo
    return [leaf.name for leaf in Phylo.read(tree, format).get_terminals()]

def aln_order_by_tree(alignment, tree, alignment_format='fasta', tree_format='newick'):
    """
    Reorder the sequences of an alignment to follow the leaves of a tree.

    Sequences missing from the tree are dropped; leaves missing from the
    alignment are ignored.

    Parameters
    ----------
    alignment : str
        Alignment file name.
    tree : str
        Tree file name.
    alignment_format, tree_format : str
        Biopython formats of the alignment and of the tree.

    Returns
    -------
    list of (str, str)
        Sequence identifier and aligned sequence, in tree order.

    Examples
    --------
    >>> from rotifer.tools.alignments import aln_order_by_tree
    >>> ordered = aln_order_by_tree("family.aln", "family.tree")
    """
    from Bio import SeqIO
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(alignment, alignment_format)}
    return [(leaf, sequences[leaf]) for leaf in tree_leaves(tree, tree_format) if leaf in sequences]
