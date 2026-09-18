# POD documentation - main docs before the code

=head1 NAME

Rotifer::DBIC::AnnotationDB::Parser::genbank - 
  load data from Genbank flat files

=head1 SYNOPSIS

  # Using this parser with rotiferDB command line app
  #
  # Note: this will automatically create a graph
  #       and relate all loaded sequences to it

  rotiferDB -if genbank te.gbk

  # Creating a new parser
  use Rotifer::DBIC::AnnotationDB::Parser;
  my $parser = Rotifer::DBIC::AnnotationDB::Parser->create("genbank");
  $parser->load(@ARGV);

=head1 DESCRIPTION

This module loads annotation, database references
and sequence data from files in any of the richly annotated formats
support by BioPerl (GenBank, EMBL, SwissProt, etc).

=head1 AUTHOR

Robson Francisco de Souza

Email robfsouza at gmail.com

=head1 CONTRIBUTORS

Additional contributors names and emails here

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package Rotifer::DBIC::AnnotationDB::Parser::genbank;

use strict;
use warnings;
use autodie qw(:all);
use Bio::SeqIO;
use Carp::Clan qr/^Rotifer::DBIC::AnnotationDB/;
use Moose;
with 'Rotifer::DBIC::AnnotationDB::Role::ParserRole';

=head2 ATTRIBUTES

=head2 seqio_options

 Usage   : $parser->seqio_options
 Function: pass options to Bio::SeqIO->new()
 Value   : hash reference

=cut

has 'seqio_options'   => (
    is       => "rw",
    isa      => "HashRef",
    required => 1,
    builder  => '_default_seqio_options',
    trigger  => \&_seqio_options_trigger,
    );

sub _default_seqio_options {
    return { '-format' => 'genbank' };
}

sub _seqio_options_trigger {
    my ($self, $new, $old) = @_;
    my $nkeys = scalar(keys %$new);
    $new = { %$old, %$new };
    $self->seqio_options($new) if ($nkeys != scalar(keys %$new));
}

=head2 genbank_tags_ontology

 Usage   : $parser->genbank_tags_ontology
 Function: ontology to store genbank annotation tags
 Value   : Rotifer::DBIC::AnnotationDB::Result::Ontology

=cut

has 'genbank_tags_ontology'   => (
    is       => "ro",
    isa      => "Object",
    required => 1,
    builder  => '_default_genbank_tags_ontology',
    );

sub _default_genbank_tags_ontology {
    return shift->rotiferdb->resultset('Ontology')->find_or_create
	({ name => "GenBank Annotation tags" });
}


=head2 parser_ontology

 Usage   : $parser->parser_ontology
 Function: ontology for parser-specific annotations tags
 Returns : Rotifer::DBIC::AnnotationDB::Result::Ontology
 Builder : _default_ontologies

=cut

has 'parser_ontology' => (
    is       => "rw", 
    isa      => "ArrayRef[Str]",
    lazy     => 1,
    builder  => '_default_ontologies',
    trigger  => \&_add_default_ontology,
    );

sub _default_ontologies {
    return shift->rotiferdb->resultset('Ontology')->find_or_create
	({ name => __PACKAGE__ });
}

sub _add_default_ontology {
    my ($self, $new, $old) = @_;
    unless (blessed $new) {
	my $onto = $self->rotiferdb->resultset('Ontology')->
	    find_or_create({ name => $new });
	$self->parser_ontology($new);
    }
}

=head2 tags

 Usage   : $parser->tags
 Function: annotation tags used by this parser
 Returns : Rotifer::DBIC::AnnotationDB::Result::Term
 Builder : _default_tags

 The following tags will be added to the default vocabulary
 unless present in one of the preferred ontologies:

  genbank file : sequence file in GenBank or EMBL format  

=cut

sub _default_tags {
    my ($self) = shift;

    my @tags = 
	(
	 {
	     name       => 'genbank file',
	     definition => "Text file with a list of sequence identifiers (GI or accession numbers), one per row.",
	 },
	);

    return $self->next::method(@tags, @_);
}

=head2 METHODS

=head2 load

 Title   : load
 Usage   : $full_path = $parser->load("data.gbk")
 Function: process and load data from all input files
 Returns : 
 Args    : list of file names

=cut

sub load {
    my ($self, @input) = @_;
    my $drs = $self->schema->resultset('Dataset');

    # Parsing files
    my @biodata = ();
    foreach my $file (@input) {
	my $flag = UNIVERSAL::isa($file,"GLOB") ? "-fh" : "-file";
	my $parser = Bio::SeqIO->new(%{$self->seqio_options}, $flag => $file);
	while (my $seq = $parser->next_seq) {
	    push(@biodata, $self->process_sequence_object($seq));
	}
	$parser->close;

	# Associate input with graph unless $file is a file handle
	if (!UNIVERSAL::isa($file,"GLOB") && scalar @biodata) {
	    foreach my $set (@{ $self->schema->datasets }) {
		my $rank = $drs->search_related({ term_id => $self->tags->{fasta}->id })->get_column("rank")->max;
		$set->create_related("dataset_attribute_values",
				     { term => $self->tags->{"genbank file"}, value => $file, rank => defined $rank ? ++$rank : 0 });
	}
    }

    return @biodata;
}

=head2 process_sequence_object

 Title   : process_sequence_object
 Usage   : @biodata = $self->process_sequence_object(@objs)
 Function: Bio::Seq::RichSeq to Rotifer::DBIC::AnnotationDB
 Returns : integer (number of sequences loaded)
 Args    : Bio::Seq::RichSeqI

 Note    : this method is used by load()

=cut

sub process_sequence_object {
    my ($self, $seq) = @_;

    # Tables that will be changed
    my $drs = $self->rotiferdb->resultset("Dbxref");
    my $trs = $self->rotiferdb->resultset("Term");

    # Sequence type
    my $datatype = $trs->find_or_create
	({ name        => $seq->alphabet, 
	   ontology_id => $self->sequence_ontology->id
	 });

    # Search for this sequence in RotiferDB using GI/accession
    my $ids = [ { accession => $seq->accession,
		  version   => $seq->seq_version,
		  acctype   => 'ACCESSION' },
		{ accession => $seq->primary_id,
		  acctype   => 'GI' } ];
    my @biodata = (); my %biodata;
#    foreach my $id (@ids) {
#	$xref  = $drs->search({ '-or' => $id })->search_related('node_dbxrefs')->search_related->('node');
#	@biodata = $xref->search({},{ group_by => [qw/node_id dbxref_id/] }) if ($xref->count);
#    }

    # Something in the database! Update?
    if (scalar @biodata == 1) {
	next unless ($self->update_all_histories);
    } elsif (scalar(keys %biodata) > 1) {
	croak join(" ","Sequence accessions for sequence", $seq->primary_id,
		   "were associated to more then one node! Check database...");
    }

    # Initialize node
    my $hash = {
	name        => $seq->primary_tag,
	datatype    => $datatype->id,
#	seqids      => [ @ids ],
	dbxrefs     => [],
	description => $seq->description,
	is_rawnode  => 0,
	annotation  => {},
	features    => [],
    };

    # Taxonomy
    if (!$self->process_taxonomy($hash, $seq)) {
    }

    # Copy annotations to hash
    $self->process_annotations($hash,$seq);

    # Process features
    my @needhistory = $self->process_features($hash,$seq);

    # Create biodata
    $hash = $self->create_biodata($hash);

    return $hash;
}

sub process_taxonomy {
    my ($self, $hash, $seq) = @_;
    my $crs = $self->rotiferdb->resultset('Context');

    # Taxonomy
    return 0 if (!defined $seq->species);

    # Lineage
    my $lineage = undef;
    my @lineage = $seq->species->classification;
    if (scalar(keys %{$self->preferred_taxons})) {
	@lineage = map { lc($_) } @lineage;
	@lineage = grep { exists $self->preferred_taxons->{$_} } @lineage;
	$lineage = join(">",@lineage);
    } else {
	$lineage = join(", ",@lineage);
    }

    # Context
    my $context = {
	is_organism   => 1,
	ncbi_taxon_id => $seq->species->ncbi_taxid,
	name          => $seq->species->taxon->scientific_name,
	lineage       => $lineage,
    };
    $context = $crs->update_or_create($context);

    # Add to hash
    if (defined $context) {
	$hash->{context_id} = $context->id;
    } else {
	croak "Failed to create/update organism entry for ".$hash->{name};
    }

    return 1;
}


sub process_features {
    my ($self, $hash, $seq) = @_;
    my $drs = $self->rotiferdb->resultset('Dbxref');
    my $srs = $self->rotiferdb->resultset('Biosequence');

    # Some minor hacks
    my %ignore    = map { ($_,1) } qw(mat_peptide misc_feature gene);
    my %ignoreTag = map { ($_,1) } qw(transl_table codon_start);
    my %rename    = qw(CDS protein); 

    # Features
    my @seqids = (); # Items to retrieve history for
    foreach my $generic ($seq->get_all_SeqFeatures) {
	# Datatype
	my $type = $generic->primary_tag;
	next if (exists $ignore{$type});                   # I don't like this one...
	$type = $rename{$type} if (exists $rename{$type}); # I like this one... but not the name.
	$type = $self->find_or_create
	    ({ name => $type, ontology_id => $self->sequence_ontology->id });

	# Initialize data structure for this feature
	my $feat = $hash;
	$feat = {
	    datatype   => $type->id,
	    context_id => $hash->{context_id},
	    dbxrefs    => [],
	    annotation => {},
	    is_rawnode => 0,
	    start      => $generic->start,
	    end        => $generic->end,
	    direction  => $generic->strand,
	} if ($type ne 'source');

	# Annotation collection-based interface
#	$self->process_annotation($feat,$generic) if $generic->can("annotation");

	# Tag-based interface
	foreach my $tag (sort $generic->get_all_tags) {
	    next if (exists $ignoreTag{$tag} || !$generic->has_tag($tag));
	    foreach my $value ($generic->each_tag_value($tag)) {
		# Database references
		if ($tag eq 'dbxref' && defined $value) {
		    if ($value =~ /GI:(\d+)/) {
			$feat->{name} = $1;
			push(@{$feat->{seqids}}, { accession => $1,
						   acctype   => 'GI',
						   version   => 0 });
		    } else {
			my ($dbname,$accession) = split(/:/,$value,2);
			my $acctype = "";
			push(@{$feat->{dbxrefs}}, $drs->find_or_create
			     ({ 
				 accession => $accession,
				 dbname    => $dbname,
				 version   => 0,
				 acctype   => $acctype,
			      }));
		    }
		}

		# Accession
		elsif ($tag eq 'protein_id') {
		    my ($acc, $vers) = split(/\./,$value);
		    $feat->{name} = $value unless (exists $feat->{name});
		    push(@{$feat->{seqids}}, { accession => $acc,
					       version   => 0,
					       acctype   => "ACCESSION",
			 });
		}

		# Protein sequence
		elsif ($tag eq 'translation') {
		    $value =~ s/\s+//g;
		    my $nlt = $srs->find_or_create
			({
			    term_id  => $seq->datatype,
			    residues => $value
			 });
		    $feat->{seq} = $nlt;
		}

		# Pseudogene
		elsif ($tag eq 'pseudo') {
		    $value = 1;
		}

		# Other tags
		else {
		    push(@{$hash->{annotation}{$tag}}, $value);
		}
	    } # foreach my $value ($generic->each_tag_value($tag))
	} # foreach my $tag (sort $generic->get_all_tags)

	$feat->{name} = Data::UUID->new->create_str unless (exists $feat->{name});
	push(@{$hash->{features}},$feat);
	push(@seqids, $feat) if (exists $feat->{seqids});
    } # foreach my $generic ($seq->get_all_SeqFeatures)

    return @seqids;
}

__PACKAGE__->meta->make_immutable;
1;
