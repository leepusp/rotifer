-----------------------------------------------------------------------
-- Main tables: permanent storage
-----------------------------------------------------------------------

-- Table: assembly_summary
-- DROP TABLE assembly_summary;
CREATE TABLE assembly_summary
(
    assembly character varying NOT NULL,
    bioproject character varying(64) ,
    biosample character varying(64) ,
    wgs_master character varying(64) ,
    refseq_category character varying(64) ,
    taxid integer,
    species_taxid integer,
    organism_name text ,
    infraspecific_name text ,
    isolate text ,
    version_status character varying(64) ,
    assembly_level character varying(64) ,
    release_type character varying(64) ,
    genome_rep character varying(64) ,
    seq_rel_date character varying(64) ,
    asm_name text ,
    submitter text ,
    gbrs_paired_asm character varying(64) ,
    paired_asm_comp character varying(64) ,
    ftp_path text ,
    excluded_from_refseq text ,
    relation_to_type_material text ,
    obsolete boolean NOT NULL,
    lineage text,
    classification text,
    CONSTRAINT assembly_summary_pkey PRIMARY KEY (assembly),
);

-- DROP INDEX asummary_assembly_idx;
CREATE UNIQUE INDEX asummary_assembly_idx
    ON assembly_summary USING btree
    (assembly );

-- DROP INDEX asummary_species_taxid_idx;
CREATE INDEX asummary_species_taxid_idx
    ON assembly_summary USING btree
    (species_taxid);

-- DROP INDEX asummary_taxid_idx;
CREATE INDEX asummary_taxid_idx
    ON assembly_summary USING btree
    (taxid);

-- Table: ftable
-- DROP TABLE ftable;
CREATE TABLE ftable
(
    feature character varying(64) ,
    feature_class character varying(64) ,
    assembly character varying(64) NOT NULL,
    assembly_unit character varying(64) ,
    seq_type character varying(64) ,
    chromosome character varying(64) ,
    genomic_accession character varying(64) NOT NULL,
    gstart integer,
    gend integer,
    strand character(1) NOT NULL DEFAULT '+'::bpchar,
    product_accession character varying(64) ,
    non_redundant_refseq character varying(64) ,
    related_accession character varying(64) ,
    product_name text ,
    symbol character varying(64) ,
    geneid character varying(64) ,
    locus_tag character varying(64) ,
    feature_length integer,
    product_length integer,
    product_attributes text ,
    class_order   integer NOT NULL,
    feature_order integer NOT NULL,
    genomic_order integer NOT NULL,
    feature_class_modified boolean,
    locus_tag_modified boolean,
    ftable_id bigint NOT NULL,
    CONSTRAINT ftable_id_pk PRIMARY KEY (ftable_id),
    CONSTRAINT ftable_assembly_fkey FOREIGN KEY (assembly)
        REFERENCES assembly_summary (assembly) MATCH SIMPLE
        ON UPDATE NO ACTION
        ON DELETE NO ACTION
);

-- DROP INDEX ftable_assembly_idx;
CREATE INDEX ftable_assembly_idx
    ON ftable USING btree
    (assembly );

-- DROP INDEX ftable_genomic_accession_idx;
CREATE INDEX ftable_genomic_accession_idx
    ON ftable USING btree
    (genomic_accession );

-- DROP INDEX ftable_feature_idx;
CREATE INDEX ftable_feature_idx
    ON ftable USING btree
    (feature );

-- DROP INDEX ftable_feature_class_idx;
CREATE INDEX ftable_feature_class_idx
    ON ftable USING btree
    (feature_class );

-- DROP INDEX ftable_feature_class_idx;
CREATE INDEX ftable_agfc_idx
    ON ftable USING btree
    (assembly, genomic, feature, feature_class);

-- DROP INDEX ftable_pacc_idx;
CREATE INDEX ftable_pacc_idx
    ON ftable USING btree
    (product_accession );

-- DROP INDEX ftable_locus_idx;
CREATE INDEX ftable_locus_idx
    ON ftable USING btree
    (locus_tag );

-- DROP INDEX ftable_class_order_idx;
CREATE INDEX ftable_class_order_idx
    ON ftable USING brin
    (class_order);

-- DROP INDEX ftable_feature_order_idx;
CREATE INDEX ftable_feature_order_idx
    ON ftable USING brin
    (feature_order);

-- DROP INDEX ftable_genomic_order_idx;
CREATE INDEX ftable_genomic_order_idx
    ON ftable USING brin
    (genomic_order);

-----------------------------------------------------------------------
-- Functions
-----------------------------------------------------------------------

-- Find upstream neighbors of the same kind as the query
DROP FUNCTION IF EXISTS first_upstream_neighbor;
CREATE OR REPLACE FUNCTION first_upstream_neighbor(in gacc varchar(64), in gorder integer, in n integer)
RETURNS integer AS $$
    SELECT min(genomic_order)
    FROM (
        SELECT genomic_accession, genomic_order
        FROM ftable AS ft
        WHERE
              ft.genomic_accession = $1
              AND ft.genomic_order <= $2
              AND ft.feature = (SELECT feature FROM ftable WHERE genomic_accession = $1 AND genomic_order = $2)
              AND ft.feature_class = (SELECT feature_class FROM ftable WHERE genomic_accession = $1 AND genomic_order = $2)
        ORDER BY genomic_order DESC
        LIMIT $3+1
    ) as ftbefore
$$
LANGUAGE SQL;

-- -- Function to find downstream neighbors of the same kind as the query
DROP FUNCTION IF EXISTS last_downstream_neighbor;
CREATE OR REPLACE FUNCTION last_downstream_neighbor(in gacc varchar(64), in gorder integer, in n integer)
RETURNS integer AS $$
    SELECT max(genomic_order)
    FROM (
        SELECT genomic_order
        FROM ftable AS ft
        WHERE
              ft.genomic_accession = $1
              AND ft.genomic_order >= $2
              AND ft.feature = (SELECT feature FROM ftable WHERE genomic_accession = $1 AND genomic_order = $2)
              AND ft.feature_class = (SELECT feature_class FROM ftable WHERE genomic_accession = $1 AND genomic_order = $2)
        ORDER BY genomic_order ASC
        LIMIT $3+1
    ) as ftbefore
$$
LANGUAGE SQL;
