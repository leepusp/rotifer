--
-- ClickHouse schema for UniProt's idmapping.dat
--
-- idmapping.dat is a tab separated file with three columns and no
-- header, described in the README of
-- ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/
--
--   1. UniProtKB-AC : the UniProtKB accession
--   2. ID_type      : name of the cross-referenced database
--   3. ID           : the identifier in that database
--
-- The 2026_01 release is 90,236,609,999 bytes and holds
-- 2,647,104,040 rows, spread over 100 distinct values of ID_type and
-- 203,130,941 UniProtKB accessions, so the layout below is built
-- around the two lookups that matter:
--
--   forward : "which identifiers does this UniProtKB accession have?"
--             answered by the table's own sorting key
--   reverse : "which UniProtKB accession has this identifier?"
--             answered by the by_id projection, a second copy of the
--             data sorted by id, which ClickHouse picks automatically
--
-- Measured on a 50 million row sample of the 2026_01 release, this
-- layout stores 4.56 GiB of raw text in 599 MiB on disk, 182 MiB for
-- the base table and 417 MiB for the projection, which extrapolates
-- to roughly 31 GiB for the whole file. Sorting by accession makes
-- that column compress 65-fold, and the projection is the larger
-- half precisely because sorting by id gives up that advantage.
--
-- Placeholders are filled in by rotifer.db.uniprot.clickhouse:
-- {dbname}, {table} and {release}.
--

CREATE DATABASE IF NOT EXISTS {dbname}
;

CREATE TABLE IF NOT EXISTS {dbname}.{table}
(
    `accession` String
        COMMENT 'UniProtKB accession, column 1 of idmapping.dat'
        CODEC(ZSTD(3)),

    `id_type` LowCardinality(String)
        COMMENT 'Name of the cross-referenced database, column 2 of idmapping.dat'
        CODEC(ZSTD(1)),

    `id` String
        COMMENT 'Identifier in the cross-referenced database, column 3 of idmapping.dat'
        CODEC(ZSTD(3)),

    `release` LowCardinality(String) DEFAULT '{release}'
        COMMENT 'UniProt release this row was loaded from, e.g. 2026_01'
        CODEC(ZSTD(1)),

    -- Reverse lookups that also name the database they search are
    -- common enough to deserve a cheap skipping index on the base
    -- table, so they can be answered without the projection.
    INDEX idx_id_type id_type TYPE set(0) GRANULARITY 4,

    PROJECTION by_id
    (
        SELECT accession, id_type, id, release
        ORDER BY id, id_type, accession
    )
)
ENGINE = MergeTree
PARTITION BY release
ORDER BY (accession, id_type, id)
SETTINGS index_granularity = 8192
;

--
-- Notes on loading
--
-- rotifer.db.uniprot.clickhouse.IdMappingCursor.load() runs the
-- equivalent of the command below, which streams the flat file
-- straight into the server and stamps every row with the release:
--
--   clickhouse client --query "
--     INSERT INTO {dbname}.{table}
--     SELECT c1, c2, c3, '{release}'
--     FROM input('c1 String, c2 String, c3 String')
--     FORMAT TabSeparated" < idmapping.dat
--
-- The projection doubles the work of the initial load. To load as
-- fast as possible, drop it from the CREATE TABLE above, load the
-- data, then build it in the background with:
--
--   ALTER TABLE {dbname}.{table}
--     ADD PROJECTION by_id (SELECT accession, id_type, id, release ORDER BY id, id_type, accession);
--   ALTER TABLE {dbname}.{table} MATERIALIZE PROJECTION by_id;
--
-- Because the table is partitioned by release, an obsolete release
-- is removed in one atomic, near instantaneous operation:
--
--   ALTER TABLE {dbname}.{table} DROP PARTITION '2024_06';
--
