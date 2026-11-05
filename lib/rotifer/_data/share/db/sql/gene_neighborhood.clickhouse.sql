--
-- This code is a prototype of the code for searching gene neighborhoods
-- in our ClickHouse database (table genomes)
--
-- The final implementation will require random names for temporary tables
--

--
-- First step is to load a list of queries
--
-- Interactively for testing purposes use:
--
-- clickhouse-client -d rotifer -m
-- DROP TABLE IF EXISTS _input;
-- CREATE TABLE _input (query String, type String) Engine=Memory;
--
-- The final, non-interactive implementation, should use:
--
--DROP TEMPORARY TABLE IF EXISTS _input;
--CREATE TEMPORARY TABLE _input (query String, type String) Engine=Memory;

--
-- Load the data into _input using clickhouse-local or another client
--
-- Interactively for testing purposes use:
-- tfilter -f '$F[1]="CDS";1' /home/rfsouza/projects/ncyclases/work/20191020/h0.acc | clickhouse-client -d rotifer -q 'INSERT INTO _input FORMAT TabSeparated'
--
-- The final, non-interactive implementation should use:
-- INSERT INTO _input ('list_of_queries');

--
-- Find all target features and save the initial data to _temp1
--
DROP TEMPORARY TABLE IF EXISTS _temp1;
CREATE TEMPORARY TABLE _temp1 (
  nuc_asm       String, -- Assembly identifier
  nucleotide    String, -- Nucleotide accession
  topology      String, -- Chromosome topology: circular or linear
  type          String, -- Feature type. E.g. CDS, tRNA, etc
  query         String, -- Query identifier (protein ID, gene name, locus tag, etc.)
  feature_order Int64,  -- 1-based and type-sensitive feature number or id
  minfo         Int64,  -- Lower limit of feature_order for the block around this feature
  maxfo         Int64,  -- Upper limit of feature_order for the block around this feature
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory
;

--INSERT INTO _temp1
--WITH 10 as neighbors_up, 10 as neighbors_down
--SELECT
--  nuc_asm, nucleotide, topology, type, query, feature_order,
--  feature_order - neighbors_up as minfo, feature_order + neighbors_down as maxfo,
--  organism, taxonomy
--FROM genomes as g INNER JOIN _input as i ON (g.pid = i.query AND g.type = i.type)
--;

--
-- Collect internal_id and feature_order ranges, per nucleotide
-- This data (iimax and fomax) could already be in the 'genomes' table!
-- The minimum values should always be 1
--
DROP TEMPORARY TABLE IF EXISTS _temp2;
CREATE TEMPORARY TABLE _temp2 (
  nuc_asm       String, -- Assembly identifier
  nucleotide    String, -- Nucleotide accession
  type          String, -- Feature type. E.g. CDS, tRNA, etc
  iimin         Int64,  -- Minimum internal_id for features of this 'type' in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id for features of this 'type' in this 'nucleotide' sequence
  fomin         Int64,  -- Minimum feature_order of this 'type' in this 'nucleotide' sequence
  fomax         Int64   -- Maximum feature_order of this 'type' in this 'nucleotide' sequence
) ENGINE=Memory;

INSERT INTO _temp2
SELECT nuc_asm, nucleotide, type, iimin, iimax, fomin, fomax
FROM (
  SELECT nuc_asm, nucleotide, type, min(feature_order) AS fomin, max(feature_order) AS fomax
  FROM genomes as g INNER JOIN _temp1 as t1 ON (g.nuc_asm = t1.nuc_asm AND g.nucleotide = t1.nucleotide and g.type = t1.type)
  GROUP BY nuc_asm, nucleotide, type
) as t1 INNER JOIN (
  SELECT nuc_asm, nucleotide, min(internal_id) AS iimin, max(internal_id) AS iimax
  FROM genomes as h INNER JOIN _temp1 as t1 ON (h.nuc_asm = t1.nuc_asm AND h.nucleotide = t1.nucleotide)
  GROUP BY nuc_asm, nucleotide
) as t2 ON (t1.nuc_asm = t2.nuc_asm AND t1.nucleotide = t2.nucleotide)
;

--
-- Add limits (_temp2 data) to the table of features, sort and mark 
-- consecutive features as belonging to either the same or different
-- blocks (bchange)
--
DROP TEMPORARY TABLE IF EXISTS _temp3;
CREATE TEMPORARY TABLE _temp3 (
  nuc_asm       String, -- Assembly identifier
  nucleotide    String, -- Nucleotide accession
  topology      String, -- Chromosome topology: circular or linear
  iimin         Int64,  -- Minimum internal_id in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id in this 'nucleotide' sequence
  fomin         Int64,  -- Minimum feature_order of this 'type' in this 'nucleotide' sequence
  fomax         Int64,  -- Maximum feature_order of this 'type' in this 'nucleotide' sequence
  type          String, -- Feature type. E.g. CDS, tRNA, etc
  query         String, -- Query identifier (protein ID, gene name, locus tag, etc.)
  feature_order Int64,  -- 1-based and type-sensitive feature number or id
  minfo         Int64,  -- Lower limit of feature_order for the block around this feature
  maxfo         Int64,  -- Upper limit of feature_order for the block around this feature
  bchange       Int64,  -- Whether this feature belongs to the same block as the previous feature (boolean)
  rid           Int64,  -- Row number
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

-- A subselect is used to make sure consecutive rows are correctly sorted
INSERT INTO _temp3
WITH 0 as minimum_block_distance
SELECT nuc_asm, nucleotide, topology, iimin, iimax, fomin, fomax, type, query, feature_order, minfo, maxfo,
       CASE WHEN nuc_asm = neighbor(nuc_asm, -1) AND nucleotide = neighbor(nucleotide, -1) AND type = neighbor(type, -1) AND
	         minfo - neighbor(maxfo, -1, minfo - minimum_block_distance - 1) - 1 <= minimum_block_distance
            THEN 0 ELSE 1 END as bchange,
       rowNumberInAllBlocks() + 1 as rid,
       organism, taxonomy
FROM 
(
 WITH 0 as minimum_block_distance
 SELECT 
   t1.nuc_asm, t1.nucleotide, t1.topology, t2.iimin, t2.iimax, t2.fomin, t2.fomax, t1.type, t1.query, t1.feature_order,
   CASE
     WHEN t1.topology = 'linear'   AND t1.minfo < t2.fomin THEN t2.fomin
     WHEN t1.topology = 'circular' AND (t2.fomax - t1.maxfo) + t1.minfo - 1 <= minimum_block_distance THEN t2.fomin
     WHEN t1.topology = 'circular' AND t1.maxfo > t2.fomax THEN t1.minfo - t2.fomax
     ELSE t1.minfo
   END as minfo,
   CASE
     WHEN t1.topology = 'linear'   AND t1.maxfo > t2.fomax THEN t2.fomax
     WHEN t1.topology = 'circular' AND (t2.fomax - t1.maxfo) + t1.minfo - 1 <= minimum_block_distance THEN t2.fomax
     WHEN t1.topology = 'circular' AND t1.maxfo > t2.fomax THEN t1.maxfo - t2.fomax
     ELSE t1.maxfo
   END as maxfo,
   t1.organism, t1.taxonomy
 FROM _temp1 as t1 INNER JOIN _temp2 AS t2 ON (t1.nuc_asm = t2.nuc_asm AND t1.nucleotide = t2.nucleotide AND t1.type = t2.type)
 ORDER BY taxonomy, organism, nuc_asm, nucleotide, type, minfo
) as t
;

DROP TEMPORARY TABLE IF EXISTS _temp1;
DROP TEMPORARY TABLE IF EXISTS _temp2;

--
-- Assign initial block IDs and merge neighborhoods that are too close
--
-- Blocks are too close if the minimum number of interviening features
-- of its own kind is smaller than the minimum required. This definition
-- also encompasses overlapping blocks.
--
DROP TEMPORARY TABLE IF EXISTS _temp4;
CREATE TEMPORARY TABLE _temp4 (
  nuc_asm       String,
  nucleotide    String,
  topology      String,
  iimin         Int64,  -- Minimum internal_id in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id in this 'nucleotide' sequence
  fomin         Int64,
  fomax         Int64,
  type          String,
  query         Array(String),
  feature_order Array(Int64),
  minfo         Int64,
  maxfo         Int64,
  block_id      Int64,
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

INSERT INTO _temp4
WITH (SELECT groupArrayMovingSum(bchange) FROM (SELECT bchange, rid FROM _temp3 ORDER BY rid)) as bid
SELECT
  nuc_asm,
  nucleotide,
  any(topology) as topology,
  min(iimin) as iimin,
  max(iimax) as iimax,
  min(fomin) as fomin,
  max(fomax) as fomax,
  type,
  arraySort(groupArray(query)) as query,
  arraySort(groupArray(feature_order)) as feature_order,
  min(minfo) as minfo,
  max(maxfo) as maxfo,
  bid[rid] AS block_id,
  organism,
  taxonomy
FROM (SELECT * FROM _temp3 ORDER BY rid) as t1
GROUP BY taxonomy, organism, nuc_asm, nucleotide, type, block_id
ORDER BY taxonomy, organism, nuc_asm, nucleotide, type, minfo
;

-- Should I save memory?
-- No information is lost as data moves from _temp3 to _temp4
DROP TEMPORARY TABLE IF EXISTS _temp3;

--
-- Parse assembly ids and collect assembly statistics
--
DROP TEMPORARY TABLE IF EXISTS _temp5;
CREATE TEMPORARY TABLE _temp5 (
  adiv          String, -- Assembly section: GCF, GCA or same as nuc_asm
  aacc          String, -- Core (numeric part) of the assembly ID or same as nuc_asm
  aver          Int64,  -- Assembly version or same as nuc_asm
  acat          Int64,  -- Assembly division category
  nuc_asm       String, -- Assembly identifier
  nqfound       Int64,  -- Number of queries found in each assembly
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

INSERT INTO _temp5
SELECT
  CASE
    WHEN substring(nuc_asm,1,3) = 'GCA' OR substring(nuc_asm,1,3) = 'GCF' THEN splitByChar('_',nuc_asm)[1]
    ELSE 'GenBank'
  END as adiv,
  CASE
    WHEN substring(nuc_asm,1,3) = 'GCA' OR substring(nuc_asm,1,3) = 'GCF' THEN splitByChar('.',splitByChar('_',nuc_asm)[2])[1]
    ELSE splitByChar('.',nuc_asm)[1]
  END as aacc,
  toInt64(splitByChar('.',nuc_asm)[2]) as aver,
  CASE
    WHEN adiv = 'GCA' THEN -1
    WHEN adiv = 'GCF' THEN  1
    ELSE 0
  END as acat,
  nuc_asm,
  sum(length(query)) as nqfound,
  organism,
  taxonomy
FROM _temp4
GROUP BY taxonomy, organism, nuc_asm
;

--
-- Detect older/redundant versions of the same assembly ond/or nucleotide
--
DROP TEMPORARY TABLE IF EXISTS _assembly;
CREATE TEMPORARY TABLE _assembly (
  nass          Int64,         -- Number of related assemblies
  aid           String,        -- selected non-redundant assembly
  adiv          String,        -- Assembly section: GCF, GCA or same as nuc_asm
  aacc          String,        -- Core (numeric part) of the assembly ID or same as nuc_asm
  aver          Int64,         -- Assembly version or same as nuc_asm
  acat          Int64,         -- Assembly division category
  nuc_asm       String,        -- Assembly identifier
  nqfound       Int64,         -- Number of queries found in each assembly
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

--
-- Add number of redundant assemblies and assembly score
--
INSERT INTO _assembly
SELECT nass, aid, t5.*
FROM _temp5 as t5 INNER JOIN (
  SELECT
    aacc,
    count(distinct nuc_asm) as nass,
    groupUniqArray(adiv) as aadiv,
    CASE
      WHEN length(aadiv) = 1 THEN argMax(nuc_asm, aver)
      WHEN has(aadiv,'GCF')  THEN argMax(nuc_asm, acat*aver)
      WHEN has(aadiv,'GCA')  THEN argMin(nuc_asm, acat*aver)
      ELSE argMax(nuc_asm, aver)
    END as aid
  FROM _temp5
  GROUP BY aacc
) as t1 ON (t5.aacc = t1.aacc)
ORDER BY nass DESC, aacc, adiv DESC, aver DESC
;

DROP TEMPORARY TABLE IF EXISTS _temp5;

--
-- 1: Fix blocks that run through the origin of circular replicons 
--    by comparing coordinates of the first and last blocks (t3)
-- 2: Merge last and first blocks, if last necessary (t4)
-- 3: Fix blocks running through the origin (minfo < fomin)
--    We break these blocks coordinates in two pairs:
--     [ [minfo+fomax,fomax], [fomin,maxfo] ]
--    and expand (arrayJoin) the rows to separate the pairs (t5)
-- 4: Add block change markers (bchange) and whether the assembly
--    is valid (keep), i.e. non-redundant
--
DROP TEMPORARY TABLE IF EXISTS _temp6;
CREATE TEMPORARY TABLE _temp6 (
  nuc_asm       String,
  nucleotide    String,
  topology      String,
  iimin         Int64,  -- Minimum internal_id in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id in this 'nucleotide' sequence
  fomin         Int64,
  fomax         Int64,
  type          String,
  query         Array(String),
  feature_order Array(Int64),
  minfo         Int64,
  maxfo         Int64,
  block_id      Int64,
  origin        Int64,
  keep          Int64,
  rid           Int64,
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

INSERT INTO _temp6
SELECT
  nuc_asm, nucleotide, topology, iimin, iimax, fomin, fomax, type, query, feature_order,
  coord[1] as minfo, coord[2] as maxfo, block_id,
  origin, keep, rowNumberInAllBlocks() + 1 as rid,
  organism, taxonomy
FROM (
 SELECT
   nuc_asm, nucleotide, topology, iimin, iimax, fomin, fomax, type, query, feature_order,
   arrayJoin(CASE WHEN minfo < fomin THEN [[minfo+fomax,fomax],[fomin,maxfo]] ELSE [[minfo,maxfo]] END) as coord,
   CASE WHEN minfo < fomin THEN 1 ELSE 0 END as origin,
   CASE WHEN nuc_asm IN (SELECT DISTINCT aid FROM _assembly) THEN 1 ELSE 0 END as keep,
   block_id, organism, taxonomy
 FROM (
  SELECT
    nuc_asm, nucleotide, any(topology) as topology,
    min(iimin) as iimin, max(iimax) as iimax,
    min(fomin) as fomin, max(fomax) as fomax,
    type,
    arraySort(groupUniqArrayArray(query)) as query, arraySort(groupArrayArray(feature_order)) as feature_order,
    min(minfo) as minfo, max(maxfo) as maxfo, block_id, taxonomy, organism
  FROM (
   SELECT
     t2.nuc_asm, t2.nucleotide, t2.topology, t2.iimin, t2.iimax, t2.fomin, t2.fomax, t2.type, t2.query, t2.feature_order,
     CASE WHEN t2.block_id = t1.oldbid THEN minfo - fomax ELSE minfo       END as minfo,
     CASE WHEN t2.block_id = t1.oldbid THEN maxfo - fomax ELSE maxfo       END as maxfo,
     CASE WHEN t2.block_id = t1.oldbid THEN t1.newbid     ELSE t2.block_id END as block_id,
     t2.organism, t2.taxonomy
   FROM _temp4 as t2 LEFT JOIN (
     WITH 0 as minimum_block_distance
     SELECT nuc_asm, nucleotide, type, max(block_id) as oldbid, min(block_id) as newbid
     FROM _temp4
     WHERE topology = 'circular'
     GROUP BY nuc_asm, nucleotide, type
     HAVING count() > 1 AND max(fomax) + min(minfo) - max(maxfo) - 1 <= minimum_block_distance
   ) as t1 ON (t2.nuc_asm = t1.nuc_asm AND t2.nucleotide = t1.nucleotide AND t2.block_id = t1.oldbid)
   ORDER BY taxonomy, organism, nuc_asm, nucleotide, type, minfo
  ) as t3
  GROUP BY taxonomy, organism, nuc_asm, nucleotide, type, block_id
  ORDER BY taxonomy, organism, nuc_asm, nucleotide, type, minfo
 ) as t4
 ORDER BY keep DESC, taxonomy, organism, nuc_asm, nucleotide, type, block_id, coord[1] DESC
) as t5
ORDER BY keep DESC, taxonomy, organism, nuc_asm, nucleotide, type, block_id, minfo DESC
;

-- Should I save memory?
-- No information is lost as data moves from _temp4 to _temp6
DROP TEMPORARY TABLE IF EXISTS _temp4;

--
-- Collect internal IDs
--
DROP TEMPORARY TABLE IF EXISTS _temp7;
CREATE TEMPORARY TABLE _temp7 (
  nuc_asm       String,
  nucleotide    String,
  topology      String,
  iimin         Int64,  -- Minimum internal_id in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id in this 'nucleotide' sequence
  fomin         Int64,
  fomax         Int64,
  type          String,
  query         Array(String),
  feature_order Array(Int64),
  minfo         Int64,
  maxfo         Int64,
  minii         Int64,
  maxii         Int64,
  start         Int64,
  end           Int64,
  block_id      Int64,
  origin        Int64,
  keep          Int64,
  rid           Int64,
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

INSERT INTO _temp7
SELECT
  any(nuc_asm) as nuc_asm, any(nucleotide) as nucleotide, any(topology) as topology,
  min(iimin) as iimin, max(iimax) as iimax, min(fomin) as fomin, max(fomax) as fomax,
  any(type), groupUniqArrayArray(query) as query, any(feature_order) as feature_order,
  min(minfo) as minfo, max(maxfo) as maxfo, min(internal_id) as minii, max(internal_id) as maxii,
  min(start) as start, max(end) as end,
  any(block_id) as block_id, max(origin) as origin, max(keep) as keep, rid,
  any(organism) as organism, any(taxonomy) as taxonomy
FROM (
 SELECT
   CASE WHEN origin = 1 AND minfo = fomin THEN iimin ELSE g.internal_id     END as internal_id,
   CASE WHEN origin = 1 AND minfo = fomin THEN 1     ELSE toUInt64(g.start) END as start,
   toUInt64(g.end) as end,
   nuc_asm, nucleotide, t1.topology as topology,
   iimin, iimax, fomin, fomax, type, query, t1.feature_order as feature_order,
   minfo, maxfo, block_id, origin, keep, rid, t1.organism as organism, t1.taxonomy as taxonomy
 FROM genomes as g INNER JOIN _temp6 as t1 ON (g.nuc_asm = t1.nuc_asm) AND (g.nucleotide = t1.nucleotide) AND (g.feature_order = t1.minfo) AND (g.type = t1.type)
 UNION ALL
 SELECT
   CASE WHEN origin = 1 AND maxfo = fomax THEN iimax           ELSE g.internal_id   END as internal_id,
   toUInt64(g.start) as start,
   CASE WHEN origin = 1 AND maxfo = fomax THEN toUInt64(g.end) ELSE toUInt64(g.end) END as end,
   nuc_asm, nucleotide, t2.topology as topology,
   iimin, iimax, fomin, fomax, type, query, t2.feature_order as feature_order,
   minfo, maxfo, block_id, origin, keep, rid, t2.organism as organism, t2.taxonomy as taxonomy
 FROM genomes as g INNER JOIN _temp6 as t2 ON (g.nuc_asm = t2.nuc_asm) AND (g.nucleotide = t2.nucleotide) AND (g.feature_order = t2.maxfo) AND (g.type = t2.type)
) as t3
GROUP BY rid
ORDER BY rid
;

DROP TEMPORARY TABLE IF EXISTS _temp6;

--
-- Reset block_id
--
DROP TEMPORARY TABLE IF EXISTS _blocks;
CREATE TEMPORARY TABLE _blocks (
  nuc_asm       String,
  nucleotide    String,
  topology      String,
  iimin         Int64,  -- Minimum internal_id in this 'nucleotide' sequence
  iimax         Int64,  -- Maximum internal_id in this 'nucleotide' sequence
  fomin         Int64,
  fomax         Int64,
  type          String,
  query         Array(String),
  feature_order Array(Int64),
  minfo         Int64,
  maxfo         Int64,
  minii         Int64,
  maxii         Int64,
  start         Int64,
  end           Int64,
  block_id      Int64,
  origin        Int64,
  keep          Int64,
  rid           Int64,
  organism      String, -- Organism name
  taxonomy      String  -- Organism taxonomic data
) ENGINE=Memory;

INSERT INTO _blocks
WITH (
  SELECT arrayCumSum(groupArray(CASE WHEN block_id = neighbor(block_id,-1,-1) THEN 0 ELSE 1 END)) as bchange
  FROM (
    SELECT block_id
    FROM _temp7
    ORDER BY keep, rid
  )
) as bid
SELECT
  nuc_asm, nucleotide, topology, iimin, iimax, fomin, fomax, type, query, feature_order,
  minfo, maxfo, minii, maxii, start, end, bid[rid] as block_id, origin, keep, rid, organism, taxonomy
FROM _temp7
ORDER BY keep, rid
;

-- Should I save memory?
-- No information is lost as data moves from _temp7 to _blocks
DROP TEMPORARY TABLE IF EXISTS _temp7;

--
-- Find and store missing queries
--
DROP TEMPORARY TABLE IF EXISTS _missing;
CREATE TEMPORARY TABLE _missing (nuc_asm String, nucleotide String, query String) ENGINE=Memory;

INSERT INTO _missing
SELECT '' as nuc_asm, '' as nucleotide, query
FROM _input
WHERE query NOT IN (SELECT distinct(arrayJoin(query)) FROM _blocks)
;

INSERT INTO _missing
SELECT *
FROM (
  SELECT DISTINCT nuc_asm, nucleotide, arrayJoin(query) as xquery
  FROM _blocks
  WHERE keep = 0
  GROUP BY nuc_asm, nucleotide, xquery
) as t1
WHERE xquery NOT IN (SELECT distinct(arrayJoin(query)) FROM _blocks WHERE keep = 1)
;

--
-- Prepare list of rows to retrieve from genomes tables
--
DROP TABLE IF EXISTS _temp8;
CREATE TABLE _temp8 (
  nuc_asm       String,
  nucleotide    String,
  internal_id   String,
  query         Array(String),
  block_id      Int64,
  origin        Int64,
  rid           Int64
) ENGINE=Join(ALL,INNER, nuc_asm, nucleotide, internal_id);

INSERT INTO _temp8
SELECT
  nuc_asm, nucleotide,
  toString(arrayJoin(range(toUInt64(minii),toUInt64(maxii)+1))) as internal_id,
  query, block_id, origin, rid
FROM _blocks
WHERE keep = 1
;

--
-- Output table (unsorted)
DROP TEMPORARY TABLE IF EXISTS _temp9;
CREATE TEMPORARY TABLE _temp9 (
  nucleotide    String,
  start         Int64,
  end           Int64,
  strand        String,
  block_id      Int64,
  rid           Int64,
  query         Int64,
  pid           String,
  type          String,
  plen          String,
  locus         String,
  seq_type      String,
  assembly      String,
  gene          String,
  origin        Int64,
  topology      String,
  product       String,
  organism      String,
  taxonomy      String
) ENGINE = Memory;

--
-- Dump selected rows to stdout
--
INSERT INTO _temp9
WITH ['CDS','PSE','assembly_gap','gap','misc_RNA','misc_binding','misc_feature','mobile_element','ncRNA','operon','rRNA','regulatory','repeat_region','tRNA','tmRNA'] as included
SELECT 
  nucleotide,
  toInt64(start) as start,
  toInt64(end) as end,
  strand,
  block_id,
  rid,
  has(query,pid) as query,
  pid,
  type,
  plen,
  locus,
  seq_type,
  nuc_asm as assembly,
  gene,
  origin,
  topology,
  product,
  organism,
  taxonomy as classification
FROM genomes as g INNER JOIN _temp8 as t8 ON (g.nuc_asm = t8.nuc_asm AND g.nucleotide = t8.nucleotide AND toString(g.internal_id) = t8.internal_id)
WHERE has(included,type)
;

DROP TABLE _temp8;

SELECT *
FROM _temp9
ORDER BY taxonomy, organism, assembly, nucleotide, block_id, rid
FORMAT TSVWithNames
;

DROP TEMPORARY TABLE IF EXISTS _temp9;

--
-- Using the calculated feature_order limits to search for a specific neighborhood
-- Can I use it for all neighborhoods at once with Join?
--
--SELECT *
--FROM genomes
--WHERE
--  nuc_asm = 'GCF_000614305.2' AND nucleotide = 'NZ_JHKG01000021.1'
--  AND internal_id BETWEEN (
--       SELECT internal_id as minii
--       FROM genomes
--       WHERE nuc_asm = 'GCF_000614305.2' AND nucleotide = 'NZ_JHKG01000021.1' AND type = 'CDS' AND feature_order = 10
--      ) AND (
--       SELECT internal_id as maxii
--       FROM genomes
--       WHERE nuc_asm = 'GCF_000614305.2' AND nucleotide = 'NZ_JHKG01000021.1' AND type = 'CDS' AND feature_order = 30
--      )
--  AND type NOT IN ('gene')
--ORDER BY internal_id
--;

