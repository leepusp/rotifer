-- CREATE TEMPORARY TABLE queries (id TEXT, uuid TEXT);
-- INSERT INTO queries VALUES ('FD01876845_01788','{uuid}');
-- INSERT INTO queries VALUES ('FD01876175_02260','{uuid}');
-- INSERT INTO queries VALUES ('NP_448379.1','{uuid}');
-- CREATE INDEX IF NOT EXISTS uidx ON queries (uuid);
SELECT
  -SUM(iid) OVER (ORDER BY representative) as id, w.ipg_source,
  w.nucleotide, w.start, w.end as stop, w.strand, w.pid, w.description,
  w.ipg_organism, w.strain, w.assembly,
  SUM(1) OVER (ORDER BY representative, is_query DESC) as 'order',
  w.is_query, w.representative
FROM (
  SELECT
    CASE
      WHEN representative = LAG(representative) OVER (ORDER BY representative, is_query DESC)
      THEN 0
      ELSE 1
    END AS iid,
    '{path}' as ipg_source,
    t.*
  FROM (
      SELECT f.nucleotide, f.start, f.end, f.strand,
             f.pid, f.product as description,
             f.organism as ipg_organism, 'none' as strain,
             f.assembly,
             z.is_query, z.representative
      FROM (
       SELECT
         c2.{identical_column} as representative, c2.pid,
         max(CASE WHEN c1.pid = c2.pid THEN 1 ELSE 0 END) as is_query
       FROM
         queries as q
         inner join {identical} as c1 on (q.id == c1.pid)
         inner join {identical} as c2 using ({identical_column})
       WHERE q.uuid = '{uuid}'
       GROUP BY representative, c2.pid
      ) as z inner join features as f on (f.pid = z.pid)
  ) as t
  ORDER BY representative, is_query DESC, assembly, nucleotide
) as w
ORDER BY id DESC, 'order'
