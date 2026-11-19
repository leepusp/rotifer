from rotifer.db import ncbi
from rotifer.devel.alpha import gian_func as gf
import pandas as pd
input_file = snakemake.input["acc"]
ipg_pickle= snakemake.input["ipg"]
output_file = snakemake.output["out_file"]
batch_size = snakemake.config.get("batch_size", 100)
after = int(snakemake.config.get("after", 5))
before =int(snakemake.config.get("before", 5))

with open(input_file) as f:
    pid = [line.strip() for line in f if line.strip()]

gnc = ncbi.GeneNeighborhoodCursor(batch_size=batch_size, after=after, before=before)
ipg = pd.read_pickle(ipg_pickle)
neighborhood_df = gnc.fetchall(pid,ipgs=ipg.query('superkingdom !="Eukaryota"'))
euk = ipg.query('superkingdom=="Eukaryota"').rename({'position':'block_id', 'stop':'end', 'description':'product'}, axis=1).eval('query =1').eval('type ="CDS"')
euk = euk.reindex(columns=neighborhood_df.columns)
neighborhood_df = pd.concat([neighborhood_df, euk], ignore_index=True)
#try:
#    neighborhood_df = gf.single_tax2ndf(neighborhood_df)
#except ValueError:
#    print(' Putative missing taxID, updating ETE3 taxonomoy database')
#    from ete4 import NCBITaxa
#    ncbi = NCBITaxa()
#    ncbi.update_taxonomy_database()
#    print('ETE3 taxonomy database updated')
#    #try: removed this try is breaking the code!
#    neighborhood_df = gf.single_tax2ndf(neighborhood_df)
    


neighborhood_df.to_pickle(output_file)
