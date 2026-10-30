from pandas.compat import numpy
from rotifer.db import ncbi
import pandas as pd
input_file = snakemake.input["acc"]
output_file = snakemake.output["out_file"]
batch_size = snakemake.config.get("ipg_batch_size",60)
thread = snakemake.config.get("ipg_thread",3)
print(thread)
print(batch_size)
with open(input_file) as f:
    pid = [line.strip() for line in f if line.strip()]

ic = ncbi.IPGCursor(batch_size=batch_size,  threads=thread)
ipgtable = ic.fetchall(pid)
#i1 = ipgtable.query('pid in @pid')
#missing = list(set(pid) - set(i1.pid))
#i2 = ipgtable.query('representative in @missing')
#i2.pid = i2.representative
#ipg = pd.concat([i1,i2]).drop_duplicates('pid')
ipgtable.to_pickle(output_file)
