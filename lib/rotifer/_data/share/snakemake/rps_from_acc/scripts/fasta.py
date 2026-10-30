from pandas.compat import numpy
from rotifer.devel.beta import sequence as rdbs
import pandas as pd
input_file = snakemake.input["acc"]
output_file = snakemake.output["out_file"]

sequence = rdbs.sequence(input_file)
sequence.to_pickle(output_file)
