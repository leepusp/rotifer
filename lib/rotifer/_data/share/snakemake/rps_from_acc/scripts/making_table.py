import pandas as pd

# Read all input .pkl files and collect into list
cluster = snakemake.input.cluster
pfam = snakemake.input.pfam
profiledb = snakemake.input.profiledb
arch = snakemake.input.arch

output = pd.read_csv(cluster, sep='\t', names=['c80e3', 'pid'])[['pid','c80e3']]
output = output.merge(
        pd.read_csv(arch,
                    sep='\t',
                    names=['pid', 'arch'],
                    usecols=[0, 1]
                    ),
        how='left')

output = output.merge(
        pd.read_csv(profiledb,
                    sep='\t',
                    names=['pid', 'profiledb'],
                    usecols=[0, 1]
                    ),
        how='left')

output = output.merge(
        pd.read_csv(pfam,
                    sep='\t',
                    names=['pid', 'pfam'],
                    usecols=[0, 1]
                    ),
        how='left')


output.to_pickle('final.pkl')
