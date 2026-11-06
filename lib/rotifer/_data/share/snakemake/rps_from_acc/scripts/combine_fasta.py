import pandas as pd

# Read all input .pkl files and collect into list
dfs = []

for pkl_file in snakemake.input:
    df = pd.read_pickle(pkl_file)
   # sample_name = pkl_file.split("/")[-1].replace(".pkl", "")
   # df["sample"] = sample_name  # optional, useful for tracking
    dfs.append(df)

# Concatenate into single DataFrame
combined_df = pd.concat(dfs, ignore_index=True)

# Save to pickle
combined_df.to_pickle(snakemake.output[0])
