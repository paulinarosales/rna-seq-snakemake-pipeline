# Define packages
import pandas as pd

# File paths
## Input
file_paths = snakemake.input['subdata']

## Output
output_path = snakemake.output[0]

# Script
df = []
for subdata_file in file_paths:
    subdf = pd.read_csv(subdata_file, sep='\t')
    df.append(subdf)

df = pd.concat(df, ignore_index=True)
df.to_csv(output_path, sep='\t', index=False)