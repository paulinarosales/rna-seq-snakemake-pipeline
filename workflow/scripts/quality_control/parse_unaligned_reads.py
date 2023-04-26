# Import packages
import os

# File paths
## Input
raw_reads_path = snakemake.input[0]

## Output
parsed_reads_path = snakemake.output[0]


# Script
filename = os.path.basename(raw_reads_path)
with open(raw_reads_path, 'r') as file_in:
    with open(parsed_reads_path, 'w') as file_out:
        for i, row in enumerate(file_in, start=1):
            row = row.strip()
            count, __, seq = row.partition(' ')
            file_out.write(f'>{filename}|Unaligned Top {i}|n={count}\n')
            file_out.write(seq)
            file_out.write('\n')