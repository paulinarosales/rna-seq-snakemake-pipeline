from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
import pandas as pd
import io

def fix_blast_result(result):
    result_text = result.read()
    # Remove this weird "CREATE_VIEW" that sometimes appears
    result_text = result_text.replace('CREATE_VIEW', '')

    f = io.StringIO(result_text)
    return f

def query_ncbi(fasta_file):

    with open(fasta_file) as file_:
        fasta_string = file_.read()

    result = NCBIWWW.qblast('blastn', 'nt', fasta_string)
    result = fix_blast_result(result)

    records = list(NCBIXML.parse(result))

    return records

def best_match(record):

    best = None
    best_e = None

    for alignment in record.alignments:
        for hsp in alignment.hsps:

            e_score = hsp.expect
            if best_e is None or e_score < best_e:
                best = alignment.title
                best_e = e_score

    return best, best_e


def process_ncbi_results(records):

    ans = []
    for record in records:

        query = record.query
        bm, e_score = best_match(record)

        ans.append([query, bm, e_score])

    ans = pd.DataFrame(ans, columns=['Sequence', 'Best match', 'E value'])

    return ans


def process_sequences(fasta_file, bam_filename, output_tsv):

    records = query_ncbi(fasta_file)

    ans_df = process_ncbi_results(records)

    ans_df['Filename'] = bam_filename

    ans_df = ans_df[['Filename', 'Sequence', 'Best match', 'E value']]

    ans_df.to_csv(output_tsv, sep='\t', index=False)

# ------- Own additions to make it work with the new wildcards ----------------

_sample_type = snakemake.wildcards.sample_type
_treatment = snakemake.wildcards.treatment
_bio_rep = snakemake.wildcards.bio_rep

identifier = f'{_sample_type}_{_treatment}_Bio_rep_{_bio_rep}'

process_sequences(str(snakemake.input.unalignedFA), identifier, snakemake.output.blastTSV)
#process_sequences(str(snakemake.input.unalignedFA), snakemake.wildcards.sample, snakemake.output.tsv)
