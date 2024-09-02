import re
from collections import defaultdict
import numpy as np

from tqdm import tqdm


def get_hits(path):
    headers = list()
    seqs = list()
    with open(path) as f:
        hits = re.findall(r'>(.+)\n(.+)', f.read())
        for gid, seq in hits:
            headers.append(gid.strip())
            seqs.append(seq.strip())

    msa_length = len(seqs[0])
    seq_count = len(seqs)
    array = np.empty((seq_count, msa_length), dtype=bytes)

    print('Populating the array')
    allowed = {'A', 'T', 'G', 'C', '-', 'N'}
    for idx, seq in tqdm(enumerate(seqs), total=len(seqs)):
        seq_list = list()
        for char in seq:
            if char in allowed:
                seq_list.append(char)
            else:
                seq_list.append('N')
        array[idx] = np.array(seq_list, dtype=bytes)
    return headers, array





def main():
    max_prop_gaps_in_column = 0.75
    max_prop_gaps_in_sequence = 0.75

    input_path = f'/Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_filtered_sina.fna'
    output_path = f'/Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_filtered_sina_trim_manual.fna'

    # Read the MSA
    print('Reading MSA')
    headers, array = get_hits(input_path)

    # Get the count of each amino acid in each column
    print('Performing column-wise filtering')
    cols_to_keep = list()
    n_cols = array.shape[1]
    n_seqs = array.shape[0]
    for col_idx in tqdm(range(n_cols), total=n_cols):
        col = array[:, col_idx]

        # Get the counts
        col_values, col_counts = np.unique(col, return_counts=True)
        d_value_to_prop = {v.decode(): (int(c) / n_seqs) for v, c in zip(col_values, col_counts)}
        all_values_in_col = set(d_value_to_prop.keys())
        res_in_col = all_values_in_col - {'-', 'N'}
        prop_gaps_in_column = d_value_to_prop.get('-', 0) + d_value_to_prop.get('N', 0)

        # Exclude columns that have one or fewer residues
        if len(res_in_col) <= 1:
            continue

        # Exclude columns with too many gaps present
        if prop_gaps_in_column >= max_prop_gaps_in_column:
            continue

        # Otherwise, keep the column
        cols_to_keep.append(col_idx)

    # Filter the array on the columns
    print(f'Kept {len(cols_to_keep):,} columns out of {n_cols:,} ({(len(cols_to_keep) / n_cols):.2%})')
    array = array[:, cols_to_keep]

    # Filter on sequences
    print('Filtering based on sequence composition')
    seqs_to_keep = list()
    for seq_idx in tqdm(range(n_seqs), total=n_seqs):
        cur_seq = array[seq_idx]
        seq_values, seq_counts = np.unique(cur_seq, return_counts=True)
        d_value_to_prop = {v.decode(): (int(c) / n_cols) for v, c in zip(seq_values, seq_counts)}
        prop_gaps_in_seq = d_value_to_prop.get('-', 0) + d_value_to_prop.get('N', 0)

        # If there are too many gaps, ignore it
        if prop_gaps_in_seq >= max_prop_gaps_in_sequence:
            continue

        # Otherwise, keep it
        seqs_to_keep.append(seq_idx)


    print('Writing the output')
    n_excluded = 0
    with open(output_path, 'w') as f:

        for header, seq_idx in tqdm(zip(headers, range(n_seqs)), total=n_seqs):
            cur_seq = array[seq_idx]
            seq_values, seq_counts = np.unique(cur_seq, return_counts=True)
            d_value_to_prop = {v.decode(): (int(c) / n_cols) for v, c in zip(seq_values, seq_counts)}
            prop_gaps_in_seq = d_value_to_prop.get('-', 0) + d_value_to_prop.get('N', 0)

            # If there are too many gaps, ignore it
            if prop_gaps_in_seq >= max_prop_gaps_in_sequence:
                n_excluded += 1
                continue

            # Otherwise, keep it
            f.write(f'>{header}\n{cur_seq.tostring().decode()}\n')

    print(f'Excluded {n_excluded:,}/{n_seqs} sequences')


if __name__ == '__main__':
    main()
