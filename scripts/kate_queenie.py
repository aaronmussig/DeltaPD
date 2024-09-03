

import re


def filter_by_length():
    path = '/Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps.fna'
    path_fna_out = '/Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_filtered.fna'
    path_mta_out = '/Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_metadata.tsv'

    re_filter = re.compile(r'>(.+)\n(.+)')
    re_headers = re.compile(r'(.+?) (.+?) \[locus_tag=(.+?)\] \[location=(\d+)\.\.(\d+)\] \[ssu_len=(\d+)\] \[contig_len=(\d+)\]')
    with open(path, 'r') as f:
        hits = re_filter.findall(f.read())
    out_fasta = dict()
    out_meta = list()
    ssu_lengths = list()
    for header, seq in hits:

        hits_headers = re_headers.match(header)
        gid = hits_headers.group(1)
        tax = hits_headers.group(2)
        locus_tag = hits_headers.group(3)
        loc_from = hits_headers.group(4)
        loc_to = hits_headers.group(5)
        ssu_len = hits_headers.group(6)
        contig_len = hits_headers.group(7)

        seq = seq.strip()

        if 'GB_GCA_015062235.1' in gid:
            print()

        if len(seq) < 1000:
            if 'g__CAZU01' in tax:
                print(f'Filtered out g__CAZU01! {len(seq)}')
            continue

        assert gid not in out_fasta
        out_fasta[gid] = seq

        out_meta.append(f'{gid}\t{locus_tag}\t{tax}\t{loc_from}\t{loc_to}\t{contig_len}\t{ssu_len}')

        ssu_lengths.append(len(seq))

    with open(path_fna_out, 'w') as f:
        for gid, seq in out_fasta.items():
            f.write(f'>{gid}\n{seq}\n')

    out_meta = sorted(out_meta, key=lambda x: x[2])
    with open(path_mta_out, 'w') as f:
        f.write('sequence_id\tlocus_tag\ttaxonomy\tlocation_from\tlocation_to\tcontig_length\tssu_length\n')
        f.write('\n'.join(out_meta))

    print(f'Originally there were {len(hits):,} SSU sequences.')
    print(f'This has been filtered down to {len(out_fasta):,} sequences.')

    print(f'Mean ssu length: {sum(ssu_lengths) / len(ssu_lengths):.2f}')
    print(f'Min ssu length: {min(ssu_lengths)}')
    print(f'Max ssu length: {max(ssu_lengths)}')
    print(f'Median ssu length: {sorted(ssu_lengths)[len(ssu_lengths) // 2]}')

    return


def main():


    filter_by_length()


    return



if __name__ == '__main__':
    main()


"""
scp /Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_filtered.fna uqamussi@hurley.ace.uq.edu.au:/srv/home/uqamussi/projects/deltapd/queenie_kate

conda activate sina-1.7.2
sina --in=/srv/home/uqamussi/projects/deltapd/queenie_kate/bac120_ssu_reps_filtered.fna --out=/srv/home/uqamussi/projects/deltapd/queenie_kate/bac120_ssu_reps_filtered_sina.fna --fasta-write-dna --db=/srv/db/silva/138.1/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --turn


sina --in=/tmp/am/bac_ssu.fna --out=/tmp/am/bac_ssu_sina.fna --fasta-write-dna --db=/srv/db/silva/138.1/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --turn

sina --in=/tmp/am/arc_ssu.fna --out=/tmp/am/arc_ssu_sina.fna --fasta-write-dna --db=/srv/db/silva/138.1/SILVA_138.1_SSURef_NR99_12_06_20_opt.arb --turn

zcp /srv/home/uqamussi/projects/deltapd/queenie_kate/bac120_ssu_reps_filtered_sina.fna .

then run trim_msa.py

scp /Users/aaron/phd/DeltaPDNew/data/bac120_ssu_reps_filtered_sina_trim_manual.fna uqamussi@hurley.ace.uq.edu.au:/srv/home/uqamussi/projects/deltapd/queenie_kate

conda activate trimal-1.4.1
trimal -in bac120_ssu_reps_filtered_sina_trim_manual.fna -out bac120_ssu_reps_filtered_sina_trim_trimal.fasta -keepheader -automated1

trimal -in arc_ssu_sina_trim.fna -out arc_ssu_sina_trim_trimal.fasta -keepheader -automated1



conda activate fasttree-2.1.11
export OMP_NUM_THREADS=10
/usr/bin/time -v FastTreeMP -nt -out bac120_tree.tree bac120_ssu_reps_filtered_sina_trim_manual.fna > foo.log



conda activate fasttree-2.1.11
export OMP_NUM_THREADS=40
/usr/bin/time -v FastTreeMP -nt -out arc_ssu_sina_trim_min1200.tree arc_ssu_sina_trim_min1200.fna  > foo.log



"""