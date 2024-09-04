from pathlib import Path
import re
RE_HEADER = re.compile(r'([GBRS]{2}_[GCAF]{3}_\d{9}\.\d) (.+?) \[locus_tag=(.+?)\] \[location=(\d+?)\.\.(\d+?)\] \[ssu_len=(\d+?)\] \[contig_len=(\d+)\]')
RE_SEQ = re.compile(r'>(.+)\n(.+)')
"""
Create some variations on max length etc for a given alignment
"""
def main():

    MIN_LENGTH = 1000

    path_ssu = Path('/tmp/am/bac_ssu_sina_trim.fna')

    with path_ssu.open() as f:
        hits = RE_SEQ.findall(f.read())

    with open(f'/tmp/am/bac_ssu_sina_trim_min{MIN_LENGTH}.fna', 'w') as f:
        n_filtered = 0
        for hit, seq in hits:
            seq = seq.replace('\n', '')
            num_nucleotides = len(seq) - seq.count('N') + seq.count('-')
            if num_nucleotides < MIN_LENGTH:
                n_filtered += 1
                print(f'Filtering {hit} with {num_nucleotides:,} nucleotides')
                continue
            header_hits = RE_HEADER.match(hit)
            gid = header_hits.group(1)
            tax = header_hits.group(2)
            locus_tag = header_hits.group(3)
            loc_from = header_hits.group(4)
            loc_to = header_hits.group(5)
            ssu_len = header_hits.group(6)
            contig_len = header_hits.group(7)

            f.write(f'>{gid}__{locus_tag}\n{seq}\n')

    print(f'Filtered {n_filtered:,} sequences')


    return


if __name__ == "__main__":
    main()
