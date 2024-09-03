import re
from collections import defaultdict
from pathlib import Path
from tqdm import tqdm

GENOME_DIRS = Path('/srv/db/gtdb/genomes/ncbi/release220/genome_dirs.tsv')

RE_SEQ = re.compile(r'>(.+)\n(.+)')


def read_genome_dirs():
    out = dict()
    with GENOME_DIRS.open() as f:
        for line in f.readlines():
            gid, path, _ = line.strip().split('\t')
            if gid.startswith('GCA_'):
                gid = f'GB_{gid}'
            elif gid.startswith('GCF_'):
                gid = f'RS_{gid}'
            else:
                raise Exception("???")
            out[gid] = Path(path)
    return out


def collect_ssu(ssu_path):
    out = dict()
    if ssu_path.exists():
        with ssu_path.open() as f:
            hits = RE_SEQ.findall(f.read())
        for header, seq in hits:
            out[header] = seq.replace('\n', '')
    return out

def get_gids_from_taxfile(path):
    out = dict()
    with path.open() as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            out[gid] = tax
    return out


"""
>GB_GCA_021654395.1 d__Archaea;p__Micrarchaeota;c__Micrarchaeia;o__Micrarchaeales;f__Micrarchaeaceae;g__ARM-1;s__ARM-1 sp021654395 [locus_tag=AP024486.1] [location=160070..161123] [ssu_len=1054] [contig_len=814439]
"""



def get_ssu_for_domain(d_gid_to_dir, gids, out_path, domain):

    with out_path.open('w') as f:
        for gid, taxonomy in tqdm(sorted(gids.items())):
            # print(gid)
            gid_dir = d_gid_to_dir[gid]

            summary_path = gid_dir / 'rna_silva_138.1' / 'ssu.hmm_summary.tsv'
            ssu_path = gid_dir / 'rna_silva_138.1' / 'ssu.fna'

            summary = read_hmm_summary_file(summary_path)
            ssu = collect_ssu(ssu_path)

            for header, seq in ssu.items():
                seqinfo = summary[header]
                if seqinfo['hmm'] == domain:
                    f.write(f'>{gid} {taxonomy} [locus_tag={header}] [location={seqinfo["start_hit"]}..{seqinfo["end_hit"]}] [ssu_len={seqinfo["ssu_len"]}] [contig_len={seqinfo["contig_len"]}]\n{seq}\n')
    return


def read_hmm_summary_file(path):
    out = dict()
    if path.exists():
        with path.open() as f:
            header = {k: i for i, k in enumerate(f.readline().strip().split('\t'))}

            for line in f.readlines():
                values = line.strip().split('\t')
                seq_id = values[header['Sequence Id']]
                hmm = values[header['HMM']]
                start_hit = values[header['Start hit']]
                end_hit = values[header['End hit']]
                ssu_len = values[header['SSU gene length']]
                contig_len = values[header['Sequence length']]

                out[seq_id] = {
                    'hmm': hmm,
                    'start_hit': int(start_hit),
                    'end_hit': int(end_hit),
                    'ssu_len': int(ssu_len),
                    'contig_len': int(contig_len)
                }

    return out


def main():
    # arc_path = Path('/tmp/am/arc_ssu.fna')
    bac_path = Path('/tmp/am/bac_ssu.fna')

    d_gid_to_dir = read_genome_dirs()

    # arc_taxonomy = get_gids_from_taxfile(Path('/tmp/am/ar53_taxonomy_r220.tsv'))
    bac_taxonomy = get_gids_from_taxfile(Path('/tmp/am/bac120_taxonomy_r220.tsv'))

    # arc_ssu = get_ssu_for_domain(d_gid_to_dir, arc_taxonomy, arc_path, domain='ar')
    bac_ssu = get_ssu_for_domain(d_gid_to_dir, bac_taxonomy, bac_path, domain='bac')

    return


if __name__ == '__main__':
    main()
