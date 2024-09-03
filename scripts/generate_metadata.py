from pathlib import Path
import dendropy
from deltapd.tree import Tree


def get_taxa_from_tree(path: Path) -> set[str]:
    tree = Tree(path)
    return {x.label for x in tree.tree.taxon_namespace}

def read_extra(path: Path):
    out = dict()
    with path.open() as f:
        for line in f.readlines():
            gid, _, _, contig_len = line.strip().split('\t')
            out[gid] = int(contig_len)
    return out

def read_metadata(path):
    out = dict()
    with path.open() as f:
        header = f.readline().strip().split('\t')
        rep_col = header.index('gtdb_genome_representative')
        for line in f.readlines():
            data = line.strip().split('\t')
            gid = data[0]
            rep = data[rep_col]
            out[gid] = rep
    return out


def main():
    path_ref = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53.tree')
    path_meta = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/ar53_metadata_r220.tsv')
    path_qry = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200.tree')
    path_out = Path('/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53/arc_ssu_sina_trim_min1200_metadata.tsv')

    d_gid_to_rep = read_metadata(path_meta)
    qry_taxa = get_taxa_from_tree(path_qry)

    lines = [('sequence_id', 'genome_id', 'contig_length')]
    for taxon in sorted(qry_taxa):
        contig_len = 1234
        rep_gid = d_gid_to_rep[taxon.split('__')[0]]
        lines.append((taxon, rep_gid, contig_len))

    with path_out.open('w') as f:
        for line in lines:
            f.write('\t'.join(map(str, line)) + '\n')

    return



if __name__ == '__main__':
    main()