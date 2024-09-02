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


def main():
    path_ref = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_r220.tree')
    path_qry = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/non_bs.tree')
    path_out = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/ar53_non_bs_metadata.tsv')

    path_extra = Path('/Users/aaron/phd/DeltaPD/examples/ar53_r220_ssu/metadata.tsv')

    extra = read_extra(path_extra)
    ref_taxa = get_taxa_from_tree(path_ref)
    qry_taxa = get_taxa_from_tree(path_qry)

    lines = [('sequence_id', 'genome_id', 'contig_length')]
    for taxon in sorted(qry_taxa):
        contig_len = extra[taxon]
        lines.append((taxon, taxon, contig_len))

    with path_out.open('w') as f:
        for line in lines:
            f.write('\t'.join(map(str, line)) + '\n')

    return



if __name__ == '__main__':
    main()