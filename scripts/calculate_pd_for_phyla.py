import re
from collections import defaultdict, deque
from pathlib import Path
from typing import Set, Dict

import dendropy


def read_tax(path: Path):
    """
    Read the GTDB taxonomy file and extract a dict keyed on phylum name, with a set of GIDs in the phylum.
    """
    out = defaultdict(set)
    with path.open() as f:
        for line in f.readlines():
            gid, tax = line.strip().split('\t')
            tax_split = tax.split(';')
            out[tax_split[1]].add(gid)
    return out


def calculate_named_node_to_leaf(d_phylum_to_node: dict[str, dendropy.Node]):
    """
    Find the total sum of branch lengths per-phylum.
    """
    out = dict()
    for phylum, node in d_phylum_to_node.items():

        # We should not skip counting the stem of the named node as this contributes to its length
        queue = deque([node])

        # Traverse the tree and sum the branch lengths
        total_len = 0
        while len(queue) > 0:
            cur_node = queue.popleft()
            total_len += cur_node.edge_length
            queue.extend(cur_node.child_nodes())

        # Store the result
        out[phylum] = total_len
    return out


def get_phylum_named_nodes(tree: dendropy.Tree, singletons: dict[str, str]):
    """
    Obtain the named node for each phylum.
    """
    # Pattern matching tax strings with a phylum in it (p__)
    re_compiled = re.compile(r'.*(p__[^\s;]+).*')

    # Go over each node in the tree
    out = dict()
    for node in tree.nodes():

        # Check if this node has a tax string in it
        if node.label and 'p__' in node.label:
            # Using the regex pattern, extract the phylum name
            re_hits = re_compiled.match(node.label)
            phylum = re_hits.group(1)
            out[phylum] = node

    # Go over the singleton cases (as they will not have had a tax string in the tree)
    for phylum, gid in singletons.items():
        if phylum not in out:

            # Find the named node that matches the representative genome id
            node = [x for x in tree.leaf_node_iter() if x.taxon.label == gid]

            # If there were no matches then this is the RED scaled tree, where shortened GIDs are used
            if len(node) == 0:
                gid = f'G{gid[7:-2]}'

            # Find the named node that matches the representative genome id
            node = [x for x in tree.leaf_node_iter() if x.taxon.label == gid]
            out[phylum] = node[0]
    return out


def calculate_for_tree(path_tree: Path, singletons: Dict[str, str]) -> Dict[str, float]:
    """
    Read the provided tree and calculate the sum of the branch lengths for each phylum.
    Returns per-phylum branch lengths.
    :param path_tree: The path to the newick tree.
    :param singletons: A dictionary of phylum to representative genome ID.
    """
    # Read the tree
    tree = dendropy.Tree.get(path=str(path_tree), schema='newick', preserve_underscores=True)

    # For each phylum, find the corresponding dendropy Node
    d_phylum_to_node = get_phylum_named_nodes(tree, singletons)

    # For each named phylum, return the total sum of all branch lengths
    named_to_leaf = calculate_named_node_to_leaf(d_phylum_to_node)
    return named_to_leaf


def get_gids_in_tree(path: Path) -> Set[str]:
    """
    Return a set of all genome IDs in the tree.
    """
    tree = dendropy.Tree.get(path=str(path), schema='newick', preserve_underscores=True)
    return {x.taxon.label for x in tree.leaf_node_iter()}


def calculate_for_domain(path_tax, path_tree, path_tree_scaled, path_out):
    # Read the taxonomy file and extract the phylum to GID mapping
    d_phylum_to_gids = read_tax(path_tax)

    # Read the tree file and read the representative genomes from the tree
    gids_in_tree = get_gids_in_tree(path_tree)

    # Filter the taxonomy mapping file to contain only those GIDs that are representatives
    d_phylum_to_gids = {k: v.intersection(gids_in_tree) for k, v in d_phylum_to_gids.items()}

    # Calculate the singletons (i.e. phyla with only one genome)
    singletons = {k: v.pop() for k, v in d_phylum_to_gids.items() if len(v) == 1}

    # For the both the scaled (RED scaled), and unscaled (GTDB) trees, calculate the sum of the branch lengths
    scaled_named_to_leaf = calculate_for_tree(path_tree_scaled, singletons)
    unscaled_named_to_leaf = calculate_for_tree(path_tree, singletons)

    # Calculate the total sum of the branch lengths for each tree
    scaled_total_sum = sum(scaled_named_to_leaf.values())
    unscaled_total_sum = sum(unscaled_named_to_leaf.values())

    # Write the results to a file
    with path_out.open('w') as f:
        f.write('phylum\tscaled_pd\tunscaled_pd\n')

        # Keep a running total of the relative percentages to ensure they sum to 100%
        total_scaled_normed = 0
        total_unscaled_normed = 0

        # Go over each phylum in alphabetical order
        for phylum in sorted(d_phylum_to_gids.keys()):
            # Calculate the sum of branch lengths for this phylum,
            # normalised by the sum of all branch lengths belonging to a phylum in the tree
            scaled_normed = 100 * scaled_named_to_leaf[phylum] / scaled_total_sum
            unscaled_normed = 100 * unscaled_named_to_leaf[phylum] / unscaled_total_sum

            # Update the sanity check counter
            total_unscaled_normed += unscaled_normed
            total_scaled_normed += scaled_normed

            # Write the results to the file
            f.write(f'{phylum}\t{scaled_normed:.4f}\t{unscaled_normed:.4f}\n')
            # print('{ taxon: "' + phylum + '", scaled: ' + str(round(scaled_normed_2, 2)) + ', unscaled: ' + str(round(unscaled_normed_2, 2)) + ', abs: ' + str(abs_delta) + ' },')

        # Complete, this should add up to 100%
        print(f'Total scaled: {total_scaled_normed:.4f}')
        print(f'Total unscaled: {total_unscaled_normed:.4f}')


def main():
    """
    Configure the program here.
    """

    path_tax = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53_taxonomy.tsv')
    path_tree = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/ar53.tree')
    path_tree_scaled = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/gtdb_r220_ar53_curation_scaled.tree')
    path_out = Path(f'/tmp/phylum_pd_ar53.tsv')

    path_tax = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/bac120_taxonomy.tsv')
    path_tree = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/bac120.tree')
    path_tree_scaled = Path(f'/Users/aaron/phd/DeltaPDNew/data/gtdb_r220/gtdb_r220_bac120_curation_scaled.tree')
    path_out = Path(f'/tmp/phylum_pd_bac120.tsv')

    calculate_for_domain(path_tax, path_tree, path_tree_scaled, path_out)

    return


if __name__ == '__main__':
    main()
