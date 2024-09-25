from pathlib import Path
from typing import Collection, Tuple

import dendropy


class DeltaPdTree:

    def __init__(self, path: Path):
        self.path = path
        self.tree = self.read(path)

    @staticmethod
    def read(path: Path) -> dendropy.Tree:
        return dendropy.Tree.get(path=path, schema='newick', preserve_underscores=True)

    def get_nodes_and_edges_for_deltapd(self):
        d_node_to_idx = dict()
        taxa = set()
        edges = set()

        # Collect the nodes and taxa
        for node in self.tree.postorder_node_iter():

            # Store this node
            cur_node_id = d_node_to_idx.get(node)
            if cur_node_id is None:
                cur_node_id = len(d_node_to_idx)
                d_node_to_idx[node] = cur_node_id

            # Check if this is the seed node
            if node.parent_node is None:
                parent_idx = None
            else:
                parent_idx = d_node_to_idx.get(node.parent_node)
                if parent_idx is None:
                    parent_idx = len(d_node_to_idx)
                    d_node_to_idx[node.parent_node] = parent_idx

            # If the taxon is present, save it
            if node.taxon is not None:
                if node.taxon.label in taxa:
                    raise Exception("Duplicate taxon found!")
                taxa.add((node.taxon.label, cur_node_id))

            # If this is not the seed node, store the edge
            if parent_idx is not None:
                edges.add((parent_idx, cur_node_id, node.edge_length))

        taxa, edges = tuple(taxa), tuple(edges)
        return taxa, edges

    def subset_to_taxa(self, taxa: Collection[str]):
        return self.tree.extract_tree_with_taxa_labels(taxa)

    @staticmethod
    def to_file(out_path: Path, taxa: Tuple[Tuple[str, int]], edges: Tuple[Tuple[int, int, float]]):
        """
        Write the nodes and edges to a file for later import in Rust.
        """
        print(f'Saving distance matrix to: {out_path}')
        with out_path.open('w') as f:
            f.write('#TAXA\n')
            for taxon, taxon_id in taxa:
                f.write(f'{taxon}\t{taxon_id}\n')
            f.write('#EDGES\n')
            for parent, child, length in edges:
                f.write(f'{parent}\t{child}\t{length}\n')
        return
