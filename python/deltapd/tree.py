from pathlib import Path
from typing import Collection

import dendropy


class Tree:

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
            child_idx = d_node_to_idx.get(node)
            if child_idx is None:
                child_idx = len(d_node_to_idx)
                d_node_to_idx[node] = child_idx
        
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
                    raise Exception("??")
                taxa.add((node.taxon.label, child_idx))
                
            # If this is not the seed node, store the edge
            if parent_idx is not None:
                edges.add((parent_idx, child_idx, node.edge_length))
                
        return tuple(taxa), tuple(edges)

    def subset_to_taxa(self, taxa: Collection[str]):
        return self.tree.extract_tree_with_taxa_labels(taxa)
