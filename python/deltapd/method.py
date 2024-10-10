from pathlib import Path

from deltapd import PyDistMatrix, PyDeltaPD, PyLinearModelType, PyParams
from deltapd.model.params import CorrelationFn, ErrorFn, Direction
from deltapd.tree import DeltaPdTree
import matplotlib.pyplot as plt
import numpy as np


def run_deltapd(
        ref_path: Path,
        qry_path: Path,
        meta_path: Path,
        output_dir: Path,
        cpus: int,
        knn: int,
        direction: Direction,
        correlation_fn: CorrelationFn,
        error_fn: ErrorFn,
        only_taxa: set[str],
        debug: bool
):
    # Create the output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the reference tree
    print('Processing the reference tree.')
    ref_tree = DeltaPdTree(ref_path)
    ref_taxa, ref_edges = ref_tree.get_nodes_and_edges_for_deltapd()
    ref_dm = PyDistMatrix(taxa=ref_taxa, edges=ref_edges)

    # Load the query tree
    print('Processing the query tree.')
    qry_tree = DeltaPdTree(qry_path)
    qry_taxa, qry_edges = qry_tree.get_nodes_and_edges_for_deltapd()
    qry_dm = PyDistMatrix(taxa=qry_taxa, edges=qry_edges)

    # Save the distance matrices to disk if we're debugging
    if debug:
        ref_tree.to_file(out_path=output_dir / f'{ref_path.name}.dm', taxa=ref_taxa, edges=ref_edges)
        qry_tree.to_file(out_path=output_dir / f'{qry_path.name}.dm', taxa=qry_taxa, edges=qry_edges)

    # Create the DeltaPD object
    dpd = PyDeltaPD(qry_dm, ref_dm, str(meta_path.absolute()), '\t')

    # Create the parameter object
    dpd_params = PyParams(
        cpus=cpus,
        sample_size=1,
        replicates=1,
        knn=knn,
        direction=direction.to_rs(),
        taxa=only_taxa if only_taxa else set(),
        model=PyLinearModelType.RepeatedMedian,
        error=error_fn.to_rs(),
        corr=correlation_fn.to_rs(),
        debug=debug,
        output_dir=str(output_dir.absolute())
    )

    # Run the analysis
    results = dpd.run(dpd_params)

    print('Getting the terminal branch lengths.')
    d_taxon_to_term_branch = qry_tree.get_terminal_branch_lengths(only_taxa)

    result_path = output_dir / 'results.tsv'
    with result_path.open('w') as f:
        header = ('taxon', 'mean_error', 'median_error', 'std_error', 'terminal_branch_length')
        f.write('\t'.join(header) + '\n')
        for result in sorted(results, key=lambda x: -x.error_median):
            error_mean = round(result.error_mean * 100, 4)
            error_std = round(result.error_std * 100, 4)
            error_median = round(result.error_median * 100, 4)
            cur_term_branch = d_taxon_to_term_branch[result.taxon]

            col_values = (result.taxon, error_mean, error_median, error_std, cur_term_branch)
            f.write('\t'.join(map(str, col_values)) + '\n')

    print(f'Results are here: {result_path}')


    # If we're debugging then also plot the x/y data


    return
