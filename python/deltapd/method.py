from pathlib import Path

from deltapd import PyDistMatrix, PyDeltaPD, PyLinearModelType, PyParams
from deltapd.model.params import CorrelationFn, ErrorFn
from deltapd.tree import DeltaPdTree


def run_deltapd(
        ref_path: Path,
        qry_path: Path,
        meta_path: Path,
        output_dir: Path,
        cpus: int,
        sample_size: float,
        replicates: int,
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
        sample_size=sample_size,
        replicates=replicates,
        taxa=only_taxa if only_taxa else set(),
        model=PyLinearModelType.RepeatedMedian,
        error=error_fn.to_rs(),
        corr=correlation_fn.to_rs(),
        debug=debug,
        output_dir=str(output_dir.absolute())
    )

    # Run the analysis
    results = dpd.run(dpd_params)
    #
    # # Write the data to disk
    # path_out = output_dir / 'results.tsv'
    # with path_out.open('w') as f:
    #     f.write(f'taxon\t')
    #     for result in results:
    #         f.write(f'')
    #
    print('Done with the analysis. Writing the results to disk.')
    print(len(results))
    print('TODO!')
    return
