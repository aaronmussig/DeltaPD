from pathlib import Path

import typer
from typing_extensions import Annotated

from deltapd.method import run_deltapd
from deltapd.model.params import CorrelationFn, ErrorFn
from deltapd.util.cli import validate_cli_arguments, parse_only_taxa
from deltapd.util.logging import log
from deltapd.util.sysinfo import N_CPUS

app = typer.Typer()


@app.command()
def run(
        reference_tree: Annotated[Path, typer.Argument(
            help="Path to the reference tree (Newick format)."
        )],
        query_tree: Annotated[Path, typer.Argument(
            help="Path to the query tree (Newick format)."
        )],
        metadata: Annotated[Path, typer.Argument(
            help="Tab-delimited metadata file."
        )],
        output_dir: Annotated[Path, typer.Argument(
            help="Output directory for the results."
        )],
        cpus: Annotated[int, typer.Option(
            help="Number of CPUs to use.", min=1, max=N_CPUS
        )] = N_CPUS,
        sample_size: Annotated[float, typer.Option(
            help="Fraction of the data to sample for bootstrapping if <1, otherwise the absolute number of samples.",
            min=0.0
        )] = 0.1,
        replicates: Annotated[int, typer.Option(
            help="Number of bootstrap replicates to run.", min=1
        )] = 10,
        correlation_fn: Annotated[CorrelationFn, typer.Option(
            help="Correlation function to use.", case_sensitive=False
        )] = CorrelationFn.R2,
        error_fn: Annotated[ErrorFn, typer.Option(
            help="Error function to use.", case_sensitive=False
        )] = ErrorFn.RMSE,
        only_taxa: Annotated[str, typer.Option(
            help="Only analyse these taxa in the query tree (csv)."
        )] = None,
        debug: Annotated[bool, typer.Option(
            help="Verbose output, intermediate files are saved."
        )] = False,
):
    # Validation
    validate_cli_arguments(reference_tree, query_tree, metadata, output_dir)

    if 1 < cpus < N_CPUS:
        log(f'TODO: Implement a varying number of CPUs not min or max. Will use all CPUs for now.')

    only_taxa = parse_only_taxa(only_taxa)

    # Run the method
    run_deltapd(reference_tree, query_tree, metadata, output_dir, cpus, sample_size, replicates, correlation_fn,
                error_fn, only_taxa, debug)

    print('Done.')
    return


@app.command()
def test():
    print('test.')


if __name__ == '__main__':
    app()
