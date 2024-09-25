from pathlib import Path

import sys


def validate_cli_arguments(reference_tree: Path, query_tree: Path, metadata: Path, output_dir: Path):
    errors = list()

    # Check that the reference tree exists
    if not reference_tree.exists():
        errors.append(f'Reference tree does not exist: {reference_tree}')

    # Check that the query tree exists
    if not query_tree.exists():
        errors.append(f'Query tree does not exist: {query_tree}')

    # Check that the metadata file exists
    if not metadata.exists():
        errors.append(f'Metadata file does not exist: {metadata}')

    # Exit early if the output directory exists
    if output_dir.exists():
        print('TODO: Output directory already existing will be an error in a future version.')
        # errors.append(f'Output directory already exists: {output_dir}')

    # Terminate the program if any errors were found
    if len(errors) > 0:
        maybe_plural = 's' if len(errors) > 1 else ''
        print(f'DeltaPD is unable to run due to the following error{maybe_plural}:')
        for error in errors:
            print(f'  - {error}')
        sys.exit(1)


def parse_only_taxa(only_taxa: str) -> set[str]:
    return set(only_taxa.split(',')) if only_taxa else set()
