import os
from pathlib import Path
import shutil
import click
import importlib.metadata
import pandas as pd
import pyarrow as pa
import time
from dynaconf import (
    settings,
)
from .analyze import (
    find_variant,
    get_genomic_coordinates_by_gene_and_protein_change,
    variant_frequency_per_cancer_type,
    MUTATION_COLUMNS,
)
from .data_commands import data  # Import the data subcommand group
from .study import Study  # Assuming the Study class is in a file named study.py
from tabulate import tabulate

settings.PROCESSED_PATH = os.path.expanduser(settings.PROCESSED_PATH)


def common_options(func):
    """Decorator to add common options to Click commands."""
    func = click.option(
        "--processed-dir",
        default=None,
        help="Directory containing the processed parquet files",
    )(func)
    return func


@click.group()
def cli():
    pass


cli.add_command(data)


@cli.command()
def config():
    """Display the current configuration settings."""
    click.echo("Current Configuration:")
    click.echo(f"Processed Path: {settings.PROCESSED_PATH}")


@cli.command()
def version():
    """Display the current version of the tool."""
    try:
        version = importlib.metadata.version("cbiohub")
        click.echo(f"{version}")
    except importlib.metadata.PackageNotFoundError:
        click.echo("Package metadata not found. Ensure the project is installed.")


class CustomCommand(click.Command):
    def format_usage(self, ctx, formatter):
        formatter.write_usage(
            ctx.command_path, "[CHROM START END REF ALT] | [GENE PROTEIN_CHANGE]"
        )


@cli.command(
    cls=CustomCommand,
    help="Find a variant in the combined mutations parquet and return details. "
    "You must provide either (chrom, start, end, ref, alt) or (gene, protein_change).",
)
@click.argument("arg1", required=False)
@click.argument("arg2", required=False)
@click.argument("arg3", required=False)
@click.argument("arg4", required=False)
@click.argument("arg5", required=False)
def find(arg1, arg2, arg3, arg4, arg5):
    """Find a variant in the combined mutations parquet and return details."""
    if arg1 and arg2 and arg3 and arg4:
        # assuming chrom/pos/start/end
        exists, unique_ids = find_variant(
            chrom=arg1, start=arg2, end=arg3, ref=arg4, alt=arg5
        )
    elif arg1 and arg2:
        # assuming gene/protein_change
        exists, unique_ids = find_variant(hugo_symbol=arg1, protein_change=arg2)
    else:
        click.echo(click.style("❌ Invalid arguments.", fg="red"))
        return

    if exists:
        studies = set([id.split(":")[0] for id in unique_ids])
        click.echo(
            click.style(
                f"✅ Variant found in {len(unique_ids)} samples across {len(studies)} studies:",
                fg="green",
            )
        )
        for unique_id in unique_ids:
            click.echo(click.style(unique_id, fg="green"))
    else:
        click.echo(click.style("❌ Variant not found.", fg="red"))


@cli.command(help="Check how frequently a particular variant occurs per cancer type.")
@click.argument("chrom")
@click.argument("start", type=int)
@click.argument("end", type=int)
@click.argument("ref")
@click.argument("alt")
@click.option(
    "--clinical-attribute",
    default="CANCER_TYPE",
    help="Clinical attribute to group by (default: CANCER_TYPE)",
)
@common_options
def variant_frequency(chrom, start, end, ref, alt, clinical_attribute, processed_dir):
    """Check how frequently a particular variant occurs per cancer type (or
    other clinical sample attributes)."""
    result = variant_frequency_per_cancer_type(
        chrom, start, end, ref, alt, clinical_attribute, directory=processed_dir
    )
    if result:
        click.echo(
            click.style(f"✅ Variant frequency per {clinical_attribute}:", fg="green")
        )
        headers = ["Cancer Type", "Count"]
        table = tabulate(result, headers, tablefmt="plain")
        click.echo(table)
    else:
        click.echo(click.style("❌ No data found.", fg="red"))


@cli.command(
    help="Convert a gene and protein change to genomic coordinates and count occurrences."
)
@click.argument("gene")
@click.argument("protein_change")
@click.option(
    "--processed-dir",
    default=None,
    help="Directory containing the processed parquet files",
)
@common_options
def convert(gene, protein_change, processed_dir):
    """Convert a gene and protein change to its corresponding genomic coordinates and count occurrences."""
    try:
        results = get_genomic_coordinates_by_gene_and_protein_change(
            gene, protein_change, directory=processed_dir
        )
        if results:
            click.echo(
                click.style(
                    f"✅ Genomic coordinates for {gene} {protein_change}:", fg="green"
                )
            )
            headers = ["Chromosome", "Start", "End", "Ref", "Alt", "Frequency"]
            table = tabulate(results, headers, tablefmt="plain")
            click.echo(table)
        else:
            click.echo(click.style("❌ No data found.", fg="red"))
    except ValueError as e:
        click.echo(click.style(str(e), fg="red"))


if __name__ == "__main__":
    cli()
