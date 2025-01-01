import os
from pathlib import Path
import shutil
from cbiohub.variant import GenomicVariant, ProteinVariant
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
    variant_frequency_per_clinical_attribute,
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
    func = click.option(
        "--sql",
        is_flag=True,
        default=False,
        help="Output SQL command rather than the result",
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
@common_options
def find(arg1, arg2, arg3, arg4, arg5, processed_dir, sql):
    """Find a variant in the combined mutations parquet and return details."""
    if arg1 and arg2 and arg3 and arg4:
        # assuming chrom/pos/start/end
        exists, unique_ids = find_variant(
            chrom=arg1, start=arg2, end=arg3, ref=arg4, alt=arg5, directory=processed_dir
        )
    elif arg1 and arg2:
        # assuming gene/protein_change
        exists, unique_ids = find_variant(hugo_symbol=arg1, protein_change=arg2, directory=processed_dir)
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

def detect_variant_type(args):
    """
    Detects the variant type based on the number and structure of arguments.

    Returns:
        - GenomicVariant if args match chrom/start/end/ref/alt
        - ProteinVariant if args match gene/protein_change
    """
    if len(args) == 5:
        return GenomicVariant(*args)
    elif len(args) == 2:
        return ProteinVariant(*args)
    else:
        raise click.UsageError("Could not detect variant type. Ensure arguments match either chrom/start/end/ref/alt or gene/protein_change format.")

@cli.command(help="Check how frequently a particular variant occurs per cancer type.")
@click.argument('args', nargs=-1, required=True)
@click.option(
    "--clinical-attribute",
    default="CANCER_TYPE",
    help="Clinical attribute to group by (default: CANCER_TYPE)",
)
@click.option(
    "--group-by-study-id",
    is_flag=True,
    default=False,
    help="Group the results by cancer study id",
)
@click.option(
    "--count-samples",
    is_flag=True,
    default=False,
    help="Count samples instead of patients",
)
@common_options
def variant_frequency(args, clinical_attribute, processed_dir, sql, group_by_study_id, count_samples):
    """Check how frequently a particular variant occurs per cancer type (or
    other clinical sample attributes)."""
    variant = detect_variant_type(args)

    result = variant_frequency_per_clinical_attribute(
        variant, clinical_attribute, directory=processed_dir, sql=sql, group_by_study_id=group_by_study_id,
        count_samples=count_samples
    )
    if sql:
        return print(result)

    if result:
        click.echo(
            click.style(f"✅ Variant frequency per {clinical_attribute}:", fg="green")
        )
        headers = [clinical_attribute, "altered", "total", "freq"]
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
def convert(gene, protein_change, processed_dir, sql):
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
            headers = ["Chromosome", "Start", "End", "Ref", "Alt", "Count"]
            table = tabulate(results, headers, tablefmt="plain")
            click.echo(table)
        else:
            click.echo(click.style("❌ No data found.", fg="red"))
    except ValueError as e:
        click.echo(click.style(str(e), fg="red"))


if __name__ == "__main__":
    cli()
