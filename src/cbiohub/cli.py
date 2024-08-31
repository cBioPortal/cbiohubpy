import os
from pathlib import Path
import shutil
import click
import importlib.metadata
from tqdm import tqdm
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
from .study import Study  # Assuming the Study class is in a file named study.py
from tabulate import tabulate

settings.PROCESSED_PATH = os.path.expanduser(settings.PROCESSED_PATH)


@click.group()
def cli():
    pass


@cli.command()
@click.argument("folder_name", type=click.Path(exists=True))
def ingest(folder_name):
    """Ingest studies from the given folder and create Parquet files."""
    folder_path = Path(folder_name)

    if folder_path.is_dir():
        # Check if the folder contains multiple studies or a single study
        study_paths = [
            p for p in folder_path.iterdir() if p.is_dir() and Study.is_study(p)
        ]
        if not study_paths:
            # Single study
            if Study.is_study(folder_path):
                study_paths = [folder_path]
            else:
                click.echo(f"No valid studies found in {folder_name}.")
                return

        processed_count = 0
        already_processed_count = 0
        skipped_due_to_errors_count = 0
        skipped_due_to_missing_files_count = 0

        skipped_due_to_errors_studies = []
        skipped_due_to_missing_files_studies = []

        with tqdm(
            total=len(study_paths), desc="Processing studies", unit="study"
        ) as pbar:
            for study_path in study_paths:
                study = Study(study_path)
                pbar.set_description(f"Processing {study.name}")
                if study.is_processed():
                    already_processed_count += 1
                    pbar.update(1)
                    continue

                if study.check_integrity():
                    if study.create_parquets():
                        processed_count += 1
                    else:
                        click.echo(
                            f"⚠️ Skipped study {study.name} due to errors during Parquet creation."
                        )
                        skipped_due_to_errors_count += 1
                        skipped_due_to_errors_studies.append(study.name)
                else:
                    click.echo(
                        f"⚠️ Skipped study {study.name} due to missing required files."
                    )
                    skipped_due_to_missing_files_count += 1
                    skipped_due_to_missing_files_studies.append(study.name)
                pbar.update(1)

        click.echo(
            click.style(
                f"✅ Finished processing {processed_count} studies.", fg="green"
            )
        )
        click.echo(
            click.style(f"ℹ️ {already_processed_count} studies were already processed.")
        )
        click.echo(
            click.style(
                f"⚠️ {skipped_due_to_errors_count} studies were skipped due to errors during Parquet creation.",
                fg="red",
            )
        )
        if skipped_due_to_errors_studies:
            click.echo(
                click.style(
                    f"Skipped due to errors: {', '.join(skipped_due_to_errors_studies)}",
                    fg="red",
                )
            )
        click.echo(
            click.style(
                f"⚠️ {skipped_due_to_missing_files_count} studies were skipped due to missing required files.",
                fg="yellow",
            )
        )
        if skipped_due_to_missing_files_studies:
            click.echo(
                click.style(
                    f"Skipped due to missing files: {', '.join(skipped_due_to_missing_files_studies)}",
                    fg="yellow",
                )
            )
    else:
        click.echo(f"Error: {folder_name} is not a directory.", fg="red")


@cli.command()
def clean():
    """Remove everything in the processed path folder."""
    processed_path = Path(settings.PROCESSED_PATH)

    if processed_path.is_dir():
        for item in processed_path.iterdir():
            if item.is_dir():
                shutil.rmtree(item)
            else:
                item.unlink()
        click.echo(f"Cleaned all files and directories in {processed_path}.")
    else:
        click.echo(f"Error: {processed_path} is not a directory.")


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


@cli.command()
@click.option(
    "--output-dir",
    type=click.Path(),
    default=None,
    help="Optional output directory for combined files.",
)
def combine(output_dir):
    """Combine all processed studies into a single combined processed study."""
    processed_studies_path = Path(settings.PROCESSED_PATH) / "studies"
    combined_path = (
        Path(output_dir) if output_dir else Path(settings.PROCESSED_PATH) / "combined"
    )

    combined_path.mkdir(parents=True, exist_ok=True)

    mutation_tables = []
    clinical_patient_tables = []
    clinical_sample_tables = []

    study_paths = [p for p in processed_studies_path.iterdir() if p.is_dir()]

    with tqdm(total=len(study_paths), desc="Loading studies", unit="study") as pbar:
        for study_path in processed_studies_path.iterdir():
            if study_path.is_dir():
                study = Study(study_path)

                if not study.is_processed():
                    click.echo(f"⚠️ Skipping {study_path} (not successfully processed)")
                    pbar.update(1)
                    continue

                mutation_file = study.processed_path / "data_mutations.parquet"
                clinical_patient_file = (
                    study.processed_path / "data_clinical_patient.parquet"
                )
                clinical_sample_file = (
                    study.processed_path / "data_clinical_sample.parquet"
                )

                if mutation_file.exists():
                    table = pq.read_table(mutation_file)
                    # Select only specific columns and adjust their types
                    columns_to_include = MUTATION_COLUMNS

                    # Filter out columns that do not exist in the table schema
                    existing_columns = {
                        col: dtype
                        for col, dtype in columns_to_include.items()
                        if col in table.schema.names
                    }

                    table = table.select(list(existing_columns.keys()))
                    table = table.cast(pa.schema(existing_columns))
                    mutation_tables.append(table)
                if clinical_patient_file.exists():
                    clinical_patient_tables.append(pq.read_table(clinical_patient_file))
                if clinical_sample_file.exists():
                    clinical_sample_tables.append(pq.read_table(clinical_sample_file))

                pbar.update(1)

    if mutation_tables:
        start_time = time.time()
        combined_mutations = pa.concat_tables(mutation_tables, promote=True)
        concat_time = time.time() - start_time

        click.echo(
            click.style(f"⏱️ Concatenation time: {concat_time} seconds", fg="green")
        )

        start_time = time.time()
        pq.write_table(combined_mutations, combined_path / "combined_mutations.parquet")
        write_time = time.time() - start_time

        click.echo(
            click.style(
                f"✅ Combined mutations saved to {combined_path / 'combined_mutations.parquet'}",
                fg="green",
            )
        )
        click.echo(click.style(f"⏱️ Write time: {write_time} seconds", fg="green"))

    if clinical_patient_tables:
        start_time = time.time()
        combined_clinical_patient = pa.concat_tables(
            clinical_patient_tables, promote=True
        )
        concat_time = time.time() - start_time

        click.echo(
            click.style(f"⏱️ Concatenation time: {concat_time} seconds", fg="green")
        )

        start_time = time.time()
        pq.write_table(
            combined_clinical_patient,
            combined_path / "combined_clinical_patient.parquet",
        )
        write_time = time.time() - start_time

        click.echo(
            click.style(
                f"✅ Combined clinical patient data saved to {combined_path / 'combined_clinical_patient.parquet'}",
                fg="green",
            )
        )
        click.echo(click.style(f"⏱️ Write time: {write_time} seconds", fg="green"))

    if clinical_sample_tables:
        start_time = time.time()
        combined_clinical_sample = pa.concat_tables(
            clinical_sample_tables, promote=True
        )
        concat_time = time.time() - start_time

        click.echo(
            click.style(f"⏱️ Concatenation time: {concat_time} seconds", fg="green")
        )

        start_time = time.time()
        pq.write_table(
            combined_clinical_sample, combined_path / "combined_clinical_sample.parquet"
        )
        write_time = time.time() - start_time

        click.echo(
            click.style(
                f"✅ Combined clinical sample data saved to {combined_path / 'combined_clinical_sample.parquet'}",
                fg="green",
            )
        )
        click.echo(click.style(f"⏱️ Write time: {write_time} seconds", fg="green"))


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
@click.option(
    "--processed-dir",
    default=None,
    help="Directory containing the processed parquet files",
)
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
