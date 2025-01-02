import time
import click
from tqdm import tqdm
from pathlib import Path
import pyarrow.parquet as pq
import pyarrow as pa
from dynaconf import (
    settings,
)  # Assuming settings is a module with PROCESSED_PATH defined

from cbiohub.study import Study, MUTATION_DATA_FILES

MUTATION_COLUMNS = {
    "Chromosome": pa.string(),
    "Start_Position": pa.string(),
    "End_Position": pa.string(),
    "Reference_Allele": pa.string(),
    "Tumor_Seq_Allele1": pa.string(),
    "Tumor_Seq_Allele2": pa.string(),
    "t_ref_count": pa.string(),
    "t_alt_count": pa.string(),
    "n_ref_count": pa.string(),
    "n_alt_count": pa.string(),
    "Hugo_Symbol": pa.string(),
    "HGVSp_Short": pa.string(),
    "Tumor_Sample_Barcode": pa.string(),
    "study_id": pa.string(),
}


@click.group()
def data():
    """Commands for data ingestion, and conversion to parquet."""
    pass


@data.command()
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


@data.command()
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

                clinical_patient_file = (
                    study.processed_path / "data_clinical_patient.parquet"
                )
                clinical_sample_file = (
                    study.processed_path / "data_clinical_sample.parquet"
                )

                for mutation_file in MUTATION_DATA_FILES:
                    mutation_file = study.processed_path / mutation_file.replace("txt","parquet")
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


@data.command()
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
