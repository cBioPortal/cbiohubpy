from pathlib import Path

import pyarrow.parquet as pq
import pyarrow.dataset as ds
import duckdb
import pyarrow as pa
import pandas as pd
from dynaconf import settings

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


def get_combined_df(directory=None):
    """Get combined study data."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    mut = pd.read_parquet(
        directory / "combined_mutations.parquet", columns=MUTATION_COLUMNS
    )
    clinp = pd.read_parquet(directory / "combined_clinical_patient.parquet")
    clins = pd.read_parquet(directory / "combined_clinical_sample.parquet")

    return mut, clinp, clins


def find_samples_in_parquet(filter_expression, directory):
    dataset = ds.dataset(directory / "combined_mutations.parquet", format="parquet")
    table = dataset.to_table(
        filter=filter_expression, columns=["Tumor_Sample_Barcode", "study_id"]
    )
    if table.num_rows > 0:
        unique_identifiers = [
            f"{study_id}:{barcode}"
            for study_id, barcode in zip(
                table["study_id"], table["Tumor_Sample_Barcode"]
            )
        ]
        return True, unique_identifiers
    else:
        return False, []


def variant_exists(chrom, start, end, ref, alt, directory=None):
    """Check if a particular variant exists in the combined mutations parquet."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    filter_expression = (
        (ds.field("Chromosome") == chrom)
        &
        # TODO: treat everything as string for now
        (ds.field("Start_Position") == str(start))
        & (ds.field("End_Position") == str(end))
        & (ds.field("Reference_Allele") == ref)
        & (ds.field("Tumor_Seq_Allele2") == alt)
    )

    return find_samples_in_parquet(filter_expression, directory)


def variant_exists_by_protein_change(hugo_symbol, protein_change, directory=None):
    """Check if a particular variant exists based on Hugo symbol and protein change."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    protein_change = (
        protein_change if protein_change.startswith("p.") else f"p.{protein_change}"
    )

    filter_expression = (ds.field("Hugo_Symbol") == hugo_symbol) & (
        ds.field("HGVSp_Short") == protein_change
    )

    return find_samples_in_parquet(filter_expression, directory)


def find_variant(
    chrom=None,
    start=None,
    end=None,
    ref=None,
    alt=None,
    hugo_symbol=None,
    protein_change=None,
    directory=None,
):
    """Find a variant based on either genomic coordinates or Hugo symbol and protein change."""
    if chrom and start and end and ref and alt:
        return variant_exists(chrom, start, end, ref, alt, directory)
    elif hugo_symbol and protein_change:
        return variant_exists_by_protein_change(hugo_symbol, protein_change, directory)
    else:
        raise ValueError("Insufficient arguments provided to find a variant.")


def variant_frequency_per_cancer_type(
    chrom, start, end, ref, alt, clinical_attribute, directory=None
):
    """Check how frequently a particular variant occurs per cancer type."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    mutations_path = directory / "combined_mutations.parquet"
    clinical_path = directory / "combined_clinical_sample.parquet"

    query = f"""
    SELECT clinical.{clinical_attribute}, COUNT(*) as frequency
    FROM '{mutations_path}' AS mutations
    JOIN '{clinical_path}' AS clinical
    ON mutations.Tumor_Sample_Barcode = clinical.SAMPLE_ID
    WHERE mutations.Chromosome = '{chrom}'
    AND mutations.Start_Position = '{start}'
    AND mutations.End_Position = '{end}'
    AND mutations.Reference_Allele = '{ref}'
    AND mutations.Tumor_Seq_Allele2 = '{alt}'
    GROUP BY clinical.{clinical_attribute}
    ORDER BY frequency DESC
    """

    con = duckdb.connect()
    result = con.execute(query).fetchall()
    con.close()

    return result


def get_genomic_coordinates_by_gene_and_protein_change(
    gene, protein_change, directory=None
):
    """Convert a gene name and protein change to its corresponding genomic coordinates and count occurrences."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    mutations_path = directory / "combined_mutations.parquet"

    # Ensure the protein change starts with "p."
    protein_change = (
        protein_change if protein_change.startswith("p.") else f"p.{protein_change}"
    )

    query = f"""
    SELECT 
        Chromosome, 
        Start_Position, 
        End_Position, 
        Reference_Allele, 
        Tumor_Seq_Allele2, 
        COUNT(*) as frequency
    FROM '{mutations_path}'
    WHERE Hugo_Symbol = '{gene}'
    AND HGVSp_Short = '{protein_change}'
    GROUP BY Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2
    ORDER BY frequency DESC
    """

    con = duckdb.connect()
    result = con.execute(query).fetchall()
    con.close()

    return result
