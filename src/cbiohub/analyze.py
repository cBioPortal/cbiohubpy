from pathlib import Path
from typing import Union
import pandas as pd
import pyarrow.dataset as ds
import duckdb

from cbiohub.variant import GenomicVariant, ProteinVariant
from cbiohub.data_commands import MUTATION_COLUMNS


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


def variant_frequency_per_clinical_attribute(
    variant: Union[GenomicVariant, ProteinVariant], clinical_attribute,
    directory=None, group_by_study_id=False, count_samples=False, sql=False
):
    """Check how frequently a particular variant occurs per cancer type."""
    if directory is None:
        directory = Path(settings.PROCESSED_PATH) / "combined"
    else:
        directory = Path(directory)

    mutations_path = directory / "combined_mutations.parquet"
    clinical_sample_path = directory / "combined_clinical_sample.parquet"
    count_attribute = "SAMPLE_ID" if count_samples else "PATIENT_ID"

    if isinstance(variant, GenomicVariant):
        where_clause = f"""
            mutations.Chromosome = '{variant.chrom}'
            AND mutations.Start_Position = '{variant.start}'
            AND mutations.End_Position = '{variant.end}'
            AND mutations.Reference_Allele = '{variant.ref}'
            AND mutations.Tumor_Seq_Allele2 = '{variant.alt}'
        """
    elif isinstance(variant, ProteinVariant):
        where_clause = f"""
            mutations.Hugo_Symbol = '{variant.gene}'
            AND mutations.HGVSp_Short = '{variant.protein_change}'
        """
    else:
        raise TypeError("Unsupported variant type.")

    select_clause = f"""
        clinical.{clinical_attribute},
        COUNT(DISTINCT clinical.{count_attribute}) AS altered,
        TotalPerAttribute.total,
        ROUND((COUNT(DISTINCT clinical.{count_attribute}) * 100.0 / TotalPerAttribute.total), 1) AS freq
    """
    group_by_clause = f"clinical.{clinical_attribute}, TotalPerAttribute.total"
    total_group_by = f"clinical.{clinical_attribute}"
    total_join = f"clinical.{clinical_attribute} = TotalPerAttribute.{clinical_attribute}"
    total_select = f"clinical.{clinical_attribute}, COUNT(DISTINCT {count_attribute}) AS total"

    if group_by_study_id:
        total_select = f"clinical.STUDY_ID, {total_select}"
        total_group_by = f"clinical.STUDY_ID, {total_group_by}"
        total_join = f"clinical.STUDY_ID = TotalPerAttribute.STUDY_ID AND {total_join}"
        select_clause = f"clinical.STUDY_ID, {select_clause}"
        group_by_clause = f"clinical.STUDY_ID, {group_by_clause}"

    query = f"""
    WITH TotalPerAttribute AS (
        SELECT {total_select}
        FROM '{clinical_sample_path}' AS clinical
        GROUP BY {total_group_by}
    )
    SELECT
        {select_clause}
    FROM '{mutations_path}' AS mutations
    JOIN '{clinical_sample_path}' AS clinical
        ON mutations.Tumor_Sample_Barcode = clinical.SAMPLE_ID
    JOIN TotalPerAttribute
        ON {total_join}
    WHERE
        {where_clause}
    GROUP BY
        {group_by_clause}
    ORDER BY
        freq DESC;
    """

    if not sql:
        con = duckdb.connect()
        result = con.execute(query).fetchall()
        con.close()
        return result
    else:
        return query


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
