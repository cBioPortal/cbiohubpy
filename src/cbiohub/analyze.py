from pathlib import Path

import pandas as pd
from dynaconf import settings

MUTATION_COLUMNS = [
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele1",
    "Tumor_Seq_Allele2",
    "t_ref_count",
    "t_alt_count",
    "n_ref_count",
    "n_alt_count",
    "Hugo_Symbol",
    "Tumor_Sample_Barcode",
    "study_id",
]


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
