import os
from pathlib import Path
import pandas as pd
from dynaconf import (
    settings,
)  # Assuming settings is a module with PROCESSED_PATH defined

MUTATION_DATA_FILES = ["data_mutations.txt", "data_mutations_extended.txt", "data_nonsignedout_mutations.txt"]

class Study:
    def __init__(self, study_path: Path):
        self.study_path = study_path
        self.name = study_path.name
        self.processed_path = Path(settings.PROCESSED_PATH) / "studies" / self.name
        self.sample_data_file = "data_clinical_sample.txt"
        self.patient_data_file = "data_clinical_patient.txt"
        # mutation data can have multiple names
        self.mutation_data_files = [
            mut_file for mut_file in MUTATION_DATA_FILES if (self.study_path / mut_file).exists()
        ]
        # all mutation data files are combined into this file
        self.mutation_parquet_file = "data_mutations.parquet"
        self.sample_df = None
        self.patient_df = None
        self.mutation_df = None

    @classmethod
    def is_study(cls, path):
        """Check if the given path is a study by looking for 'meta_study.txt'."""
        return (Path(path) / "meta_study.txt").exists()

    def check_integrity(self):
        """Check if the study contains all required files."""
        required_files = [
            self.sample_data_file,
            self.patient_data_file,
            "meta_study.txt",
        ]
        # Check if at least one mutation data file exists
        if len(self.mutation_data_files) == 0:
            missing_files.append("at least one mutation data file")

        missing_files = [
            file_name
            for file_name in required_files
            if not (self.study_path / file_name).exists()
        ]
        if missing_files:
            print(
                f"Study {self.name} is missing the following files: {', '.join(missing_files)}"
            )
            return False
        return True

    def list_files(self, file_type: str):
        """List the file in the study directory matching the file type."""
        if not self.study_path.exists():
            raise FileNotFoundError(f"Study directory {self.study_path} not found.")

        if file_type == "sample":
            file_name = self.sample_data_file
        elif file_type == "patient":
            file_name = self.patient_data_file
        elif file_type == "mutation":
            file_name = self.mutation_data_file
        else:
            raise ValueError(f"Unknown file type: {file_type}")

        file_path = self.study_path / file_name
        if not file_path.exists():
            raise FileNotFoundError(
                f"File {file_name} not found in study {self.study_path}."
            )

        return [file_path]

    def create_parquet(self, file_type: str):
        """Create a Parquet file for the specified file type in the PROCESSED_PATH folder."""
        if file_type == "sample":
            file_names = [self.sample_data_file]
        elif file_type == "patient":
            file_names = [self.patient_data_file]
        elif file_type == "mutation":
            file_names = self.mutation_data_files
        else:
            raise ValueError(f"Unknown file type: {file_type}")

        for file_name in file_names:
            file_path = self.study_path / file_name
            if not file_path.exists():
                raise FileNotFoundError(
                    f"File {file_name} not found in study {self.study_path}."
                )

            try:
                # Read the data file into a DataFrame
                df = pd.read_csv(
                    file_path, sep="\t", comment="#", low_memory=False, dtype=str
                )
                # add study_id as a column
                df["study_id"] = self.name

                # Define the output directory and file path
                output_dir = self.processed_path
                output_dir.mkdir(parents=True, exist_ok=True)
                output_file = output_dir / file_name.replace(".txt", ".parquet")

                # Write the DataFrame to a Parquet file
                df.to_parquet(output_file)
            except pd.errors.ParserError as e:
                print(f"Parse error in study {self.name} for file {file_name}: {e}")
                return False

        return True

    def get_parquet(self, file_type: str):
        """Get the DataFrame for the specified file type, generating the Parquet file if necessary."""
        if file_type == "sample":
            file_names = [self.sample_data_file]
            df_attr = "sample_df"
        elif file_type == "patient":
            file_names = [self.patient_data_file]
            df_attr = "patient_df"
        elif file_type == "mutation":
            file_names = self.mutation_data_files
            df_attr = "mutation_df"
        else:
            raise ValueError(f"Unknown file type: {file_type}")

        for file_name in file_names:
            file_path = self.study_path / file_name
            output_dir = Path(settings.PROCESSED_PATH) / self.name
            output_file = output_dir / file_name.replace(".txt", ".parquet")

            # Check if the Parquet file needs to be created or updated
            if (
                not output_file.exists()
                or file_path.stat().st_mtime > output_file.stat().st_mtime
            ):
                self.create_parquet(file_type)

            # Load the DataFrame from the Parquet file if not already loaded
            if getattr(self, df_attr) is None:
                setattr(self, df_attr, pd.read_parquet(output_file))

        return getattr(self, df_attr)

    def get_sample_df(self):
        """Get the DataFrame for the sample data."""
        return self.get_parquet("sample")

    def get_patient_df(self):
        """Get the DataFrame for the patient data."""
        return self.get_parquet("patient")

    def get_mutation_df(self):
        """Get the DataFrame for the mutation data."""
        return self.get_parquet("mutation")

    def create_parquets(self):
        """Create Parquet files for sample, patient, and mutation data in the PROCESSED_PATH folder."""
        success = True
        success &= self.create_parquet("sample")
        success &= self.create_parquet("patient")
        success &= self.create_parquet("mutation")

        if success:
            # Create a success indicator file
            success_file = self.processed_path / "ingestion_success.txt"
            success_file.touch()
        return success

    def is_processed(self):
        """Check if the study has already been processed."""
        success_file = self.processed_path / "ingestion_success.txt"
        if not success_file.exists():
            return False

        # Check if any source file is newer than the success indicator file
        source_files = [
            self.study_path / self.sample_data_file,
            self.study_path / self.patient_data_file,
            *[self.study_path / mutation_data_file for mutation_data_file in self.mutation_data_files],
            self.study_path / "meta_study.txt",
        ]

        success_file_mtime = success_file.stat().st_mtime
        for source_file in source_files:
            if (
                source_file.exists()
                and source_file.stat().st_mtime > success_file_mtime
            ):
                return False

        return True

    def create_sample_parquet(self):
        """Create a Parquet file for the sample data."""
        self.create_parquet("sample")

    def create_patient_parquet(self):
        """Create a Parquet file for the patient data."""
        self.create_parquet("patient")

    def create_mutation_parquet(self):
        """Create a Parquet file for the mutation data."""
        self.create_parquet("mutation")
