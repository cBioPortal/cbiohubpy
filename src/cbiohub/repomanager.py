from pathlib import Path
from dynaconf import settings
from study import Study
import glob


class RepoManager:
    def __init__(self, base_path):
        self.base_path = Path(base_path).expanduser()
        self.studies = self._load_studies()

    def _load_studies(self):
        """Load all studies from the base directory based on the presence of meta_study.txt."""
        studies = []
        for study_path in glob.glob(str(self.base_path / "**"), recursive=True):
            if (
                Path(study_path).is_dir()
                and (Path(study_path) / "meta_study.txt").exists()
            ):
                studies.append(Study(Path(study_path)))
        return studies

    def list_studies(self):
        """List all study directories."""
        return [study.study_path.name for study in self.studies]

    def get_study(self, study_name):
        """Get a specific study by name."""
        for study in self.studies:
            if study.study_path.name == study_name:
                return study
        raise FileNotFoundError(f"Study {study_name} not found.")


# Example usage
if __name__ == "__main__":
    manager = RepoManager(settings.DATAHUB_PATH)

    # List all studies
    studies = manager.list_studies()
    print("Studies:", studies)

    # Get a specific study
    study_name = "msk_impact_2017"
    try:
        study = manager.get_study(study_name)
        print(f"Sample files in {study_name}:", study.list_files("sample"))
        print(f"Patient files in {study_name}:", study.list_files("patient"))
        print(f"Mutation files in {study_name}:", study.list_files("mutation"))
    except FileNotFoundError as e:
        print(e)

    # Get a specific file in a study
    file_name = "data_clinical_sample.txt"
    try:
        file_path = study.get_file(file_name)
        print(f"Path to {file_name}:", file_path)
    except FileNotFoundError as e:
        print(e)
