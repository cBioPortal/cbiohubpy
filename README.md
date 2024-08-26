# cbiohub

**WARNING ⚠️: This package is still under construction.**


`cbiohub` is a Python package that provides convenience functions for analyzing
data files from [cBioPortal](https://cbioportal.org). Although several Python
API clients exist, they work on slices of the cBioPortal data retrieved via the
REST rather than that they enable easy analysis of all the data files in bulk.
This package aims to provide a more user-friendly interface for accessing data
from cBioPortal like those stored in the public
[datahub](https://github.com/cBioPortal/datahub).

## Usage

### Analyze Local files

#### Step 1: Obtain data files

You can e.g. download the cBioPortal datahub files:

```sh
git clone git@github.com:cbioportal/datahub
```

### Step 2: Ingest and Combine

Now ingest them i.e. convert them into parquet files on your local machine:

```sh

cbiohub ingest ~/git/datahub/public/
```

All the data by default gets stored in `~/cbiohub/`. Now you can combine all the study data together into a single combined study:

```sh
cbiohub combine
```

### Step 3: Analyze

Now you can use the `cbiohub` package to analyze the data. For example, you can load the combined study data into a pandas DataFrame:

```python
import cbiohub

df = cbiohub.get_combined_df()
```

Or you can use the cbiohub cli to do quick analyses:

```sh
> cbiohub find BRAF V600E
✅ Variant found in 3595 samples across 117 studies:
kirp_tcga:TCGA-AL-3467-01
kirp_tcga:TCGA-UZ-A9PP-01
...
```

or search for the same BRAF V600E variant but with a specific genomic change
(A>T):

```sh
> cbiohub find 7 140453136 140453136 A T
✅ Variant found in 3571 samples across 117 studies:
kirp_tcga:TCGA-AL-3467-01
kirp_tcga:TCGA-UZ-A9PP-01
...
```

### Clean

Remove all local parquet files.

```sh
cbiohub clean
```

## Development

To set up the development environment, install the development dependencies:

```sh
poetry install
```

You can run the cli using e.g.:

```sh
poetry ingest ~/git/datahub/public/
```

and

```sh
poetry run cbiohub find BRAF V600E
```

You can also use IPython for interactive exploration:

```sh
poetry run ipython
```

## TODO

- [ ] Add github action datahub that usies cbiohub to push combined parquet data to hugging face (https://huggingface.co/datasets/cBioPortal/datahub)