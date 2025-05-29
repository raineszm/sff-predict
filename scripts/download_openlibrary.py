import os
import kagglehub
from kagglehub.exceptions import UnauthenticatedError
from snakemake.script import snakemake


# Make sure we're logged in
try:
    print(kagglehub.whoami())
except UnauthenticatedError:
    print("Not logged in, please login.")
    print("You can skip this step by following the instructions at:")
    print("https://github.com/Kaggle/kagglehub")
    kagglehub.login()
    print(kagglehub.whoami())

# Set the cache directory to be in the raw data directory
os.environ["KAGGLEHUB_CACHE"] = os.path.join(
    os.path.dirname(snakemake.output[0]), "kaggle_cache"
)

# Download latest version
path = kagglehub.dataset_download(
    "zacharymraines/open-library-works-dump-2025-01-08", path="ol_works.parquet"
)

# Link it to the data directory
os.symlink(os.path.abspath(path), os.path.abspath(snakemake.output[0]))
