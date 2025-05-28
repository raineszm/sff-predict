import shutil
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

# Download latest version
path = kagglehub.dataset_download(
    "zacharymraines/open-library-works-dump-2025-01-08", path="ol_works.parquet"
)

# Move it to the data directory
shutil.move(path, snakemake.output[0])
