from tqdm.auto import tqdm
import httpx
from snakemake.script import snakemake

DATASET_URL = "https://www.kaggle.com/api/v1/datasets/download/zacharymraines/open-library-works-dump-2025-01-08/ol_works.parquet"

with httpx.stream("GET", DATASET_URL, follow_redirects=True) as r:
    r.raise_for_status()
    with open(snakemake.output[0], "wb") as f:
        with tqdm(
            desc="Downloading OpenLibrary works dataset",
            total=int(r.headers["Content-Length"]),
            unit_scale=True,
        ) as pbar:
            for chunk in r.iter_bytes():
                f.write(chunk)
                pbar.update(len(chunk))
