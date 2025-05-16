#!/usr/bin/env python3
import requests
from tqdm.auto import tqdm
from snakemake.script import snakemake
import os
import shutil


def download_file(url, output_path):
    """Download a file with progress bar using requests and tqdm."""
    # Stream the download
    response = requests.get(url, stream=True)
    response.raise_for_status()  # Raise an exception for bad status codes

    # Get the total file size
    total_size = int(response.headers.get("content-length", 0))

    # Download with progress bar
    with open(output_path, "wb") as f:
        with tqdm.wrapattr(
            response.raw,
            "read",
            total=total_size,
            desc=f"Downloading {os.path.basename(output_path)}",
            unit="iB",
        ) as wrapped:
            shutil.copyfileobj(wrapped, f)


download_file(snakemake.params["url"], snakemake.output[0])
