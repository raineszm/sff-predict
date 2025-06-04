#!/usr/bin/env bash

set -euo pipefail

echo "Fetching prefetched NYT headlines and embeddings..."
TMPDIR="$(mktemp -d)"
if [[ ! -f "data/raw/headlines.json" ]] || [[ ! -f "data/headline_embeddings.parquet" ]]; then
    pixi run kaggle datasets download -d zacharymraines/sff-project-data-cache --path "$TMPDIR" --unzip
    mv "$TMPDIR/headlines.json" "data/raw/headlines.json" 
    mv "$TMPDIR/headline_embeddings.parquet" "data/headline_embeddings.parquet"
fi

echo "Generating monthly headlines, so that Snakemake is happy"
pixi run python -m scripts.split_cached_headlines
for file in data/raw/headlines/*.json; do
    touch "$file"
done
touch "data/raw/headlines.json"
touch "data/headline_embeddings.parquet"

echo "Done!"
