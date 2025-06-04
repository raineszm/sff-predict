#!/usr/bin/env bash

set -euo pipefail

echo "Fetching prefetched NYT headlines and embeddings..."
TMPDIR="$(mktemp -d)"
if [[ ! -f "data/raw/headlines.json" ]]; then
	pixi run kaggle datasets download -d zacharymraines/nyt-headlines-1954-2024 --path "$TMPDIR" --unzip
    echo "Generating monthly headlines, so that Snakemake is happy"
	mv "$TMPDIR/headlines.json" "data/raw/headlines.json"
    pixi run python -m scripts.split_cached_headlines
fi

if [[ ! -f "data/headline_embeddings.parquet" ]]; then
	pixi run kaggle datasets download -d zacharymraines/nyt-headline-embeddings-1954-2024 --path "$TMPDIR" --unzip
	mv "$TMPDIR/headline_embeddings.parquet" "data/headline_embeddings.parquet"
fi


bin/snakemake.sh --quiet all --touch data/raw/headlines/*.json data/raw/headlines.json data/headline_embeddings.parquet

echo "Done!"
