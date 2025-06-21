# Data Processing Pipeline

The data processing pipeline uses Snakemake to orchestrate the following steps:

## Pipeline Overview

### 1. Data Download

- **Wikidata Data**: Fetches award information from Wikidata using SPARQL
  queries at <https://query.wikidata.org>, including:
  - Award nominations and winners for major science fiction and fantasy awards
  - Author biographical information (gender, birth country, date of birth)
  - Publication metadata (publication dates, language)
  - Cross-referenced identifiers (ISBN, OpenLibrary, Goodreads)

- **OpenLibrary Data**: Downloads work metadata from OpenLibrary, including:
  - Book descriptions
  - Publication information
  - Author details

- **New York Times Bestsellers**: Downloads bestseller data from Kaggle,
  including:
  - Bestseller rankings and duration on lists
  - Publication dates and sales performance

- **New York Times Headlines**: Fetches historical headlines from NYT API to
  provide cultural context for topicality analysis

### 2. Data Processing

The pipeline processes the raw data to:

- Clean and normalize award information
- Match award entries with author data
- Generate cumulative award statistics
- Extract and organize book descriptions from multiple sources
- Match books with bestseller data
- Compute topicality scores using NLP embeddings
- Create analysis-ready datasets

### 3. NLP Processing

The pipeline includes several NLP processing steps:

- **Description Embeddings**: Book descriptions are embedded using the
  `all-MiniLM-L6-v2` model
- **Headline Embeddings**: NYT headlines are embedded for cultural context
- **Debiasing**: Description embeddings are debiased using nullspace projection
- **Topicality Computation**: Books are scored for topicality by comparing their
  embeddings with contemporaneous news headlines

### 4. Output Files

The pipeline produces several key datasets:

- `data/nominated_novels.csv`: Cleaned award nomination data with author
  information, including:
  - Title, author, and publication year
  - Award nominations and wins
  - Author biographical information (age at award, gender, birth country)
  - Previous award counts for authors

- `data/cumulative_awards.csv`: Running totals of awards per author over time

- `data/descriptions.csv`: Book descriptions collected from multiple sources:
  - OpenLibrary
  - Wikipedia
  - Google Books (via ISBN)
  - Google Books (via title and author search)

- `data/description_embeddings.parquet`: Sentence embeddings of book
  descriptions

- `data/descriptions_debiased.parquet`: Debiased embeddings for fairer analysis

- `data/headline_embeddings.parquet`: Embeddings of NYT headlines for topicality
  analysis

- `data/topicality_scores.csv`: Computed topicality scores for each book

- `data/bestseller_stats.csv`: Bestseller performance statistics for matched
  books

- `data/train_novels.csv` and `data/test_novels.csv`: Final datasets ready for
  modeling

## Data Files

### Raw Files

Located in `data/raw/`:

**These are not included in the repo, they are generated when the snakemake
workflow is run**

- `award_novels.json`: Downloaded from Wikidata using SPARQL query, contains
  information about novels nominated for major science fiction and fantasy
  awards, including:
  - Title, author, and publication year
  - Award nominations and wins
  - OpenLibrary and ISBN identifiers
  - Publication dates and language information
  - **NB** This includes only books nominated for one of the major scifi/fantasy
    awards in the best novel category (see
    `scripts/queries/wikidata_awards.sparql` for which awards are included)

- `awards_all.json`: Downloaded from Wikidata, contains comprehensive award
  information including:
  - All winners and nominees across multiple award categories
  - Award dates and categories
  - Author information
  - Covers major awards like Hugo, Nebula, World Fantasy, and more
  - **NB** This includes _winners_ in _any_ category and for a broader set of
    awards (see `scripts/queries/wikidata_winners.sparql`)

- `nominee_biographical.json`: Downloaded from Wikidata, contains biographical
  information about award nominees:
  - Author names and Wikidata IDs
  - Gender information
  - Date of birth
  - Birth country
  - Residence information

- `ol_works.parquet`: Downloaded from OpenLibrary, contains:
  - Work metadata including descriptions
  - Publication information
  - Author details

- `nyt_full.tsv`: Downloaded from Kaggle, contains:
  - New York Times bestseller data
  - Rankings and duration on lists
  - Publication dates

- `headlines.json`: Downloaded from NYT API, contains:
  - Historical headlines from 1954-2024
  - Used for topicality analysis

### Processed Files

Located in `data/`:

- `nominated_novels.csv`: Processed data from `award_novels.json`, containing:
  - Cleaned and normalized award nomination data
  - Counts of nominations and wins per work
  - Matched identifiers for cross-referencing
  - Publication information
  - Author biographical information
  - Previous award counts for authors

- `cumulative_awards.csv`: Generated from `awards_all.json`, containing:
  - Running total of awards per author as of year
  - Example: As of the end of `2023`, `Ann Leckie` had won awards listed as
    major awards on sfadb `14` times.

- `descriptions.csv`: Book descriptions collected from multiple sources:
  - OpenLibrary
  - Wikipedia
  - Google Books (via ISBN)
  - Google Books (via title and author search)

- `description_embeddings.parquet`: Sentence embeddings of book descriptions
  using all-MiniLM-L6-v2

- `descriptions_debiased.parquet`: Debiased embeddings using nullspace
  projection

- `headline_embeddings.parquet`: Embeddings of NYT headlines for topicality
  analysis

- `topicality_scores.csv`: Computed topicality scores comparing book embeddings
  with contemporaneous headlines

- `bestseller_stats.csv`: Bestseller performance statistics for books that
  appear on NYT lists

- `train_novels.csv` and `test_novels.csv`: Final datasets with all features
  merged and ready for modeling

## Running the Pipeline

```sh
# Run with a single core
pixi run snakemake --cores 1

# Or use multiple cores for faster processing
pixi run snakemake --cores 4
```

There's a script to make this slightly simpler:

```sh
bin/snakemake.sh
```

## Adding New Data Processing Steps

- Add new rules to `Snakefile`
- Place scripts in `scripts/` directory
- Place SQL or SPARQL queries in `scripts/queries`
