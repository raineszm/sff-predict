# Predicting Science Fiction and Fantasy Award Winners

## Project Overview

### Problem Statement

We want to better understand how to pick the winners of literary awards in
science fiction and fantasy from a shortlist of nominees in a given year.

This has potentially both a predictive and inferential component:

- **Predictive**: We will construct a predictive model that assigns an award
  worthiness score to novels published in a given year. When given a cohort of
  such novels the score can be used to assign a probability of winning a
  randomly chosen award for that year.

- **Inferential**: We will also investigate how the distribution of award
  winners has changed over time and try to identify trends in the makeup of the
  award winners.

### Stakeholders

The predictive models would potentially be of use to:

- **Book publishers** looking for up and coming authors
- **Media outlets** trying to decide which works to acquire rights to and where
  to allocate resources

A reasonably predictive model would position stakeholders to make use of a boost
in awareness coming from literary awards.

### Key Performance Indicators

For historical data we could compute an estimated return on investment for
awards predicted by our model based on subsequent media franchises.

## Getting Started

### Prerequisites

- [pixi](https://pixi.sh) for Python environment management

### Setup

1. **Install pixi**

   On macOS you can use [Homebrew](https://brew.sh):
   ```sh
   brew install pixi
   ```
   but other ways of installing are given on the
   [pixi website](https://pixi.sh).

2. **Setup the environment**
   ```sh
   pixi install
   ```

3. **Set up API Keys**

   The project uses the Google Books API to fetch book descriptions. You'll need
   to:

   1. Get an API key from the
      [Google Cloud Console](https://console.cloud.google.com/apis/credentials)
   2. Enable the Books API for your project
   3. Set the API key in one of two ways:
      - Create a `.env` file in the project root with:
        ```
        GOOGLE_BOOKS_API_KEY=your_api_key_here
        ```
      - Or set it as an environment variable:
        ```sh
        export GOOGLE_BOOKS_API_KEY=your_api_key_here
        ```

4. **Run the Pipeline**
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

## Project Approach

This project uses a **cohort-based ranking approach** where books are grouped
into yearly shortlists and trained to predict award winners within each cohort.
The model assigns "award worthiness" scores that can be converted to winning
probabilities.
[Project Overview and Approach](methodology/project_overview.md) - High-level

### Data Processing

For detailed information about the data pipeline, see
[Data Processing Pipeline](methodology/data_pipeline.md).

### Model Training

The model uses a multinomial cross-entropy loss function optimized for winner
selection within shortlist cohorts.

For the complete mathematical methodology and training process, see
[Model Training Methodology](methodology/model_training.md).

## Development

### Environment Management

- All dependencies are managed through pixi
- Add new dependencies with `pixi add`
- Run `pixi install` after updating dependencies

### Workflow

- Raw downloads of data go in `data/raw/`
- Output files are in the `data/` directory
- Add new processing steps to `Snakefile`
- Place scripts in `scripts/` directory
- Place SQL or SPARQL queries in `scripts/queries`
