# Predicting Science Fiction and Fantasy Award Winners

## Project Overview

### Problem Statement

This project aims to predict winners of science fiction and fantasy literary
awards from shortlisted nominees. Using machine learning models trained on
historical award data, author information, and book characteristics, we
developed a system that assigns award worthiness scores to novels. Our best
model achieves an F1 score of ~50%, more than doubling the baseline performance
of ~23%.

### Stakeholders

Our predictive models would potentially be of use to:

- **Book publishers** looking for up and coming authors
- **Media outlets** trying to decide which works to acquire rights to and where
  to allocate resources

A reasonably predictive model would position stakeholders to make use of a boost
in awareness coming from literary awards.

### Key Performance Indicators

For our KPI we compute the balanced accuracy and F1 score of predicted winners.

## Current Results

Our current implementation achieves:

- **F1 Score**: ~50% with logistic regression models
- **Performance Improvement**: More than doubles the baseline performance of
  ~23%
- **Multiple Models**: Logistic regression, decision trees, random forests, and
  XGBoost all show similar performance

For detailed model comparisons and performance analysis, see the
`notebooks/CompareModels.ipynb` notebook.

## Technical Approach

### Data Sources

Our model combines multiple data sources to predict award success:

- **Award Data**: Historical nominations and winners from major science fiction
  and fantasy awards
- **Author Information**: Demographics, previous awards, and career trajectory
- **Book Descriptions**: Text descriptions from multiple sources (OpenLibrary,
  Wikipedia, Google Books)
- **Commercial Success**: New York Times bestseller data when available
- **Cultural Context**: Historical news headlines for topicality analysis

### Machine Learning Pipeline

1. **Data Processing**: Complex ETL pipeline using Snakemake to orchestrate data
   collection, cleaning, and feature engineering
2. **NLP Processing**: Book descriptions and news headlines are embedded using
   advanced language models
3. **Topicality Scoring**: Books are scored for how well their themes align with
   contemporaneous news headlines
4. **Binary Classification**: Models predict whether a book will win any awards
   using standard classification algorithms
5. **Robust Evaluation**: Training uses year-based grouping to prevent data
   leakage and ensure reliable performance estimates

For detailed technical information, see the methodology documentation:

- [Project Overview](methodology/project_overview.md) - High-level approach and
  current performance
- [Data Pipeline](methodology/data_pipeline.md) - Detailed data processing
  workflow
- [Model Training](methodology/model_training.md) - Current modeling approach
  and results
- [Future Directions](methodology/future_directions.md) - Potential improvements
  and extensions

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

   Similarly you'll need to use the New York Time API to get the headlines

   1. Get an API key from the
      [New York Times Developer Portal](https://developer.nytimes.com)
   2. Enable the Article Search API for your project
   3. Set the API key in one of two ways:
      - Create (or add to) a `.env` file in the project root with:
        ```
        NYTIMES_API_KEY=your_api_key_here
        ```
      - Or set it as an environment variable:
        ```sh
        export NYTIMES_API_KEY=your_api_key_here
        ```

4. **Fetch Cached Data** (Optional)

   Some intermediate data files are cached to speed up pipeline runs. To fetch
   them:

   1. First, log in to Kaggle:
      ```sh
      pixi run kaggle authenticate
      ```
      You'll need a Kaggle account and API credentials from your
      [Kaggle account settings](https://www.kaggle.com/settings)

   2. Then download the cached data:
      ```sh
      bin/fetch_cache.sh
      ```

   This will save you from having to reprocess the full headlines dataset, which
   requires NYT API credentials and takes significant time.

   If at any point snakemake tries to regenerate the cached files, you can reset
   their timestamps by running `bin/fetch_cache.sh` again.

5. **Run the Pipeline**
   ```sh
   bin/snakemake.sh
   ```

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

### Key Notebooks

#### Data Exploration

- [Exploratory data analysis](notebooks/EDA.ipynb)
- [Topicality analysis and visualization](notebooks/Topicality.ipynb)

#### Model exploration

- [Logistic regression modeling](notebooks/Logistic.ipynb)
- [Tree-based model analysis](notebooks/Trees.ipynb)
- [Linear regression analysis](notebooks/LinearModels.ipynb)

#### Evaluation

- [Model comparison and performance analysis](notebooks/CompareModels.ipynb)
- [Evaluation of model performance on test set](notebooks/ModelPerformance.ipynb)
