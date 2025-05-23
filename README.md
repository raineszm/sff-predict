# Predicting Scifi and Fantasy Award Winners

## Development Setup

### Prerequisites
- [pixi](pixi.sh) for Python environment management

### Getting Started

1. **Install pixi**

   On macOS you can use [Homebrew](brew.sh):
   ```sh
   brew install pixi
   ```
   but other ways of installing are given on the [pixi website](pixi.sh).

2. **Clone and Setup**
   ```sh
   git clone <repository-url>
   cd scifi-fantasy
   pixi install
   ```

3. **Run the Pipeline**
   ```sh
   # Run with a single core
   pixi run snakemake --cores 1

   # Or use multiple cores for faster processing
   pixi run snakemake --cores 4
   ```

### Development Workflow

1. **Environment Management**
   - All dependencies are managed through pixi
   - Add new dependencies with `pixi add`
   - Run `pixi install` after updating dependencies

2. **Data Processing**
   - Configuration is in `config.yaml`
   - Adjust thresholds and paths as needed
   - Raw downloads of data go in `data/raw/`
   - Output files are in the `data/` 

3. **Adding New Data Processing steps**
   - Add new rules to `Snakefile`
   - Place scripts in `scripts/` directory

## Data Processing Pipeline

The data processing pipeline uses Snakemake to orchestrate the following steps:

### 1. Data Download
- **Goodreads Data**: Downloads raw book, author, and works data from the Goodreads dataset <https://cseweb.ucsd.edu/~jmcauley/datasets/goodreads.html>
- **Awards Data**: Fetches award information from Wikidata using SPARQL queries <https://query.wikidata.org>

### 2. Data Combination and Filtering
The `combine_data` rule processes the raw data using DuckDB to:
- Filter works based on:
  - Minimum number of ratings (`ratings_threshold`)
  - Minimum number of scifi/fantasy tags (`tag_threshold`)
- Generate three output files:
  - `selected_works_data`: Filtered works meeting the criteria
  - `augmented_works_data`: Selected works with additional metadata
  - `identifiers_data`: Mapping of various identifiers (ISBN, Goodreads IDs, etc.)

### 3. Awards Processing
- **Pivot Awards**: Transforms raw awards data into a more analysis-friendly format
- **Clean Awards**: Matches award entries with works using the identifiers data

### Output Files
The pipeline produces several key datasets:
- `data/sff_works.parquet`: Filtered science fiction and fantasy works
- `data/sff_works_augmented.parquet`: Works with additional metadata
- `data/cleaned_awards_data.parquet`: Processed award information linked to works

## Project Brainstorming

The basic idea is to try break data into years and train a ranking model on each year independently.

Schematically the idea can be depicted like the following

![Data flow diagram showing how books are grouped into yearly cohorts and scored for award worthiness](DataFlow.png)

### Potential Features
- Bibliographic data
- Publisher information
- Author information
- World events/trends for the award year
- Proxies for sales information

## Model Training Process for Predicting Award-Winning Books

The model training pipeline for predicting book award winners involves the following key stages:

### 1. **Data Ingestion**
Multiple sources contribute data:
- **Goodreads**, **Open Library**, **Library of Congress**, **Bibliographic Data**, and **Wikidata** provide rich metadata on books.
- **Award Winners & Nominees** supplies ground truth labels for supervised learning.
- **World State** adds contextual information relevant to the publication year (e.g., cultural or industry trends).

### 2. **Data Preparation**
All inputs are passed through a **Combine, Clean, & Fill Missing** step to create a unified and complete **Book Data Set**. This step resolves inconsistencies, integrates fields, and imputes missing values where necessary.

### 3. **Cohort Definition**
Books are grouped into **Yearly Cohorts**, representing the set of eligible works for each award year. This framing allows the model to compare books within a constrained competition set.

### 4. **Modeling Process**
- Each book in a cohort is assigned an **"Award Worthiness" Score** by the model.
- These scores are used to compute separate probabilities for winning and nomination:
  - **Win Probability** $p_i$ for each book being selected as a winner
  - **Nomination Probability** $q_i$ for each book being nominated
- Both probabilities are modeled using a multinomial distribution, reflecting the fixed number of winners and nominees each year.
  - **Fixed Number of Wins**: Typically, there is exactly 1 winner per award per year.
  - **Fixed Number of Nominations**: The number of nominees is predetermined, often 5-6 for major awards.

### 5. **Loss Computation**
Separate **Multinomial Cross-Entropy Loss Functions** are computed for win and nomination probabilities:

- **Win Loss**: $$L_{win} = -\sum_i \hat{p}_i \log(p_i)$$
- **Nomination Loss**: $$L_{nom} = -\sum_i \hat{q}_i \log(q_i)$$

These measure the divergence between predicted and actual probabilities for winners and nominees, respectively, while respecting the fixed number of winners and nominees each year.

#### Derivation of Cross-Entropy Loss

The cross-entropy loss can be derived from the negative log-likelihood of the observed number of wins, assuming a multinomial distribution with probabilities $p_i$ and $n = N_{\text{win}}$:

The probability mass function for a multinomial distribution is:

$$ 
P(X_1 = x_1, X_2 = x_2, \ldots, X_k = x_k) = \frac{n!}{x_1! x_2! \cdots x_k!} p_1^{x_1} p_2^{x_2} \cdots p_k^{x_k}
$$

The negative log-likelihood is:

$$ 
LL = -\log P(X_1 = x_1, X_2 = x_2, \ldots, X_k = x_k) = -\log \left( \frac{n!}{x_1! x_2! \cdots x_k!} p_1^{x_1} p_2^{x_2} \cdots p_k^{x_k} \right)
$$

This simplifies to:

$$
LL = -\log(n!) + \sum_{i=1}^{k} \log(x_i!) - \sum_{i=1}^{k} x_i \log(p_i)
$$

In machine learning, we focus on the terms that depend on the parameters $p_i$, leading to the cross-entropy loss:

$$
L = -\sum_{i=1}^{k} x_i \log(p_i)
$$

This loss function measures the divergence between the observed distribution (given by $x_i$) and the predicted probabilities $p_i$.

### 6. **Total Loss Aggregation**
Losses are aggregated **across all cohorts**, potentially with weighted contributions based on the significance of awards or nominations, to compute the **Total Loss**.

### 7. **Model Training**
The **Total Loss** guides gradient-based optimization to train the model parameters.