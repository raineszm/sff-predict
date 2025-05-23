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

- **Win Loss**: $$L_{win} = -\sum_i \frac{N_i}{N_\text{win}} \log(\hat{p}_i)$$
- **Nomination Loss**: $$L_{nom} = -\sum_i \frac{M_i}{N_\text{nom}} \log(\hat{q}_i)$$

where $N_i$ and $M_i$ are the observed counts (how many wins or nominations for that book respectively), $N_\text{win}$ and $N_\text{nom}$ are the total number of wins and nominations for that cohort, and $\hat{p}_i$ and $\hat{q}_i$ are the predicted probabilities. These measure the divergence between the observed outcomes and predicted probabilities while respecting the fixed number of winners and nominees each year.


#### Derivation of Scaled Cross-Entropy Loss for Award Predictions

Consider a cohort of $k$ books with:
- Total winners $N_{\mathrm{win}}>1$ and total nominees $N_{\mathrm{nom}}>1$.  
- For each book $i$: observed win count $N_i$, nomination count $M_i$; predicted probabilities $\hat p_i$ (win) and $\hat q_i$ (nomination).

i. **Multinomial PMF for wins**  
   Model the allocation of exactly $N_{\mathrm{win}}$ win slots among the $k$ books:
   $$
     P(\{N_i\})
     = \frac{N_{\mathrm{win}}!}{\prod_{i=1}^k N_i!}\;\prod_{i=1}^k \hat p_i^{\,N_i}.
   $$
   This gives the likelihood of observing each book's win count.

ii. **Negative log-likelihood**  
   Transform products into sums:
   $$
     -\ln P
     = -\ln(N_{\mathrm{win}}!) \;+\; \sum_{i=1}^k \ln(N_i!) \;-\; \sum_{i=1}^k N_i\,\ln(\hat p_i).
   $$
   The first two terms are constant w.r.t.\ the model.

iii. **Normalize by total winners**  
   Divide by $N_{\mathrm{win}}$ to put different years on the same scale:
   $$
     L_{\mathrm{win}}
     = -\frac{1}{N_{\mathrm{win}}}\sum_{i=1}^k N_i\,\ln(\hat p_i)
     = -\sum_{i=1}^k \frac{N_i}{N_{\mathrm{win}}}\,\ln(\hat p_i).
   $$
   Here $\frac{N_i}{N_{\mathrm{win}}}$ is the empirical win distribution over the cohort.

iv. **Nomination loss**  
   Analogously, with $N_{\mathrm{nom}}$ fixed nominations:
   $$
     L_{\mathrm{nom}}
     = -\sum_{i=1}^k \frac{M_i}{N_{\mathrm{nom}}}\,\ln(\hat q_i).
   $$

v. **Interpretation**  
   - **Normalization** by $N_{\mathrm{win}}$ or $N_{\mathrm{nom}}$ ensures comparability across years.  
   - **Minimizing** these losses aligns the model's predicted distributions $\{\hat p_i\}$, $\{\hat q_i\}$ with the empirical distributions $\{N_i/N_{\mathrm{win}}\}$, $\{M_i/N_{\mathrm{nom}}\}$, giving more weight to books with more wins or nominations.

### 6. **Total Loss Aggregation**
Losses are aggregated **across all cohorts**, potentially with weighted contributions based on the significance of awards or nominations, to compute the **Total Loss**.

### 7. **Model Training**
The **Total Loss** guides gradient-based optimization to train the model parameters.