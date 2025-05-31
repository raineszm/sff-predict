# Model Training Methodology

## Training Process Overview

The model training pipeline for predicting book award winners involves the
following key stages:

### 1. **Data Ingestion**

Multiple sources contribute data:

- **Wikidata**: Provides award nominations and winners, author biographical
  information, and publication metadata
- **OpenLibrary**: Supplies book descriptions and additional metadata
- **Wikipedia**: Offers book descriptions and contextual information
- **Google Books API**: Provides book descriptions and additional metadata
- **World State** adds contextual information relevant to the publication year
  (e.g., cultural or industry trends).

### 2. **Data Preparation**

All inputs are passed through a **Combine, Clean, & Fill Missing** step to
create a unified and complete **Book Data Set**. This step resolves
inconsistencies, integrates fields, and imputes missing values where necessary.

Key features extracted include:

- **Bibliographic data**: Title, publication year, language
- **Author information**: Gender, birth country, age at award
- **Award history**: Previous awards won by the author
- **Book descriptions**: Collected from multiple sources (OpenLibrary,
  Wikipedia, Google Books)

### 3. **Cohort Definition**

Books are grouped into **Yearly Shortlist Cohorts**, representing the set of
nominated or candidate works for each award and year. This framing allows the
model to compare books within each shortlist for winner selection.

### 4. **Modeling Process**

- Each book in a shortlist cohort is assigned an **"Award Worthiness" Score** by
  the model.
- These scores are used to compute probabilities for winning the award within
  the shortlist.
- The model assumes a fixed number of winners per award per year.

### 5. **Loss Computation**

A **Multinomial Cross-Entropy Loss Function** is computed for winner selection:

$$L_{win} = -\sum_i \frac{N_i}{N_\text{win}} \log(\hat{p}_i)$$

where $N_i$ is the observed win count for book $i$, $N_\text{win}$ is the total
number of winners for that cohort, and $\hat{p}_i$ is the predicted probability
of winning. This loss measures the divergence between the observed winner
distribution and the model's predicted probabilities, respecting the fixed
number of winners per shortlist.

## Mathematical Derivation

### Derivation of Scaled Cross-Entropy Loss for Winner Selection

Consider a shortlist cohort of $k$ books with:

- Total winners $N_{\mathrm{win}}$.
- For each book $i$: observed win count $N_i$; predicted probability $\hat p_i$.

#### i. **Multinomial PMF for wins**

Model the allocation of exactly $N_{\mathrm{win}}$ win slots among the $k$
books:

$$P(\{N_i\})
     = \frac{N_{\mathrm{win}}!}{\prod_{i=1}^k N_i!}\;\prod_{i=1}^k \hat p_i^{\,N_i}.$$

This gives the likelihood of observing each book's win count.

#### ii. **Negative log-likelihood**

Transform products into sums:

$$ -\ln P
     = -\ln(N_{\mathrm{win}}!) \;+\; \sum_{i=1}^k \ln(N_i!) \;-\; \sum_{i=1}^k N_i\,\ln(\hat p_i).$$

The first two terms are constant w.r.t.\ the model.

#### iii. **Normalize by total winners**

Divide by $N_{\mathrm{win}}$ to put different years on the same scale:

$$  L_{\mathrm{win}}
     = -\frac{1}{N_{\mathrm{win}}}\sum_{i=1}^k N_i\,\ln(\hat p_i)
     = -\sum_{i=1}^k \frac{N_i}{N_{\mathrm{win}}}\,\ln(\hat p_i).$$

Here $\frac{N_i}{N_{\mathrm{win}}}$ is the empirical win distribution over the
shortlist cohort.

#### iv. **Interpretation**

- **Normalization** by $N_{\mathrm{win}}$ ensures comparability across years.
- **Minimizing** this loss aligns the model's predicted distribution
  $\{\hat p_i\}$ with the empirical distribution $\{N_i/N_{\mathrm{win}}\}$,
  giving more weight to books with more wins.

### 6. **Total Loss Aggregation**

Losses are aggregated **across all shortlist cohorts**, potentially with
weighted contributions based on the significance of awards, to compute the
**Total Loss**.

### 7. **Model Training**

The **Total Loss** guides gradient-based optimization to train the model
parameters.
