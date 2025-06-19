# Baseline Model

## Overview

The baseline model serves as a simple, interpretable comparison point for
evaluating the performance of more sophisticated machine learning approaches. It
implements a nomination-weighted random selection strategy that reflects the
intuitive notion that books with more nominations are more likely to win awards.

## Methodology

### Core Assumption

The baseline model operates on the fundamental assumption that **a book's
probability of winning is proportional to its number of nominations within each
year cohort**. This reflects the common-sense idea that books nominated for more
awards are generally more likely to win at least one award.

### Mathematical Framework

For each year cohort with $k$ books:

- **Total nominations**: $N_{\text{total}} = \sum_{i=1}^k n_{\text{nom},i}$
- **Total winners**: $N_{\text{winners}} = \sum_{i=1}^k n_{\text{win},i}$
- **Winning probability for book $i$**:
  $p_i = \frac{n_{\text{nom},i}}{N_{\text{total}}}$

The model then performs a **multinomial draw** with:

- **Number of trials**: $N_{\text{winners}}$ (total number of wins in the
  cohort)
- **Probabilities**: $\{p_1, p_2, \ldots, p_k\}$ (proportional to nominations)

### Implementation Details

The baseline model is implemented in `models/naive.py` with two key functions:

1. **`naive_cohort_win_counts(df_year)`**: Processes a single year cohort
   - Calculates nomination-based probabilities
   - Performs multinomial sampling
   - Returns win counts for each book

2. **`naive_win_counts(X)`**: Processes the entire dataset
   - Groups data by year
   - Applies cohort processing to each year
   - Returns win predictions for all books

## Evaluation Approach

### Monte Carlo Simulation

Since the baseline model is stochastic (uses random sampling), performance is
evaluated through Monte Carlo simulation:

- **Number of runs**: 1,000 simulations
- **Metrics**: F1 score and balanced accuracy
- **Aggregation**: Mean and standard deviation across runs

### Performance Characteristics

Based on the evaluation in `notebooks/CompareModels.ipynb`:

- **Mean F1 Score**: ~0.23 (23%)
- **Standard Deviation**: ~0.025 (2.5%)
- **Balanced Accuracy**: Similar performance to F1 score

This establishes a clear benchmark: any model achieving significantly better
than 23% F1 score demonstrates meaningful predictive capability.
