# Model Training Methodology

## Training Process Overview

The model training pipeline for predicting book award winners involves the
following key stages:

### 1. **Data Preparation**

The modeling pipeline starts with the processed datasets from the Snakemake
pipeline:

- **Author Data Aggregation**: Multiple authors per book are aggregated by
  taking means for numerical features and lists for categorical features
- **Feature Engineering**:
  - Nomination counts are adjusted to include wins
  - Bestseller statistics are filled with zeros for missing values
  - Publication dates are parsed to extract month information
  - Topicality scores are imputed using year-level means
  - Author categorical features (gender, birth country) are encoded using
    one-hot encoding
  - Numerical features (awards_as_of_year, topicality) are z-scored within year
    cohorts

### 2. **Model Architecture**

The current implementation uses a pipeline approach with the following
components:

- **Data Preprocessing**: Imputation, encoding, and normalization steps
- **Feature Engineering**: Custom transformers for author data encoding and
  cohort-based normalization
- **Classification Models**: Standard binary classification models including:
  - Logistic Regression with cross-validation
  - Decision Trees with hyperparameter tuning
  - Random Forest with hyperparameter tuning
  - XGBoost with hyperparameter tuning

### 3. **Training Approach**

- **Binary Classification**: Models predict whether a book will win any awards
  (n_win > 0)
- **Grouped Cross-Validation**: Uses GroupKFold and StratifiedGroupKFold with
  year as the grouping variable to prevent data leakage
- **Threshold Tuning**: Uses TunedThresholdClassifierCV to optimize
  classification thresholds for F1 score
- **Sample Weighting**: Models are weighted by total nominations (n_nom_all) to
  account for varying nomination counts

### 4. **Evaluation Metrics**

Models are evaluated using:

- **F1 Score**: Primary metric for imbalanced classification
- **Balanced Accuracy**: Secondary metric to account for class imbalance
- **Cross-Validation**: 5-fold grouped cross-validation to ensure robust
  performance estimates

### 5. **Baseline Comparison**

Performance is compared against a naive baseline that randomly assigns wins
based on nomination counts within each year cohort.

## Current Model Performance

Based on the CompareModels notebook, the current best performing model is:

- **Logistic Regression**: F1 score of ~0.50 (doubling the baseline performance
  of ~0.23)
- **Tree-based models**: Similar performance with F1 scores around 0.47-0.49
- **All models**: Significantly outperform the naive baseline

## Feature Importance

The current feature set includes:

- **Author Demographics**: Gender, birth country, age at award
- **Author History**: Cumulative awards won as of the award year
- **Book Characteristics**: Publication month, topicality score
- **Commercial Success**: Bestseller rankings and duration (when available)

## Limitations and Future Directions

The current implementation has several limitations:

- **Binary Classification**: Treats winning as a binary outcome rather than
  modeling the ranking within cohorts
- **Feature Engineering**: Limited interaction between features
- **Temporal Dynamics**: Does not explicitly model how award preferences change
  over time
- **Ensemble Methods**: Could benefit from combining multiple models

See `future_directions.md` for discussion of potential improvements including
learning-to-rank approaches and enhanced topicality modeling.
