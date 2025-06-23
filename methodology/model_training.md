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
    one-hot encoding and then aggregated over nominated novel
  - Numerical features (awards_as_of_year, topicality) are z-scored within year
    cohorts

### 2. **Model Architecture**

- **Data Preprocessing**: Imputation, encoding, and normalization steps
- **Feature Engineering**: Custom transformers for author data encoding and
  cohort-based normalization

We implemented two distinct modeling approaches with different architectures:

#### 2.1 **Binary Classification Architecture**

The binary classification pipeline uses the following components:

- **Classification Models**: Standard binary classification models including:
  - Logistic Regression with cross-validation
  - Decision Trees with hyperparameter tuning
  - Random Forest with hyperparameter tuning
  - XGBoost with hyperparameter tuning

#### 2.2 **Regression Architecture**

The regression pipeline for predicting win counts uses:

- **Regression Models**: Linear regression models with regularization including:
  - Linear Regression
  - Ridge Regression with cross-validated alpha selection
  - Lasso Regression with cross-validated alpha selection
  - Naive baseline model using nomination-based predictions (basline)

### 3. **Training Approach**

We approached modeling in two distinct ways:

#### 3.1 **Binary Classification Training**

- **Binary Classification**: Models predict whether a book will win any awards
  (n_win > 0)
- **Grouped Cross-Validation**: Uses GroupKFold and StratifiedGroupKFold with
  year as the grouping variable to prevent data leakage
- **Threshold Tuning**: Uses TunedThresholdClassifierCV to optimize
  classification thresholds for F1 score
- **Sample Weighting**: Models are weighted by total nominations (n_nom_all) to
  account for varying nomination counts

#### 3.2 **Regression Training**

- **Win Count Prediction**: Models predict the actual number of awards won
  (n_win as continuous variable)
- **Cross-Validation**: Uses 5-fold cross-validation for model evaluation
- **Regularization**: Ridge and Lasso models use cross-validated alpha selection
  to prevent overfitting
- **Alpha Optimization**: Ridge uses RidgeCV with alphas from 10^-4 to 10^4,
  Lasso uses LassoCV with automatic alpha selection
- **Baseline Comparison**: Includes naive model that predicts wins based on
  nomination patterns within year cohorts

### 4. **Evaluation Metrics**

Models are evaluated using different metrics appropriate for each approach:

#### 4.1 **Binary Classification Metrics**

- **F1 Score**: Primary metric for imbalanced classification
- **Balanced Accuracy**: Secondary metric to account for class imbalance
- **Cross-Validation**: 5-fold grouped cross-validation to ensure robust
  performance estimates

#### 4.2 **Regression Metrics**

- **Root Mean Squared Error (RMSE)**: Primary metric for regression performance
- **R-squared (R²)**: Proportion of variance explained by the model
- **Cross-Validation**: 5-fold cross-validation to ensure robust performance
  estimates

### 5. **Baseline Comparison**

Performance is compared against appropriate baselines for each approach:

#### 5.1 **Binary Classification Baseline**

Performance is compared against a naive baseline that randomly assigns wins
based on nomination counts within each year cohort.

#### 5.2 **Regression Baseline**

Performance is compared against:

- **Naive Model**: Predicts win counts proportional to the nomination patterns
  within a year cohorts (R² ≈ -0.56, RMSE ≈ 0.59)

## Current Model Performance

Based on the CompareModels notebook and linear_ridge_lasso notebook, the current
best performing models are:

#### Binary Classification Performance

- **Logistic Regression**: F1 score of ~0.50 (doubling the baseline performance
  of ~0.23)
- **Tree-based models**: Similar performance with F1 scores around 0.47-0.49
- **All models**: Significantly outperform the naive baseline

#### Regression Performance

- **Linear Models**: All linear models show similar performance:
  - **Linear Regression**: R² ≈ 0.30, RMSE ≈ 0.41
  - **Ridge Regression**: R² ≈ 0.29, RMSE ≈ 0.40
  - **Lasso Regression**: R² ≈ 0.29, RMSE ≈ 0.40
- **All linear models**: Significantly outperform the naive baseline (R² ≈
  -0.56)

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
- **Linear Assumptions**: Linear models assume linear relationships between
  features and win counts
- **Feature Interactions**: Limited modeling of feature interactions
- **Temporal Modeling**: No explicit modeling of temporal trends
- **Non-linear Relationships**: Linear models may miss non-linear feature
  relationships that could be captured by tree-based regression models

See `future_directions.md` for discussion of potential improvements including
learning-to-rank approaches and enhanced topicality modeling.
