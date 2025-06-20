# Project Overview and Approach

## Current Features

Our model considers multiple types of features to predict award worthiness:

- **Bibliographic data**: Publication details, length, genre classifications
- **Author information**: Demographics (gender, birth country, age), previous
  awards, career trajectory
- **Commercial success**: Bestseller rankings and duration on lists (when
  available)
- **Topicality**: NLP-derived scores measuring how well book themes align with
  contemporaneous news headlines
- **Award history**: Cumulative awards won by authors as of the award year

## Data Sources

The model training pipeline ingests data from multiple sources:

- **Wikidata**: Provides award nominations and winners, author biographical
  information, and publication metadata
- **OpenLibrary**: Supplies book descriptions and additional metadata
- **Wikipedia**: Offers book descriptions and contextual information
- **Google Books API**: Provides book descriptions and additional metadata
- **New York Times API**: Historical headlines for topicality analysis
- **Kaggle**: NYT bestseller data for commercial success metrics

## Current Approach Summary

1. **Data Processing**: Complex ETL pipeline using Snakemake to orchestrate data
   collection, cleaning, and feature engineering
2. **NLP Processing**: Book descriptions and news headlines are embedded using
   sentence transformers, with debiasing applied to description embeddings
3. **Topicality Scoring**: Books are scored for topicality by comparing their
   embeddings with contemporaneous news headlines
4. **Binary Classification**: Models predict whether a book will win any awards
   (n_win > 0) using standard classification algorithms
5. **Grouped Cross-Validation**: Training uses year-based grouping to prevent
   data leakage and ensure robust performance estimates
6. **Feature Engineering**: Custom transformers handle author data aggregation,
   categorical encoding, and cohort-based normalization

## Current Performance

The current implementation achieves:

- **F1 Score**: ~0.50 with logistic regression (doubling the baseline
  performance of ~0.23)
- **Robust Evaluation**: 5-fold grouped cross-validation with year-based splits
- **Multiple Models**: Logistic regression, decision trees, random forests, and
  XGBoost all show similar performance

## Future Directions

While the current implementation provides a solid foundation, several
improvements are planned:

- **Learning-to-Rank**: Implement ranking-based approaches that model the
  relative ordering within award cohorts
- **Enhanced Topicality**: Improve topicality scoring with temporal dynamics and
  content-specific analysis
- **Deep Learning**: Explore neural network and transformer-based approaches
- **Additional Data Sources**: Incorporate book reviews, social media, and
  academic citations

See `future_directions.md` for detailed discussion of potential improvements.
