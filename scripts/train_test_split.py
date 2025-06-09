import pandas as pd
from snakemake.script import snakemake

# Load the nominee data
df = pd.read_csv(snakemake.input.nominated_novels)

# Load auxiliary data
topicality_scores = pd.read_csv(snakemake.input.topicality_scores)
bestseller_stats = pd.read_csv(snakemake.input.bestseller_stats)

# Merge the auxiliary data into the nominee data
df = df.merge(topicality_scores, on="work_qid", how="left").merge(
    bestseller_stats, on="work_qid", how="left"
)

# Filter to the years we're interested in
# Training set is 1959-2018 (inclusive)
df_train = df[(df["year"] < 2019) & (df["year"] >= 1959)]
# Test set is 2019-2024 (inclusive)
df_test = df[(df["year"] >= 2019) & (df["year"] < 2025)]

# Save the train and test sets
df_train.to_csv(snakemake.output.train_novels, index=False)
df_test.to_csv(snakemake.output.test_novels, index=False)
