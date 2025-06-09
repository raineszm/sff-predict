import pandas as pd
import spacy
from unidecode import unidecode
import splink.comparison_library as cl
import splink.blocking_rule_library as brl
from splink import SettingsCreator, Linker, DuckDBAPI
from snakemake.script import snakemake
from spacy.tokens import Doc
from spacy.language import Language


def normalize_doc(doc: Doc):
    """
    Normalize a spaCy document by lemmatizing non-stop words longer than 2 characters.

    Args:
        doc: A spaCy document to normalize

    Returns:
        A string of lemmatized words joined by spaces
    """
    return " ".join([tok.text for tok in doc if not tok.is_stop and len(tok) > 2])


def normalized_column(s: pd.Series, nlp: Language):
    """
    Normalize a pandas Series of text using spaCy.

    Args:
        s: A pandas Series of text to normalize
        nlp: A spaCy language model

    Returns:
        A list of normalized strings
    """
    return list(map(normalize_doc, nlp.pipe(s.str.casefold().apply(unidecode))))


# Load spaCy model
nlp = spacy.load("en_core_web_sm", disable=["parser", "ner"])

best_sellers = pd.read_table(snakemake.input.bestsellers)
novels = pd.read_csv(snakemake.input.nominated)

# Load up bestsellers for the years we're interested in
# and produce normalized versions of the title and author
# to compare on
sff_bestsellers = (
    (
        best_sellers.drop_duplicates(subset=["title_id"], keep="first").pipe(
            lambda df: df[df.year >= 1959]
        )
    )[["title_id", "title", "author", "year"]]
    .assign(
        title_key=lambda df: normalized_column(df.title, nlp),
        author_key=lambda df: normalized_column(df.author, nlp),
    )
    .rename(columns={"title_id": "unique_id"})
)

# Do the same for our nominated novels
sff_novels = (
    novels.groupby("work_qid")
    .agg(
        {
            "authorLabel": lambda x: " ".join(x.unique()),
            "title": "first",
            "year": "first",
        }
    )
    .assign(
        title_key=lambda x: normalized_column(x.title, nlp),
        author_key=lambda x: normalized_column(x.authorLabel, nlp),
    )
    .reset_index()
    .rename(
        columns={
            "work_qid": "unique_id",
            "authorLabel": "author",
        }
    )
)

# A custom rule that says that the year on the bestseller
# list must be the year before the award year
previous_year_rule = brl.CustomRule("l.year = r.year - 1")

# Configure splink for match
settings = SettingsCreator(
    link_type="link_only",
    blocking_rules_to_generate_predictions=[
        previous_year_rule,
    ],
    comparisons=[
        cl.JaroWinklerAtThresholds("title_key"),
        cl.JaroWinklerAtThresholds("author_key"),
    ],
    retain_matching_columns=True,
)

# Initialize and train the linker
linker = Linker([sff_bestsellers, sff_novels], settings, db_api=DuckDBAPI())

linker.training.estimate_u_using_random_sampling(max_pairs=200_000)
linker.training.estimate_parameters_using_expectation_maximisation(previous_year_rule)

pairs = linker.inference.predict(threshold_match_probability=0.9)


# Obtain a join table for the matched records
join_table = (
    pairs.as_pandas_dataframe()
    .rename(
        columns={
            "unique_id_l": "title_id",
            "unique_id_r": "work_qid",
            "year_r": "award_year",
        }
    )
    .drop_duplicates(subset=["title_id", "work_qid"])[
        ["title_id", "work_qid", "award_year"]
    ]
    .assign(title_id=lambda x: x.title_id.astype("int"))
)

# Calculate some bestseller statistics
bestseller_stats = (
    join_table.merge(
        best_sellers[["title_id", "year", "rank"]], on="title_id", how="left"
    )
    .pipe(lambda x: x[x.year < x.award_year])
    .groupby(["title_id", "work_qid", "year"])["rank"]
    .agg(["max", "count", "median"])
    .reset_index()[["work_qid", "max", "count", "median"]]
)

bestseller_stats.to_csv(snakemake.output.bestseller_stats, index=False)
