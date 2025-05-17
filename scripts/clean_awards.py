from snakemake.script import snakemake
import pandas as pd

df = pd.read_json(snakemake.input[0])

# some titles include (by [author])
# but the author isn't always written the same
# we remove this since it's redundant with the nominee field anyway
df["title"] = df.title.str.replace(r"\(by [^)]+\)", "", regex=True)

# some titles were originally written as a short story
# and then later a book
# this is written as \u201c...\u201d (book title ...)
# or
# \u201c...\u201d (expanded as ...)
# we use the book title
df["title"] = df.title.str.replace(
    r"\u201c(.*)\u201d.*", r"\1", regex=True
).str.rstrip()


# collapse books with multiple authors
df = (
    df.groupby(["title", "year", "award", "result"])
    .nominee.apply(lambda x: ", ".join(x.sort_values()))
    .reset_index()
)
df.to_json(snakemake.output[0], orient="records", lines=True)
