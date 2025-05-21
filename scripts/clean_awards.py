from snakemake.script import snakemake
import pandas as pd
import html

# Let's clean the awards data
# so we can match it with the books data
df = pd.read_json(snakemake.input[0])

# decode HTML entities
df["title"] = df.title.map(html.unescape)

# some titles were originally written as a short story
# and then later a book
# this is written as \u201c...\u201d (book title ...)
# or
# \u201c...\u201d (expanded as ...)
# we use the book title
df["title"] = df.title.str.replace(
    r"\u201c.*\u201d\s*\((?:expanded as|book title)\s+([^)]*)\)", r"\1", regex=True
).str.rstrip()

# and remove the quotes
df["title"] = df.title.str.replace(r"\u201c(.*)\u201d", r"\1", regex=True).str.rstrip()

# some titles include trailing parentheticals
# e.g.
# (by [author])
# or (series title ...)
# we remove these to normalize the titles
df["title"] = df.title.str.replace(r"(\([^)]+\)\s*)+\s*$", "", regex=True)

# Replace hyphens with spaces
df["title"] = df.title.str.replace(r"\b-\b", " ", regex=True)

# Collapse multiple hyphens into a single hyphen
df["title"] = df.title.str.replace(r"-+", "-", regex=True)

# and strip whitespace for good measure
df["title"] = df.title.str.strip()


# collapse books with multiple authors
df = (
    df.groupby(["title", "year", "award", "result"])
    .nominee.apply(lambda x: ", ".join(x.sort_values()))
    .reset_index()
)
df.to_json(snakemake.output[0], orient="records", lines=True)
