import polars as pl
from snakemake.script import snakemake

RATINGS_THRESHOLD = snakemake.config["filter"]["ratings_threshold"]
TAG_THRESHOLD = RATINGS_THRESHOLD * snakemake.config["filter"]["tag_threshold"]
SCIFI_FANTASY_REGEX = r"(?i)fantasy|sci(ence)?-fi"

print("Selecting works with at least {} ratings".format(RATINGS_THRESHOLD))
raw_works = (
    pl.scan_ndjson(snakemake.input["works"])
    .filter(
        pl.col("ratings_count").cast(int) >= RATINGS_THRESHOLD,
        pl.col("original_publication_year").cast(int, strict=False) >= 1959,
    )
    .collect(engine="streaming")
)

books = pl.scan_ndjson(snakemake.input["books"], infer_schema_length=1000).filter(
    pl.col("language_code").str.starts_with("en"),
    pl.col("work_id").is_in(raw_works["work_id"].to_list()),
    pl.col("publication_year").cast(int, strict=False) >= 1959,
)

print("Collecting tags for selected works")
work_tags = (
    books.group_by("work_id")
    .agg(
        pl.col("popular_shelves").flatten().alias("tags"),
    )
    .filter(
        pl.col("tags")
        .list.eval(
            pl.element().struct.field("name").str.contains(SCIFI_FANTASY_REGEX)
            & (pl.element().struct.field("count").cast(int) >= TAG_THRESHOLD)
        )
        .list.any()
    )
).collect(engine="streaming")

works = raw_works.join(work_tags, on="work_id", how="inner")

print(
    "Writing filtered works and tags to {}".format(snakemake.output["filtered_works"])
)
works.write_parquet(snakemake.output["filtered_works"])

print("Combining into one dataset")
complete_works = works.lazy().join(
    books, left_on="best_book_id", right_on="book_id", how="inner"
)

print(
    "Writing combined works and tags to {}".format(snakemake.output["combined_works"])
)
complete_works.sink_parquet(snakemake.output["combined_works"])
