import polars as pl
from snakemake.script import snakemake
import re

RATINGS_THRESHOLD = snakemake.config["filter"]["ratings_threshold"]
TAG_THRESHOLD = RATINGS_THRESHOLD * snakemake.config["filter"]["tag_threshold"]
SCIFI_FANTASY_REGEX = r"(?i)fantasy|sci(ence)?-fi"

print("Selecting works with at least {} ratings".format(RATINGS_THRESHOLD))
raw_works = pl.read_ndjson(snakemake.input["works"]).filter(
    pl.col("ratings_count").cast(int) >= RATINGS_THRESHOLD,
)

print("Collecting tags for selected works")
books = pl.scan_ndjson(snakemake.input["books"], infer_schema_length=1000).filter(
    pl.col("work_id").is_in(raw_works["work_id"].to_list())
)
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
    .collect()
)

works = raw_works.join(work_tags, on="work_id", how="left")

print("Writing filtered works and tags {}".format(snakemake.output["filtered_works"]))
works.write_parquet(snakemake.output["filtered_works"])

best_book_ids = works["best_book_id"].unique().to_list()

print("Combining into one dataset")
best_book_data = books.filter(pl.col("book_id").is_in(best_book_ids)).collect()
complete_works = best_book_data.join(works, on="work_id", how="inner")

print(
    "Writing combined works and tags to {}".format(snakemake.output["combined_works"])
)
complete_works.write_parquet(snakemake.output["combined_works"])
