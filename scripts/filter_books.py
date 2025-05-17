from snakemake.script import snakemake
from isal import igzip
import json
import tqdm

from helpers.files import rawbigcount

SCIFI_FANTASY_SHELVES = snakemake.config["filter"]["genres"]
RATINGS_THRESHOLD = snakemake.config["filter"]["ratings_threshold"]
TAG_THRESHOLD = RATINGS_THRESHOLD * snakemake.config["filter"]["tag_threshold"]


def shelf_matches(shelf):
    return (
        any(tag in shelf["name"] for tag in SCIFI_FANTASY_SHELVES)
        and int(shelf["count"]) >= TAG_THRESHOLD
    )


def is_scifi_fantasy(book):
    if not book["ratings_count"] or int(book["ratings_count"]) < RATINGS_THRESHOLD:
        return False
    return any(shelf_matches(shelf) for shelf in book["popular_shelves"])


total_lines = rawbigcount(snakemake.input[0])
with (
    igzip.open(snakemake.input[0], "r") as f,
    igzip.open(snakemake.output[0], "w") as out,
):
    for line in tqdm.tqdm(
        f, total=total_lines, desc="Filtering books", unit_scale=True
    ):
        book = json.loads(line)
        if is_scifi_fantasy(book):
            out.write(line)
