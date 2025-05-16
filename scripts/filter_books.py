from snakemake.script import snakemake
from isal import igzip
import json
import tqdm

from itertools import takewhile, repeat


def rawbigcount(filename):
    f = open(filename, "rb")
    bufgen = takewhile(lambda x: x, (f.raw.read(1024 * 1024) for _ in repeat(None)))
    return sum(buf.count(b"\n") for buf in bufgen if buf)


SCIFI_FANTASY_SHELVES = snakemake.config["filter"]["genres"]
RATINGS_THRESHOLD = snakemake.config["filter"]["ratings_threshold"]
TAG_THRESHOLD = RATINGS_THRESHOLD * snakemake.config["filter"]["tag_threshold"]


def shelf_matches(shelf):
    return (
        shelf["name"] in SCIFI_FANTASY_SHELVES and int(shelf["count"]) >= TAG_THRESHOLD
    )


def is_scifi_fantasy(book):
    if not book["ratings_count"] or int(book["ratings_count"]) < RATINGS_THRESHOLD:
        return False
    return any(shelf_matches(shelf) for shelf in book["popular_shelves"])


total_lines = rawbigcount(snakemake.input[0])
with igzip.open(snakemake.input[0], "r") as f:
    with open(snakemake.output[0], "w") as out:
        for line in tqdm.tqdm(f, total=total_lines):
            book = json.loads(line)
            if is_scifi_fantasy(book):
                out.write(json.dumps(book) + "\n")
