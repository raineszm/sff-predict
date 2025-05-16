from snakemake.script import snakemake
from isal import igzip
import json
import tqdm

from helpers.files import rawbigcount

book_ids = set()
with open(snakemake.input.books, "r") as f:
    for line in f:
        book_ids.add(json.loads(line)["book_id"])


total_lines = rawbigcount(snakemake.input[0])
with igzip.open(snakemake.input[0], "r") as f, open(snakemake.output[0], "w") as out:
    for line in tqdm.tqdm(
        f, total=total_lines, desc="Filtering reviews", unit_scale=True
    ):
        review = json.loads(line)
        if review["book_id"] in book_ids:
            out.write(json.dumps(review) + "\n")
