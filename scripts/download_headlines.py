from snakemake.script import snakemake
from dotenv import load_dotenv
from utils.world_state import NYTimesArchiveAPI
import itertools
from tqdm.auto import tqdm
import json
from httpx import HTTPStatusError

from loguru import logger

logger.remove()
logger.add("logs/headlines.log", level="DEBUG")

load_dotenv()


nyt = NYTimesArchiveAPI.from_envfile()

with open(snakemake.output.headlines, "w") as f:
    for year, month in tqdm(
        itertools.product(
            range(1959, 2025),
            range(1, 13),
        ),
        total=len(range(1959, 2025)) * len(range(1, 13)),
        desc="Downloading headlines",
    ):
        try:
            for headline in nyt.get_headlines(year, month):
                json.dump(headline, f)
                f.write("\n")
        except HTTPStatusError as e:
            if e.response.status_code == 429:
                msg = e.response.json()["fault"]["faultstring"]
                logger.error(
                    f"Hit rate limit downloading headlines for {year}/{month}: {msg}"
                )
                print("Hit rate limit, try again tomorrow.")
            break
