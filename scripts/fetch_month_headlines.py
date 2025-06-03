from snakemake.script import snakemake
from dotenv import load_dotenv
from utils.world_state import NYTimesArchiveAPI
import json
from httpx import HTTPStatusError

from loguru import logger

logger.remove()
logger.add("logs/headlines.log", level="DEBUG")

load_dotenv()


nyt = NYTimesArchiveAPI.from_envfile()
try:
    with open(snakemake.output.month_headlines, "w") as f:
        for headline in nyt.get_headlines(
            snakemake.params.year, snakemake.params.month
        ):
            json.dump(headline, f)
            f.write("\n")
except HTTPStatusError as e:
    if e.response.status_code == 429:
        msg = e.response.json()["fault"]["faultstring"]
        logger.error(
            f"Hit rate limit downloading headlines for {snakemake.params.year}/{snakemake.params.month}: {msg}"
        )
        print("Hit rate limit, try again tomorrow.")
        exit(1)
