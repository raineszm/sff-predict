import scrapy
import re
from w3lib.html import remove_tags
from bookdata.items import BookDataItem

NOVEL_RE = re.compile(r"\x97\s*\bnovel\b\s*\x97", re.IGNORECASE)

AWARDS = [
    "Hugo Awards",
    "Nebula Awards",
    "World Fantasy Awards",
    "John W Campbell Memorial Award",
    "Dragon Awards",
    "International Fantasy Awards",
    "Philip K Dick Award",
]

MULTI_CATEGORY_AWARDS = [
    "Hugo Awards",
    "Nebula Awards",
    "World Fantasy Awards",
    "International Fantasy Awards",
    "Dragon Awards",
]


def get_url_from_award(award):
    return f"http://www.sfadb.com/{award.replace(' ', '_')}_All_Nominees"


def get_award_from_url(url):
    return url.split("/")[-1].replace("_All_Nominees", "").replace("_", " ").title()


def is_novel_category(award, title_mid):
    return (award not in MULTI_CATEGORY_AWARDS) or NOVEL_RE.search(title_mid)


def get_title_and_result(title_mid):
    components = [c.strip() for c in remove_tags(title_mid).split("\x97")]
    if len(components) < 2:
        return False, None, None
    return True, components[0], components[1]


class SfadbSpider(scrapy.Spider):
    name = "sfadb"
    allowed_domains = ["sfadb.com"]
    start_urls = [get_url_from_award(award) for award in AWARDS]

    def parse(self, response):
        award = get_award_from_url(response.url)

        for nomineeblock in response.css(".nomineeblock"):
            nominee = nomineeblock.css(".nominee a::text").getall()
            for year, title_mid in zip(
                nomineeblock.css(".dateleftindent a::text").getall(),
                nomineeblock.css(".titlemid").getall(),
            ):
                if not is_novel_category(award, title_mid):
                    continue

                is_valid, title, result = get_title_and_result(title_mid)

                if not is_valid:
                    continue

                yield BookDataItem(
                    title=title,
                    year=year,
                    result=result,
                    award=award,
                    nominee=" ".join(nominee).strip(),
                )
