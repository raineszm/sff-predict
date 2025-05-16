import scrapy
import re
import w3lib.html

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


class SfadbSpider(scrapy.Spider):
    name = "sfadb"
    allowed_domains = ["sfadb.com"]
    start_urls = [get_url_from_award(award) for award in AWARDS]

    def parse(self, response):
        award = get_award_from_url(response.url)

        for nomineeblock in response.css(".nomineeblock"):
            for year, nominee, title_mid in zip(
                nomineeblock.css(".dateleftindent a::text").getall(),
                nomineeblock.css(".nominee::text").getall(),
                nomineeblock.css(".titlemid").getall(),
            ):
                if award in MULTI_CATEGORY_AWARDS and not NOVEL_RE.search(title_mid):
                    continue

                components = [
                    c.strip() for c in w3lib.html.remove_tags(title_mid).split("\x97")
                ]

                # skip year with no award
                if len(components) < 2:
                    continue

                if award in MULTI_CATEGORY_AWARDS:
                    title, _, result = components
                else:
                    title, result = components

                yield {
                    "title": title,
                    "year": year,
                    "result": result,
                    "award": award,
                }
