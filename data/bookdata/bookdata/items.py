# Define here the models for your scraped items
#
# See documentation in:
# https://docs.scrapy.org/en/latest/topics/items.html

import scrapy


class BookDataItem(scrapy.Item):
    title = scrapy.Field()  # Book title
    year = scrapy.Field()  # Year of the award
    result = scrapy.Field()  # Award result (winner, finalist, etc.)
    award = scrapy.Field()  # Name of the award
    nominee = scrapy.Field()  # Author(s) nominated for the award
