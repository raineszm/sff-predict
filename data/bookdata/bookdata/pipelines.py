# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: https://docs.scrapy.org/en/latest/topics/item-pipeline.html


# useful for handling different item types with a single interface
from itemadapter import ItemAdapter


class BookdataPipeline:
    def process_item(self, item, spider):
        adapter = ItemAdapter(item)
        nominee = adapter.get("nominee")
        last, first = nominee.split(",")
        adapter["nominee"] = f"{first.strip()} {last.strip()}"
        return item
