# -*- coding: utf-8 -*-

# Define here the models for your scraped items
#
# See documentation in:
# https://doc.scrapy.org/en/latest/topics/items.html

import scrapy
from scrapy import Field


class SpiderItem(scrapy.Item):
    # define the fields for your item here like:
    name_EN = Field()
    name_CH = Field()
    CAS = Field()
    MDL = Field()
    prod_num = Field()
    Smiles = Field()
    structure = Field()
    sell_info = Field()
    purity = Field()
    picture_url = Field()
