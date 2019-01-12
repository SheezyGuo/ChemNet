# -*- coding: utf-8 -*-
import re
import urllib

import scrapy
from scrapy_splash import SplashRequest

from Spider.items import SpiderItem


class CheminfoSpider(scrapy.Spider):
    name = 'ChemInfo'
    allowed_domains = ['jkchemical.com']
    site_root = "www.jkchemical.com"
    start_urls = []

    lua_script = \
        '''
        function
            main(splash, args)
            splash: set_viewport_size(1028, 10000)
            splash: go(args.url)
            local
            scroll_to = splash:jsfunc("window.scrollTo")
            scroll_to(0, 4000)
            splash: wait(5)
            return {
                html = splash:html()
            }
        end
        '''

    def __init__(self):
        search_url_pattern = 'http://www.jkchemical.com/CH/products/search/fulltextsearch/{:s}.html'
        key_words = ['醛']
        print(search_url_pattern)
        self.start_urls = [search_url_pattern.format(urllib.parse.quote(key_word, safe="/")) for key_word in key_words]
        pass

    def parse(self, response):
        trows = response.xpath(r'//*[@id="ctl00_ContentPlaceHolder1_up_Product_GridView1"]/tbody/tr')
        for row in trows:
            item = SpiderItem()
            item["name_EN"] = row.xpath("./td/div/div[1]/a/text()").extract()[0]
            item["name_EN"] = self.remove_space(item["name_EN"])
            if re.search(r"\,[0-9\.]{1,5}%", item["name_EN"]):
                item["name_EN"] = re.sub(r"\,[0-9\.]{1,5}%", "", item["name_EN"])

            item["name_CH"] = row.xpath("./td/div/div[1]/a/span[1]/text()").extract()[0]
            item["name_CH"] = self.remove_space(item["name_CH"])

            item["purity"] = row.xpath("./td/div/div[2]/div[2]/ul/li[1]/text()").extract()[0]
            item["purity"] = self.remove_space(item["purity"]).split("：")[1]

            item["CAS"] = row.xpath("./td/div/div[2]/div[2]/ul/li[2]/a/text()").extract()[0]
            item["CAS"] = self.remove_space(item["CAS"])

            item["MDL"] = row.xpath("./td/div/div[2]/div[2]/ul/li[3]/a/text()").extract()[0]
            item["MDL"] = self.remove_space(item["MDL"])

            item["prod_num"] = row.xpath("./td/div/div[2]/div[2]/ul/li[4]/text()").extract()[0]
            item["prod_num"] = self.remove_space(item["prod_num"]).split("：")[1]

            item["structure"] = row.xpath("./td/div/div[2]/div[2]/ul/li[5]").xpath("string(.)").extract()[0]
            item["structure"] = self.remove_space(item["structure"]).split("：")[1]

            item['picture_url'] = row.xpath("./td/div/div[2]/div[1]/img/@src").extract()[0]
            item["picture_url"] = "{:s}{:s}".format(self.site_root, self.remove_space(item["picture_url"]))

            sell_table = row.xpath('.//div[@class="PRODUCT_box"]/div[contains(@id,"up_PackAge_dv_PackAge")]/div/table/tbody')
            sell_info_list = []
            for info_table_row in sell_table.xpath("./tr")[1:]:
                tds = info_table_row.xpath("./td")
                sell_info = dict()
                sell_info["package"] = tds[0].xpath("string(.)").extract()[0]
                sell_info["unit_price"] = tds[1].xpath("string(.)").extract()[0]
                sell_info["arrive_time"] = tds[4].xpath("string(.)").extract()[0]
                if sell_info["arrive_time"] == '现货\xa0详情':
                    sell_info["arrive_time"] = "现货"
                sell_info_list.append(sell_info)
            item["sell_info"] = sell_info_list
            yield item
        next_url = response.xpath('//*[@id="ctl00_ContentPlaceHolder1_Next"]/@href').extract()[0]
        if next_url:
            next_url = "http://{:s}{:s}".format(self.site_root, next_url)
            next_url = urllib.parse.quote(next_url, safe=":/")
            yield SplashRequest(url=next_url, callback=self.parse, endpoint='execute', args={"lua_source": self.lua_script, "wait": 2.0, "url": next_url}, encoding='utf-8')

    def start_requests(self):

        for url in self.start_urls:
            yield SplashRequest(url=url, callback=self.parse, endpoint='execute', args={"lua_source": self.lua_script, "wait": 5.0}, encoding='utf-8')

    def remove_space(self, str):
        return re.sub("\s", "", str)
