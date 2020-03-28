# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that was last used in Jan 2020 to generate all crystal
# data record from Aflow. The use of this code is not required to use this package
# unless user would like to abtain the most updated data again. However this code is not
# maintained and not quaranteed to work if the website it queries is changed afterwards.

from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import json


# instantiate a chrome options object so you can set the size and headless preference
chrome_options = Options()
# chrome_options.add_argument("--headless")
chrome_options.add_argument("--window-size=1920x1080")

driver = webdriver.Chrome(chrome_options=chrome_options)

for num in range(1, 231):
    driver.get('https://www.cryst.ehu.es/cryst/get_wp.html')

    inputElement = driver.find_element_by_name("gnum")
    inputElement.send_keys(num)
    driver.find_element_by_name('standard').click()

    html = driver.page_source
    soup = BeautifulSoup(html, 'html.parser')

    table = soup.find('center')
    table = table.find('tbody')

    Wyckoff_positions = {}
    i = 0
    for row in table.find_all('tr', recursive=False)[2:]:
        Wyckoff_site = row.find_all('td')[1].text
        position = row.find_all('td')[3].text.split(' ')[0].strip("()").split(',')
        Wyckoff_positions[Wyckoff_site] = position
        i += 1
    with open('space_group_{}_Wyckoff_site_data.json'.format(num), 'w') as f:
        json.dump(Wyckoff_positions, f)
