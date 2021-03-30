# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause
# License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that was last used in May 2020 to generate
# all PlaneGroup info The use of this code is not required to use this package
# unless user would like to abtain the most updated data again. However this
# code is not maintained and not quaranteed to work if the website it queries is
# changed afterwards.

import pickle
import re

import numpy as np
from bs4 import BeautifulSoup
from selenium import webdriver
from selenium.webdriver.chrome.options import Options

# instantiate a chrome options object so you can set the size and headless
# preference
chrome_options = Options()
# chrome_options.add_argument("--headless")
chrome_options.add_argument("--window-size=1920x1080")

driver = webdriver.Chrome(chrome_options=chrome_options)

planegroup_dict = {}
for num in range(1, 18):
    driver.get("https://www.cryst.ehu.es/plane/get_plane_gen.html")

    inputElement = driver.find_element_by_name("gnum")
    inputElement.send_keys(num)
    driver.find_element_by_name("list").click()

    html = driver.page_source
    soup = BeautifulSoup(html, "html.parser")

    table = soup.find("center")
    table = table.find("tbody")

    n = len(table.find_all("tr", recursive=False)[2:])
    rotations = np.zeros((n, 2, 2))
    translations = np.zeros((n, 2))
    i = 0
    for row in table.find_all("tr", recursive=False)[2:]:
        value_list = re.findall(r"[-+0-9./]+", row.find_all("td")[0]["id"])
        value_list = [eval(item) for item in value_list]
        rotations[i, 0, :] = value_list[0:2]
        rotations[i, 1, :] = value_list[3:5]
        translations[i, 0] = value_list[2]
        translations[i, 1] = value_list[5]
        i += 1
    planegroup_dict[num] = {
        "rotations": rotations,
        "translations": translations,
    }


with open("plane_group_info.pickle", "wb") as f:
    pickle.dump(planegroup_dict, f)
