# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that was last used in June 2020 to generate all PointGroup info
# The use of this code is not required to use this package
# unless user would like to abtain the most updated data again. However this code is not
# maintained and not quaranteed to work if the website it queries is changed afterwards.


from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import re
import numpy as np
import pickle
import rowan
import json


# instantiate a chrome options object so you can set the size and headless preference
chrome_options = Options()
# chrome_options.add_argument("--headless")
chrome_options.add_argument("--window-size=1920x1080")

driver = webdriver.Chrome(chrome_options=chrome_options)

point_group_rotation_matrix_dict = {}
for num in range(1, 33):
    driver.get('https://www.cryst.ehu.es/cgi-bin/cryst/programs/'
               'nph-point_genpos?num={}'.format(num))

    html = driver.page_source
    soup = BeautifulSoup(html, 'html.parser')

    table = soup.find_all('center')[1]
    table = table.find('tbody')

    n = len(table.find_all('tr', recursive=False)[2:])
    rotations = np.zeros((n, 3, 3))
    i = 0
    for row in table.find_all('tr', recursive=False)[2:]:
        table = row.find('table')
        value_list = re.findall(r'[-+0-9./]+', table.find('pre').string)
        value_list = [eval(item) for item in value_list]
        rotations[i, 0, :] = value_list[0:3]
        rotations[i, 1, :] = value_list[3:6]
        rotations[i, 2, :] = value_list[6:9]
        i += 1
    point_group_rotation_matrix_dict[num] = {'rotations': rotations}


with open('point_group_rotation_matrix_dict.pickle', 'wb') as f:
    pickle.dump(point_group_rotation_matrix_dict, f)

point_group_list = ['1', '-1', '2', 'm', '2/m', '222', 'mm2', 'mmm',
                    '4', '-4', '4/m', '422', '4mm', '-42m', '4/mmm', '3',
                    '-3', '32', '3m', '-3m', '6', '-6', '6/m', '622',
                    '6mm', '-6m2', '6/mmm', '23', 'm-3', '432', '-43m', 'm-3m']

num = range(1, 33)
point_group_name_dict = dict(zip(num, point_group_list))

with open('point_group_name_mapping.json', 'w') as f:
    json.dump(point_group_name_dict, f)

point_group_quat_dict = {}
for key, item in point_group_rotation_matrix_dict.items():
    quats = []
    n = item['rotations'].shape[0]
    for i in range(0, n):
        qtemp = rowan.from_matrix(item['rotations'][i, :, :], require_orthogonal=False)
        quats.append(qtemp.tolist())

    point_group_quat_dict[key] = quats

with open('point_group_quat_dict.json', 'w') as f:
    json.dump(point_group_quat_dict, f)
