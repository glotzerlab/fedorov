# Copyright (c) 2019-2020 The Regents of the University of Michigan
# This file is part of the fedorov project, released under the BSD 3-Clause License.

# Maintainer: Pengji Zhou

# NOTE: this is the code for record that was last used in Jan 2020 to generate all crystal
# data record from Aflow. The use of this code is not required to use this package
# unless user would like to abtain the most updated data again. However this code is not
# maintained and not quaranteed to work if the website it queries is changed afterwards.

import pandas as pd
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from bs4 import BeautifulSoup
import re


# instantiate a chrome options object so you can set the size and headless preference
chrome_options = Options()
chrome_options.add_argument("--headless")
chrome_options.add_argument("--window-size=1920x1080")

driver = webdriver.Chrome(chrome_options=chrome_options)
driver.get('http://aflowlib.org/CrystalDatabase/prototype_index.html')

html = driver.page_source
soup = BeautifulSoup(html)

table = pd.read_html(str(soup.table))

df_source = table[0]
df_source.index.name = 'prototype_index'
# data clean up
df_source.iloc[24, 0] = 'CaCO3'
df_source.iloc[25, 0] = 'FeB'
df_source.iloc[54, 0] = 'alpha-Pa'
df_source.iloc[60, 0] = 'eta-Fe2C'
df_source.iloc[62, 0] = 'alpha-Ga'
df_source.iloc[80, 0] = 'HgCl2'
df_source.iloc[165, 0] = 'SnS'
df_source.iloc[187, 0] = 'alpha-CO'
df_source.iloc[190, 0] = 'beta-Po'
df_source.iloc[204, 0] = 'alpha-Hg'
df_source.iloc[235, 0] = 'alpha-As'
df_source.iloc[253, 0] = 'PbCl2'
df_source.iloc[261, 0] = 'beta-O'
df_source.iloc[377, 0] = 'SeO2'
df_source.iloc[418, 0] = 'SrH2'
df_source.iloc[491, 0] = 'FeAs'
df_source.iloc[565, 0] = 'PbO'

df_source.head(10)

df_source.to_csv('Aflow_raw_data.csv')


df = pd.DataFrame(columns=['id', 'Pearson_symbol', 'space_group_number', 'Wyckoff_site',
                           'lattice_params_list', 'basis_params_list', 'lattice_params_value_list',
                           'basis_params_value_list'])
df.head()


def process(labels):
    output = ''
    n = len(labels)
    i = 0
    while i < n:
        if labels[i].isnumeric():
            if labels[i + 1].isnumeric():
                output += int(labels[i:i + 2]) * labels[i + 2]
                i += 3
            else:
                output += int(labels[i]) * labels[i + 1]
                i += 2
        else:
            output += labels[i]
            i += 1
    return output


def extract_info(string):
    wyckoff_site = []
    labels = string.split('_')[3:]
    for site in labels:
        wyckoff_site.append(process(site))
    return wyckoff_site


for index, row in df_source.iterrows():
    df = df.append({'id': (row['Pearson Symbol'] + '-' + row['Prototype'] + '-'
                    + str(row['Space Group Number'])), 'Pearson_symbol': row['Pearson Symbol'],
                    'space_group_number': row['Space Group Number'],
                    'Wyckoff_site': extract_info(row['AFLOW Prototype'])}, ignore_index=True)


def process_parameters(parameters):
    print(parameters)
    parameters = parameters.split(',')
    params1 = []
    params2 = []
    for para in parameters:
        if para == 'a':
            params1.append('a')
        elif para == 'b/a':
            params1.append('b/a')
        elif para == 'c/a':
            params1.append('c/a')
        elif para == '\\alpha':
            params1.append('alpha')
        elif para == '\\beta':
            params1.append('beta')
        elif para == '\\gamma':
            params1.append('gamma')
        else:
            params2.append(''.join(re.findall(r'[0-9a-z]', para)))
    return params1, params2


i = 0
# get info
while i < 590:
    name = df_source['AFLOW Prototype'][i]
    element = df_source['Prototype'][i]
    print(name, element)
    try:
        driver.get('http://aflowlib.org/CrystalDatabase/%s.html' % name)
        html = driver.page_source
        soup = BeautifulSoup(html, 'html.parser')
        info = str(soup.find_all('script')[-3])
    except BaseException:
        try:
            driver.get('http://aflowlib.org/CrystalDatabase/%s.%s.html' % (name, element))
            html = driver.page_source
            soup = BeautifulSoup(html, 'html.parser')
            info = str(soup.find_all('script')[-3])
        except BaseException:
            driver.get('http://aflowlib.org/CrystalDatabase/%s-%s.html' % (name, element))
            html = driver.page_source
            soup = BeautifulSoup(html, 'html.parser')
            info = str(soup.find_all('script')[-3])

    parameters = re.findall(r'var aflow_params_str = "\s*(.*?)";', info)[0]
    values = re.findall(r'aflow_parameter_values = \s*(.*?);', info)[0]

    params1, params2 = process_parameters(parameters)
    exec('values = {}'.format(values))
    df['lattice_params_list'][i] = params1
    df['basis_params_list'][i] = params2
    df['lattice_params_value_list'][i] = values[:len(params1)]
    df['basis_params_value_list'][i] = values[len(params1):]
    i += 1

# duplication check
print(df['id'].is_unique)
dup = df[df.duplicated(['id'])]
print(dup)

# sanity check
for index, row in df.iterrows():
    number_of_Wyckoff = len(''.join(row['Wyckoff_site']))
    try:
        number_in_para = int(row['basis_params_list'][-1][1:])
    except BaseException:
        number_in_para = number_of_Wyckoff
    try:
        assert(number_of_Wyckoff == number_in_para)
    except BaseException:
        print(number_of_Wyckoff, number_in_para, index)
        print(row)

# mannual corrections
df.loc[39, 'basis_params_list'] = "['x2', 'x3', 'x4', 'x5', 'y5']"
df.loc[73, 'basis_params_list'] = "['x1', 'y1', 'x2', 'y2']"
df.loc[77, 'basis_params_list'] = "['x1', 'y1', 'x2', 'y2']"
df.loc[132, 'basis_params_list'] = "['x1', 'x2']"
df.loc[208, 'basis_params_list'] = """['x2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5', 'x6', 'y6', 'x7',
                                       'y7', 'x8', 'y8', 'x9', 'y9', 'x10', 'y10', 'x11', 'y11',
                                       'x12', 'y12', 'z12', 'x13', 'y13', 'z13', 'x14', 'y14',
                                       'z14', 'x15', 'y15', 'z15']"""
df.loc[434, 'basis_params_list'] = """['x1', 'x2', 'x3', 'y3', 'x4', 'y4', 'x5', 'y5', 'x6', 'y6',
                                       'x7', 'y7', 'x8', 'y8', 'x9', 'y9', 'z9']"""
df.loc[495, 'basis_params_list'] = "['x1', 'x2', 'y2']"
df.index.name = 'prototype_index'

df.to_csv('Aflow_processed_data.csv')
