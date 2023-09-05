#!/bin/bash
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 17 8:58:00 2023
@author: baotong
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import re

def extract_numbers(string_list):
    numbers = []
    pattern = r'([\d.]+)(?=\()'  # 正则表达式模式，匹配数字和左括号之间的内容
    for string in string_list:
        match = re.search(pattern, string)  # 使用re.search()查找第一个匹配的结果
        if match:
            number = match.group()  # 提取括号前的数字
            numbers.append(float(number))  # 将提取的数字添加到列表中
        else:
            try:
                numbers.append(float(string))
            except ValueError:
                numbers.append(None)

    return numbers


# strings = ['123(456)', '789(012)', '345(678)']
# result = extract_numbers(strings)
# print(result)
def readPSDSS():
    file_path = '/Users/baotong/Desktop/period_terzan5/tabledata.csv'
    data = pd.read_csv(file_path)
    Pinday=data['Period_in_d'].head(326).tolist()
    period=extract_numbers(Pinday)
    period=np.array(period)*24

    return period

if __name__=='__main__':
    P_SDSS=readPSDSS()