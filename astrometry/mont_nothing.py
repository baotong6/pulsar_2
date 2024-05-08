'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2023-12-05 09:08:52
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2023-12-05 09:14:29
FilePath: /pulsar/astrometry/mont_nothing.py
Description: 

Copyright (c) 2023 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
A1=0;B1=0
for a in range(36):
    for b in range(36):
        for c in range(36):
            d = 35 - (a + b + c)
            if 0 <= d <= 35:
                # print(f"a={a}, b={b}, c={c}, d={d}")
                A1+=1
                if np.max([a,b,c,d])==9:
                    print(a,b,c,d)
                    B1+=1
print(A1,B1)
            

