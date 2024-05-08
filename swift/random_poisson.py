'''
Author: baotong && baotong@smail.nju.edu.cn
Date: 2024-02-01 16:45:26
LastEditors: baotong && baotong@smail.nju.edu.cn
LastEditTime: 2024-02-01 16:46:43
FilePath: /pulsar/swift/random_poisson.py
Description: 

Copyright (c) 2024 by baotong, All Rights Reserved. 
'''
import numpy as np
import matplotlib.pyplot as plt

a=np.random.poisson(10,1000)
b=np.random.poisson(100,1000)
c=np.random.poisson(1000,1000)
plt.hist(a/10,bins=50,histtype='step')
plt.hist(b/100,bins=50,histtype='step')
plt.hist(c/1000,bins=50,histtype='step')
plt.show()