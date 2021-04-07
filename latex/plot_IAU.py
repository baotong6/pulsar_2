import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 1.25*(2/(1+np.e**(-4*x))-1)-x**(0.5)

# x=np.linspace(0,10,1000)
# y=f(x)
# plt.plot(x,y)
# plt.show()

a=5
def dp(r2):
    return 1/np.pi*np.arcsin(1-2*r2/a)

x=np.linspace(0,2.5,1000)
y=dp(x)
plt.xlabel('r2')
plt.scatter(x,y)
plt.show()
