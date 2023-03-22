import numpy as np
import matplotlib.pyplot as plt

path='/Users/baotong/Downloads/hzq2/'
time=np.load(path+'time.npy')
flux=np.load(path+'flux.npy')
index_high=np.where((time>1.92)&(time<2.0))

time_high=time[index_high]
flux_high=flux[index_high]
print(time[3]-time[2])
# np.save(path+'time_h.npy',time_high)
# np.save(path+'flux_h.npy',flux_high)
plt.plot(time_high,flux_high)
plt.plot(time_high,2e3*np.sin(time_high*2*np.pi/0.01282-0.01)+3.7e4)
plt.show()