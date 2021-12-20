import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import julia

### fuck, still too slow! garbage###

jl = julia.Julia(compiled_modules=False)
jl.include("GL.jl")
# jl.eval("importall julia_util")
# jl.include("test.jl")
# jl.print_test()
path = '/Users/baotong/Desktop/period_Tuc/txt_all_obs_p50/'
dataname='402'
filename = str(dataname) + '.txt'
data_file = path + str(dataname) + '.txt'
epoch_file = path + 'epoch_src_{0}.txt'.format(dataname)
src_evt=np.loadtxt(data_file)
epoch_info=np.loadtxt(epoch_file)
gtis=epoch_info[:,0:2]
Tlist=src_evt[:,0]

starttime = datetime.datetime.now()
w_range=2*np.pi*np.arange(1./10000,1./3000,1.e-6)
w_range=list(w_range);Tlist=list(Tlist)


if __name__=="__main__":
    res=jl.compute_GL(Tlist,w_range)
    endtime = datetime.datetime.now()
    print(endtime-starttime)
    print(res)