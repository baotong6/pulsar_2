import numpy as np
import matplotlib.pyplot as plt
from tkinter import _flatten
#cand_id=np.loadtxt('cand_id.txt')
def merge_result_cluster():
    #cand_id=np.loadtxt('cand_id.txt')
    #cand_id=np.linspace(1,518,518)
    cand_id=np.loadtxt('G_src.txt')
    #cand_id = np.linspace(1, 847, 847)
    cand_id=cand_id.astype(int)
    result_srcid = []
    result_runtime = []
    result_Prob = []
    result_wpeak = []
    result_period=[]
    result_mopt = []
    result_wconf_lo = []
    result_wconf_hi = []
    #path='/Users/baotong/Desktop/period/result_NSC_I_all_obs/'
    #path='/Users/baotong/Desktop/period_gc/result_gc_all_obs/'
    path = '/Users/baotong/Desktop/period/result_NSC_merge_IG_all_obs/'
    #path = '/Users/baotong/Desktop/period/qsub/result_G/'
    #path = '/Users/baotong/Desktop/period_LW/result_LW_all_obs/'

    for i in range(16):
        # if i==0:
        #     res=np.loadtxt(path+'result_1h_{0}_{1}.txt'.format(int(cand_id[int(i)*30]),int(cand_id[int((i+1)*30-1)])))
        # if i ==1:
        #     res = np.loadtxt(path + 'result_1h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[int((i + 1) * 30 - 11)])))
        # if i>1:
        #     res = np.loadtxt(path + 'result_1h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30-10]), int(cand_id[int((i + 1) * 30 - 11)])))

        #just for merge 1-3h,history reason
        # if i <3:
        #     res = np.loadtxt(path + 'result_1h_3h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[int((i + 1) * 30 - 1)])))
        # elif i==3:
        #     res = np.loadtxt(path + 'result_1h_3h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[int((i + 1) * 30 - 21)])))
        # elif i>3 and i<15:
        #     res=np.loadtxt(path+'result_1h_3h_{0}_{1}.txt'.format(int(cand_id[int(i)*30-20]),int(cand_id[int((i+1)*30-21)])))
        # elif i==15:
        #     res=np.loadtxt(path + 'result_1h_3h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30-20]),int(cand_id[-1])))
        if i<15:
           res = np.loadtxt(path + 'result_3h_10h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[int((i + 1) * 30 - 1)])))
        else:
           res = np.loadtxt(path + 'result_3h_10h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]),int(cand_id[-1])))

        #for LW
        # if i<28:
        #     res = np.loadtxt(path + 'result_1h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[int((i + 1) * 30 - 1)])))
        # elif i==28:
        #     res = np.loadtxt(path + 'result_1h_{0}_{1}.txt'.format(int(cand_id[int(i) * 30]), int(cand_id[-1])))

        #print(res[:,0])
        result_srcid.append(list(res[:,0]))
        result_runtime.append(list(res[:,1]))
        result_Prob.append(list(res[:,2]))
        result_wpeak.append(list(res[:,3]))
        result_period.append(list(res[:,4]))
        result_mopt.append(list(res[:,5]))
        result_wconf_lo.append(list(res[:,6]))
        result_wconf_hi.append(list(res[:,7]))

    result_srcid=list(_flatten(result_srcid))
    result_runtime=list(_flatten(result_runtime))
    result_Prob=list(_flatten(result_Prob))
    result_wpeak=list(_flatten(result_wpeak))
    result_period=list(_flatten(result_period))
    result_mopt=list(_flatten(result_mopt))
    result_wconf_lo=list(_flatten(result_wconf_lo))
    result_wconf_hi=list(_flatten(result_wconf_hi))
    result = np.column_stack((result_srcid, result_runtime, result_Prob, result_wpeak, result_period, result_mopt,result_wconf_lo, result_wconf_hi))
    print(result_srcid)
    np.savetxt(path+'all_result_3h_10h.txt',result, fmt='%10d %10.2f %10.5f %10.5f %10.5f %10d %10.5f %10.5f')


merge_result_cluster()

