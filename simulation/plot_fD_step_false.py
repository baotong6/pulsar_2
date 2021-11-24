import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import sys
import os
import string
import pandas as pd
###----------------read table---------------###
path_table = '/Users/baotong/Desktop/period/table/'
result_NSC_IG=pd.read_excel(path_table + 'final_all_del.csv', 'result_NSC_IG')
ID_NSC_IG=result_NSC_IG['seq']
P_NSC_IG=result_NSC_IG['P']
counts_NSC_IG=result_NSC_IG['counts']
###----------------read table---------------###

def get_cts_obs():
    path_txt='/Users/baotong/Desktop/period/txt_all_obs_I/'
    ID_num=np.linspace(1,3619,3619)#源的数目
    # ID_num=np.linspace(1,3619,3619)
    # [87.48, 180.45, 274.98, 369.05, 461.18, 555.41, 649.47, 744.78, 836.88]
    cts = []
    for id in ID_num:
        file_name = str(int(id)) + '.txt'
        file_info = np.loadtxt(path_txt + file_name)
        cts.append(len(file_info[:, 0]))
    return cts

path='/Users/baotong/Desktop/period/simulation/simulation_NSC_I_step/'
#path='/Users/baotong/Desktop/period/simulation/simulation_NSC_I_step_VI2/'
threshold=0.9
cts_range=[1.0,2.0,3.0,4.0,5.0,6.0,8.0,10.0,12.0,13.0,15.0]
amp=0.5  ## useless
period_range=['20k','30k','40k','50k','60k','70k']
sim_N=100
sim_id_range=np.linspace(1,sim_N,sim_N)
sim_id_range=sim_id_range.astype(int)
font1 = {'family': 'Normal',
         'weight': 'normal',
         'size': 18, }
##src 2148-false-detectionr-rate=0.23

def plot_CRbased_contour():
    ##默认每个CR对应相同的cts,取值为cts_mean
    detect_rate=[]
    counts_mean=[]
    for cts_rate in cts_range:
        res_detect = []
        counts=0
        for period in period_range:
            detect=0
            period_get=[]
            for i in sim_id_range:
                temp_info=np.loadtxt(path+'result_{0}_{1}_{2}/result_sim_{3}.txt'.format(str(cts_rate),str(amp),period,str(i)))
                # if temp_info[2]>threshold:
                if temp_info[2]>threshold :
                        # or 1.98*period_real<temp_info[4] <2.01*period_real\
                        # or 2.97*period_real<temp_info[4] <3.03*period_real :
                    detect+=1
                counts+=temp_info[-1]

            res_detect.append(detect/sim_N)
        counts/=sim_N*len(period_range)
        counts_mean.append(counts)
        detect_rate.append(res_detect)

    lab = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    X = [1.5e4, 2.5e4, 3.5e4, 4.5e4, 5.5e4, 6.5e4]
    Y = counts_mean
    Z = detect_rate
    plt.contourf(X, Y, Z, levels=lab)
    plt.xlabel('Period')
    plt.ylabel('Counts')
    # plt.tick_params(labelsize=18)
    cbar = plt.colorbar()
    cbar.set_label('False detection rate')
    cbar.set_ticks(np.linspace(0, 1.0, 11))
    # cbar.set_ticklabels(('0.0','0.1','0.2','0.4','0.3','0.5','0.6','0.7','0.8','0.9','1.0'))
    plt.show()
#plot_CRbased_contour()
def plot_CRbased_contour_2():
    ##不同CR的结果放在一起取cts bin
    X = [1.5e4, 2.5e4, 3.5e4, 4.5e4, 5.5e4, 6.5e4]
    cts_bin = np.linspace(200, 5000, 11)
    cts_bin=np.concatenate((cts_bin,[15000]))
    false_out=[]
    for k in range(len(period_range)):
        period_lab=period_range[k]
        detect_prob=[]
        counts=[]
        false_index=[]
        for cts_rate in cts_range:
            period_get=[]
            for i in sim_id_range:
                temp_info=np.loadtxt(path+'result_{0}_{1}_{2}/result_sim_{3}.txt'.format(str(cts_rate),str(amp),period_lab,str(i)))
                if temp_info[2] > threshold:
                    detect_prob.append(1)
                else:
                    detect_prob.append(0)
                counts.append(temp_info[-1])
        # plt.hist(counts,bins=10)
        false_all=np.column_stack((counts,detect_prob))
        false_all = false_all[np.lexsort(false_all[:, ::-1].T)]
        index=plt.hist(false_all[:,0],bins=cts_bin)[0]
        plt.close()
        false_index.append(np.mean(false_all[:,1][0:int(index[0])]))
        for j in range(len(index)-1):
            false_index.append(np.mean(false_all[:,1][int(np.sum(index[0:j])):int(np.sum(index[0:j])+index[j+1])]))

        false_out.append(false_index)

    false_out=np.array(false_out)
    false_out=false_out.T

    lab = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    X = [1.5e4, 2.5e4, 3.5e4, 4.5e4, 5.5e4, 6.5e4]
    Y = cts_bin[0:-1]
    Z = false_out
    plt.figure(1)
    plt.contourf(X, Y+cts_bin[1]-cts_bin[0], Z, levels=lab)
    plt.xlabel('Period')
    plt.ylabel('Counts')
    # plt.tick_params(labelsize=18)
    cbar = plt.colorbar()
    cbar.set_label('False detection rate')
    cbar.set_ticks(np.linspace(0, 1.0, 11))
    # cbar.set_ticklabels(('0.0','0.1','0.2','0.4','0.3','0.5','0.6','0.7','0.8','0.9','1.0'))
    plt.show()

    Z = np.array(Z)
    ctsperbin = plt.hist(get_cts_obs(), cts_bin)[0]
    plt.close()
    false_result = ctsperbin * Z[:, 3]
    print(ctsperbin)
    print(Z[:, 3])
    print(false_result)

    obs_longpCV_cts=counts_NSC_IG[43:]
    obs_longpCV_P= P_NSC_IG[43:]
    # plt.hist(obs_longpCV_cts,bins=cts_bin,color='red')
    plt.figure(2)
    plt.bar(np.array(cts_bin[0:-1])+0.5*(cts_bin[1]-cts_bin[0]),false_result,width=cts_bin[1]-cts_bin[0],color="w",edgecolor="k",linewidth=5.)
    plt.hist(obs_longpCV_cts, bins=cts_bin, color='red')
    plt.show()
#plot_CRbased_contour_2()
def plot_pCV_fD():
    sim_N=50
    path='/Users/baotong/Desktop/period_Tuc/simulation/'
    #path = '/Users/baotong/Desktop/period_terzan5/simulation/'
    # ID_50k=['1674','442','1525','2344','2238','1538','1677','3564',
    #     '2199','2338','3107','1634','214','3401','1628','790','1514',
    #     '1941','1769','3596','2574','973','2187','307']
    # ID_40k=['1529','3483','2157','6']
    # ID_30k = ['1514', '2525']
    # ID_20k=['147','1084','1538','2730','1133','3120','1853','1487','2422','1219','3357','2508']
    #ID_50k=['28','41','52','225']
    # ID_40k = ['49']
    threshold = 0.99
    #ID_10k=['2020','2313','2355','687','2478']
    #ID_50k=['84','10','83','106','117','118','241','104','116','1','82','74']
    # ID_30k=[15,57,121,81,24,10]
    ID_50k=[232,258,283,345,223,261,350,567]
    period_det_50k=[37268.,32407.55744,44622.93619,48837.66361,46151.00609,23917.72303,22026.91689,38008.36184]
    # ID_50k=[84,14,37,294,85,128,78,372,374,66,292,55]
    # period_det_50k=[10431.44455,15797.78831,22805.0171,24826.21649,27987.68542,29877.50224,31486.1461,43591.97908,44404.97336,44483.98577,48473.09743,49382.71605]
    # period_det_30k=[28650.01146,11015.27819,10734.68161,19549.95992,11434.452,27638.81595]
    ID=ID_50k;period_det=period_det_50k

    detect_rate=np.zeros(len(ID))
    for k in range(len(ID)):
        for i in range(sim_N):
            if os.path.exists(path+'result_all_{0}_50k/result_sim_{1}.txt'.format(ID[k],i+1)):
                res=np.loadtxt(path+'result_all_{0}_50k/result_sim_{1}.txt'.format(ID[k],i+1))
                if res[2] > threshold and (np.abs(res[4]-period_det[k]))>0.0001*period_det[k]:
                    detect_rate[k]+=1
            else:
                print(ID[k])
    detect_rate/=sim_N+0.
    for i in range(len(ID)):
        print("ID={0},FD={1}".format(ID[i],detect_rate[i]))

plot_pCV_fD()

# print(counts_mean)
# for i in range(len(detect_rate)):
#     x=[5e4]
#     y=detect_rate[i]
#     plt.plot(x,y,marker='v')
# plt.xlabel('period')
# plt.ylabel('false detection rate')
# #plt.legend(['C~300','C~1200','C~1500','C~2500','C~3000','C~3600','C~4500'])
# plt.show()
# print(detect_rate)
