import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import functools
import datetime
import GL_algorithm_esass as GL
import funcs_sim_timeseries as funcs_sim
import funcs_timing as funcs_time
import sys

def write_LS_result(path,cts_num,amp_num,sim_N=100):
    path_file='/Users/baotong/eSASS/data/raw_data/47_Tuc/txt/'
    epoch_file = path_file + 'epoch_47Tuc_700011.txt'
    epoch_info=np.loadtxt(epoch_file)
    if epoch_info.ndim==1:epoch_info=np.array([epoch_info])
    FP_all = [];
    out_period_all = [];
    srcID_output = [];
    counts_all = []
    cts_rate_standard=0.01
    amp_standard=1.0
    ifvary_value=10
    varydelta=0
    bin_len=100
    for i in range(sim_N):
        time = funcs_sim.get_epoch_time_series(cts_rate = cts_rate_standard*cts_num, period =7200.,
                                               amp =amp_standard*amp_num,epoch=epoch_info,model = 'sin',varydelta=0)
        lc=funcs_time.get_hist(time, len_bin=bin_len,tstart=0,tstop=0)
        x = lc.time;
        flux = lc.counts
        counts=np.sum(flux)
        print(i)
        T_tot = lc.time[-1] - lc.time[0]
        freq = np.arange(1 / T_tot, 0.5 / bin_len, 1 / (10 * T_tot))
        freq = freq[np.where(freq > 1 / 20000.)]
        [FP, out_period, max_NormLSP] = funcs_time.get_LS(x, flux, freq, outpath=path, outname=str(i+1), save=False,show=False)
        FP_all.append(FP);
        out_period_all.append(out_period);
        srcID_output.append(i+1);
        counts_all.append(counts)

    FP_all = np.array(FP_all);
    out_period_all = np.array(out_period_all);
    srcID_output = np.array(srcID_output);
    counts_all = np.array(counts_all)
    LS_info = np.column_stack((srcID_output, FP_all, out_period_all, counts_all))
    path_out=path+'result_{0}_{1}/'.format(cts_num,amp_num)
    np.savetxt(path_out + 'sim_LS_bin{0}_VD{1}_N{2}.txt'.format(bin_len,varydelta,sim_N), LS_info, fmt='%10d %10.5f %10.5f %10d')

if __name__=='__main__':
    path='/Users/baotong/eSASS/data/raw_data/47_Tuc/simulation/const_delta_25k_LS/period_2h/'
    cts_range = [0.5,0.7,1.0,1.5,2.0]
    # amp_range = [0.2,0.3,0.4,0.5,0.6,0.7,0.8]
    amp_range = [0.7, 0.8]
    for x in range(len(cts_range)):
        cts_num=cts_range[x]
        for y in range(len(amp_range)):
            amp_num=amp_range[y]
            write_LS_result(path, cts_num=cts_num, amp_num=amp_num, sim_N=100)