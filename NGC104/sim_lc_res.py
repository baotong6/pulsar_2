import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange
from scipy import interpolate
from astropy.io import fits
from scipy import optimize as op
import hawkeye as hawk


def get_info(path,period,sim_N,cts_range,amp_range,threshold=0.9):
    period_real=period
    sim_id_range=np.linspace(1,sim_N,sim_N)
    sim_id_range=sim_id_range.astype(int)
    detect_rate=[]
    mean=[]
    var=[]
    std=[]
    for cts_rate in cts_range:
        res_detect = []
        res_mean = []
        res_var = []
        res_std = []
        for amp in amp_range:
            detect=0
            period_get=[]
            for i in sim_id_range:
                temp_info=np.loadtxt(path+'result_{0}_{1}/result_sim_{2}.txt'.format(str(cts_rate),str(amp),str(i)))

                if temp_info[2]>threshold and 0.9*period_real<temp_info[3]<1.1*period_real:
                    detect+=1
                    period_get.append(temp_info[4])
            period_get_array=np.array(period_get)

            if len(period_get_array)==0:
                period_mean=0
                period_var=0
                period_std=0
            else:
                period_mean=np.mean(period_get_array)  ##均值
                period_var=np.var(period_get_array)    ##方差
                period_std=np.std(period_get_array,ddof=1)  ##标准差

            res_detect.append(detect/sim_N)
            res_mean.append(period_mean)
            res_var.append(period_var)
            res_std.append(period_std)

        detect_rate.append(res_detect)
        mean.append(res_mean)
        var.append(res_var)
        std.append(res_std)

    return [detect_rate,mean,var,std]

def make_plot(path,period_real,cts_range,amp_range,save=0,show=1,figurepath=None):
    [detect_rate,mean,var,std]=get_info(path=path,period=period_real,
                                        sim_N=100,cts_range=cts_range,amp_range=amp_range)
    print(detect_rate)
    plt.figure(1)
    # plt.title('detection')
    for i in range(len(cts_range)):
        plt.plot(amp_range,detect_rate[i],marker='v',label=f'Counts={int(cts_range[i]*100)}')
    plt.legend()
    #plt.legend(['cr=1','cr=2','cr=3','cr=4','cr=5'])
    #plt.ylim(ymax=1.01)
    plt.xlabel('Amplitude',hawk.font1)
    plt.ylabel('Detection rate',hawk.font1)
    plt.tick_params(labelsize=18)
    # plt.text(0.7,0.4,'P=5540s')
    #plt.savefig('detection.eps')
    yticks=np.linspace(0,1,11)
    plt.yticks(yticks)
    if save:
        plt.savefig(figurepath + f'DR_sin_{period_real}.pdf', bbox_inches='tight', pad_inches=0.05)
    if show:
        plt.show()

def est_detnum(path_in,cts_range,amp_range):
    path='/Users/baotong/Desktop/period_Tuc/'
    src_info=np.loadtxt(path+'src_info.txt')
    bright_in_src_index=np.where((src_info[:,3]>99)&(src_info[:,3]<300)&(src_info[:,6]<3.17*60))[0]
    print(src_info[:,3][bright_in_src_index])
    [detect_rate,mean,var,std]=get_info(path=path_in,period=period_real,
                                        sim_N=100,cts_range=cts_range,amp_range=amp_range)
    mean_DR=[np.mean(detect_rate[i]) for i in range(len(detect_rate))]
    cts_range=np.array(cts_range)
    f=interpolate.interp1d(cts_range*100,mean_DR,kind='quadratic')
    print(np.sum(f(src_info[:,3][bright_in_src_index])))

def est_bkg():
    path = '/Users/baotong/Desktop/period_Tuc/'
    src_info=np.loadtxt(path+'src_info.txt')
    bins=np.logspace(0,np.log10(300),30)
    bkg_scal_cts=src_info[:,4]
    dist_all=src_info[:,7]
    bright_in_src_index = np.where((src_info[:, 3] > 99) & (src_info[:, 3] < 300) & (src_info[:, 7] < 3.17 * 60))[0]
    plt.hist(bkg_scal_cts[bright_in_src_index],bins=bins,histtype='step')
    plt.semilogx()
    plt.show()

def plot_cts_L():
    path = '/Users/baotong/Desktop/period_Tuc/'
    catalog='cheng2019_Tuc.fit'
    src_info=np.loadtxt(path+'src_info.txt')
    catalog_data=fits.open(path+catalog)
    luminosity=catalog_data[1].data['L0_5-8']
    counts=src_info[:,3];bkg_cts=src_info[:,4]
    net_cts=counts-bkg_cts;net_cts[np.where(net_cts<=0)]=1
    bright_index=np.where(net_cts>50)[0]
    def linear(x, A, B):
        return A * x + B

    def curvefit_linear(x, y):
        # param_bounds = ((period * 0.95, 0, 8, 8), (period * 1.05, 0.1, 10, 10))
        popt, pcov = op.curve_fit(linear, x, y)
        perr = np.sqrt(np.diag(pcov))
        return (popt, perr)
    (popt,perr)=curvefit_linear(np.log10(net_cts[bright_index]),np.log10(luminosity[bright_index]))
    print(popt)
    print(linear(np.log10(400),popt[0],popt[1]))
    y1 = linear(np.log10(net_cts), popt[0], popt[1])
    plt.scatter(net_cts,luminosity)
    plt.plot(net_cts[bright_index],10**y1[bright_index],'--','r')
    plt.loglog()
    plt.ylim(1e29,1e35)
    # plt.xlim(-10,10000)
    plt.show()

if __name__ == '__main__':
    period_real=5400;
    path=f'/Users/baotong/Desktop/period_Tuc/simulation/result_sin/result_{period_real}s/'
    make_plot(path=path,period_real=period_real,cts_range=[1.0,2.0,3.0,4.0],amp_range=[0.3,0.4,0.5,0.6,0.7,0.8],save=0,show=1)
    # est_detnum(path_in=path,cts_range=[1.0,2.0,3.0],amp_range=[0.3,0.4,0.5,0.6,0.7,0.8,0.9])
    # est_bkg()
    # plot_cts_L()