import os, subprocess
from multiprocessing import Pipe, Process
import numpy as np
from astropy.io import fits
from xspec import *

ratioinfo = ('E6.7 keV 6.7 6.6 6.6 6.8 6.8 1e-3',
            'I6.7 "" 1e-3 0. 0. 100. 100. 1e-6',
            'E7.0 keV 6.97 6.87 6.87 7.07 7.07 1e-3',
            'ratio "" 0.1 0. 0. 10. 10. 1e-6')

def feratio(engs, params, flux):
    I70 = params[0]*params[1]*params[3]/params[2]
    f67 = []
    f70 = []
    callModelFunction('gaussian', engs, [params[0], 0], f67)
    callModelFunction('gaussian', engs, [params[2], 0], f70)
    f67 = np.array(f67)
    f70 = np.array(f70)
    f = f67*params[1]+f70*I70
    for i in range(len(flux)):
        flux[i] = f[i]

def init():
    AllModels.addPyMod(feratio, ratioinfo, 'add')

#以上是模型初始化和模型代码，拟合是的模型名字是feratio

def load_spectra(spec):
    global AllData
    for s in spec:
        f = fits.open(s)
        if f['spectrum'].header['exposure'] <= 1:
            f.close()
            continue
        g = fits.open(s[:-6]+'bkg.fits')
        if 'exposure' not in g['spectrum'].header:
            g.close()
            continue
        AllData += s
        f.close()
        g.close()
    AllData.ignore('**-5. 8.-**')
    AllData.ignore('bad')

def load_model():
    m = Model('phabs(apec+gaussian+feratio)')
    m.feratio.norm.frozen = True
    m.phabs.nH.values = 1, 1e-3, 1e-4, 1e-4, 100, 100
    m.apec.Abundanc = 0
    m.apec.kT.values = 10, 1e-3, 1, 1, 86.1738, 86.1738
    m.apec.norm.values = 1e-2, 1e-3, 0, 0, 100, 100
    m.gaussian.LineE = 6.4
    m.gaussian.LineE.frozen = True
    m.gaussian.Sigma = 0
    m.gaussian.Sigma.frozen = True

def fit_spectra(spec, conn):
    f = fits.open(spec)
    name = f['primary'].header['object'].strip()
    f.close()
    if AllData.nSpectra == 0:
        conn.send((name, -1, 0, 0, 0, 0))
        return
    if Fit.dof <= 50:
        conn.send((name, -2, 0, 0, 0, 0))
        return
    Fit.method = ['leven', 100, 1e-6]
    Fit.query = 'yes'
    Fit.perform()
    if 0.5 <= Fit.statistic/Fit.dof <= 2:
        Fit.error('12')
    m = AllModels(1)
    os.system('rm ratio58.xcm')
    Xset.save('ratio58.xcm')
    conn.send((name, m.feratio.ratio.values[0], m.feratio.ratio.error[0]-m.feratio.ratio.values[0], m.feratio.ratio.error[1]-m.feratio.ratio.values[0], Fit.statistic/Fit.dof, Fit.dof))

def clear():
    AllData.clear()

def process(num, conn):
    obs = '%010d' % num
    os.chdir('./%s/analysis' % obs)
    r = subprocess.run('ls *_g.fits', shell= True, stdout=subprocess.PIPE)
    spec = r.stdout.decode('utf8').split('\n')[:-1]
    if spec == []:
        r = subprocess.run('ls *.evt', shell= True, stdout=subprocess.PIPE)
        evt = r.stdout.decode('utf8').split('\n')[:-1]
        if evt != []:
            fevt = fits.open(evt[0])
            name = fevt['primary'].header['object'].strip()
            conn.send((name, -3, 0, 0, 0, 0))
        else:
            conn.send(('no event', -4, 0, 0, 0, 0))
        return
    load_spectra(spec)
    load_model()
    fit_spectra(spec[0], conn)
    clear()

result = []
fl = os.listdir('.')
recv, sender = Pipe()
init()
for f in fl:
    try:
        num = int(f)
        p = Process(target=process, args=(num, sender))
        p.start()
        r = recv.recv()
        result.append((num, r))
        p.join()
    except Exception:
        pass

with open('result58', 'w') as f:
    f.write('Source\tObsId\tFe ratio\tError down\tError up\tChi2(dof)\n')
    for r in result:
        try:
            f.write('{1}\t{0}\t{2:.3f}\t{3:.3f}\t{4:.3f}\t{5:.2f}({6})\n'.format(r[0], r[1][0], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5]))
        except Exception:
            f.write('{1}\t{0}\t{2}\t{3}\t{4}\t{5}({6}) error\n'.format(r[0], r[1][0], r[1][1], r[1][2], r[1][3], r[1][4], r[1][5]))

