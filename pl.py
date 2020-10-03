import pickle
import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy
from scipy.fftpack import fftfreq, fftshift
import os

pi2 = np.pi*2
folder = '8_single_hr'
# folder = '8_single_hr_no_pump'
folder = '7_pp_long_hr'
folder = '9_Ne'
files = glob.glob(f'{folder}/*.pickle')

def nmax(x):
    return np.max(np.abs(x))

def fft(t,f, inverse=False):
    T = max(t) - min(t)
    Nt = len(t)
    dt = T/Nt

    xw = np.arange(Nt) * 2 * np.pi / T


    if inverse:
        xw = np.concatenate([xw[:Nt//2]+2*np.pi/dt, xw[Nt//2:]])
        idx = np.argsort(xw)
        xw = xw[idx]
        f = f[idx]
        fw = scipy.ifft(f, axis=0)/np.sqrt(Nt)*Nt
    else:
        fw = scipy.fft(f, axis=0)/np.sqrt(Nt)
        xw = np.concatenate([xw[:Nt//2], xw[Nt//2:]-2*np.pi/dt])
        idx = np.argsort(xw)
        xw = xw[idx]
        fw = fw[idx]
    return xw, fw

def next(name):
    i = 0
    while f'{name}_{i}.pdf' in os.listdir(f'{folder}/figs/'):
        i += 1
    return f'{folder}/figs/{name}_{i}.pdf'

res = []
for f in files:
    reader = open(f,'rb')
    try:
        while True:
            a = pickle.load(reader)
            res.append(a)

    except:
        reader.close()


def values(key):
    vals = []
    for r in res:
        vals.append(r[key])
    return list(dict.fromkeys(vals))

def sel(first=False, **kwargs):
    ret = []
    for r in res:
        cond = True
        for key, value in zip(kwargs,kwargs.values()):
            if r[key]!=value:
                cond=False
                break
        if cond:
            ret.append(r)
            if first:
                return r
    return ret

def title(r):
    freq = r['w']
    tau = np.round(r['tau'],2)
    T = r['T']
    g = r['g']
    g1 = g[0]
    if len(g)>1:
        g2 = g[1]
    else:
        g2=0
    A0 = r['A0']
    A0_pr = r['A0_pr']
    w_pr = np.round(r['w_pr'],2)
    tau_pr = np.round(r['tau_pr'],2)
    t_delay = np.round(r['t_delay']/(2*np.pi),2)
    plt.title(f"$\omega={freq}, \\tau={tau}, \gamma_1={g1},$\n $\gamma_2={g2}, T={T}, \Delta t/2\pi={t_delay}$ \n $A_0={A0}, A_0\'={A0_pr}, \\tau\'={tau_pr}, \omega\'={w_pr}$")

def getE(r, t=None):
    A0 = r['A0']
    A0_pr = r['A0_pr']
    tau = r['tau']
    tau_pr = r['tau_pr']
    w = r['w']
    w_pr = r['w_pr']
    if (not hasattr(t, "__len__")) and t==None:
        t = r['t']
    t_delay = r['t_delay']

    efield = A0*np.exp(-t**2/(2*tau**2))\
        *(t/tau**2*np.cos(w*t) + w*np.sin(w*t))\
            +A0_pr*np.exp(-(t-t_delay)**2/(2*tau_pr**2)) \
        * ((t-t_delay)/tau_pr**2*np.cos(w_pr*(t-t_delay)) + w*np.sin(w_pr*(t-t_delay)))

    A = A0*np.exp(-t**2/(2*tau**2))*np.cos(w*t) \
        +  A0_pr*np.exp(-(t-t_delay)**2/(2*tau_pr**2))*np.cos(w_pr*(t-t_delay))

    return A, efield

def plotHiggs(r, t0_by_tau = 2.5):
    d = r['d_2'].real
    t = r['t']
    tau = r['tau']
    e = r['efield']
    w = r['w']
    d_eq = r['d_eq'][:,0,0]
    g = r['g']

    t_ = t[t>t0_by_tau*tau]
    d_ = d[t>t0_by_tau*tau]
    d_ -= np.mean(d_, axis=0)
    plt.figure()
    plt.subplot(121)
    plt.plot(t/pi2,d)
    plt.axvline(t0_by_tau, c='gray', lw=2,ls='--')
    plt.xlabel(f'$t/2 \pi$')
    plt.ylabel(f'Re$\delta \Delta$')
    title(r)

    plt.subplot(122)
    lw = np.copy(g)
    lw[lw>0.01] = 2
    lw[lw<0.01] = 0.8
    plt.axvline(d_eq[0]*2, c='gray', lw=lw[0])
    if len(d_eq)>1:
        plt.axvline(d_eq[1]*2, c='gray', lw=lw[1])
    w_, dw_ = fft(t_, d_)
    plt.plot(w_, np.abs(dw_))
    plt.xlim((0,4))
    plt.xlabel(f'$\omega$')
    plt.tight_layout()


def plotHiggsLast(r, t0_by_tau = 2.5):
    d = r['d_2'].real
    t = r['t']
    tau = r['tau']
    e = r['efield']
    w = r['w']
    d_eq = r['d_eq'][:,0,0]
    g = r['g']

    t_ = t[t>t0_by_tau*tau]
    d_ = d[t>t0_by_tau*tau]
    d_ -= np.mean(d_, axis=0)
    plt.figure()
    plt.subplot(121)
    plt.plot(t/pi2,d)
    plt.axvline(t0_by_tau, c='gray', lw=2,ls='--')
    plt.xlabel(f'$t/2 \pi$')
    plt.ylabel(f'Re$\delta \Delta$')
    title(r)

def plotJ(r, order=1, t0_by_tau = 2.5):
    if order == 1:
        jd = r['jd_1']
        jp = r['jp_1']
        lb='$j_1$'
    else:
        jd = r['jd_3']
        jp = r['jp_3']
        lb = '$j_3$'
    j = jd + jp
    t = r['t']
    tau = r['tau']
    tau_pr = r['tau_pr']
    freq = r['w']
    # e = r['efield']
    A, e = getE(r)
    d_eq = r['d_eq'][:,0,0]
    tmin = np.min(t)
    tmax = np.max(t)
    Nt = len(t)
    t_delay = r['t_delay']
    g = r['g']

    T = tmax - tmin

    plt.figure()
    plt.subplot(121)
    begin = t_delay-3*tau_pr
    plt.axvline(begin/(2*np.pi), c='gray', lw=1,ls='--')
    # sl = np.logical_and(t>begin, t<10*2*np.pi)
    sl = t>begin

    title(r)
    plt.plot(t/pi2, e/nmax(e[sl])*nmax(j[sl]),label='$E$')
    plt.plot(t/pi2,j, label=lb)
    plt.xlabel(f'$t/2 \pi$')
    plt.legend()
    mm = max(np.abs(min(j[sl])),max(j[sl])) * 5
    plt.ylim((-mm*1.1,mm*1.1))


    w, jw = fft(t[sl],j[sl])
    w, ew = fft(t[sl],e[sl])

    #calculate conductivity
    s=jw/ew
    s[w==0] = 0 + 1j*np.inf

    #this is here to find the max of the conductivity so the plots could be normalized
    subs=s[np.abs(w.real-2)<2]
    smax = np.max(subs.real)
    smax2 = nmax(subs.imag)
    plt.subplot(122)
    plt.plot(w.real, np.abs(s.real), label = r'$\Re( \sigma)$')
    plt.plot(w.real, np.abs(s.imag)/smax2*smax, label = r'$\Im( \sigma)$')
    plt.legend()
    lw = np.copy(g)
    lw[lw>0.01] = 2
    lw[lw<0.01] = 0.8
    plt.axvline(d_eq[0]*1, c='gray', lw=lw[0], ls='-.')
    plt.axvline(d_eq[0]*2, c='gray', lw=lw[0], ls='-.')
    plt.axvline(d_eq[0]*3, c='gray', lw=lw[0], ls='-.')
    if len(d_eq)>1:
        plt.axvline(d_eq[1]*1, c='gray', lw=lw[1], ls='-')
        plt.axvline(d_eq[1]*2, c='gray', lw=lw[1])
        plt.axvline(d_eq[1]*3, c='gray', lw=lw[1], ls='-')
    plt.xlim(0,4)
    plt.ylim(0,1.2*smax)
    plt.xlabel(f'$\omega$')

    plt.tight_layout()


def plotLegett(r, t0_by_tau = 3.5):
    d_eq = r['d_eq'][:,0,0]
    d = r['d_2']
    dphase = d[:,0].imag/d_eq[0] - d[:,1].imag/d_eq[1]
    t = r['t']
    tau = r['tau']
    w = r['w']
    g = r['g']

    t_ = t[t>t0_by_tau*tau]
    dp_ = dphase[t>t0_by_tau*tau]
    dp_ -= np.mean(dp_, axis=0)
    plt.figure()
    plt.subplot(121)
    plt.plot(t/pi2,dphase)
    plt.axvline(t0_by_tau, c='gray', lw=2,ls='--')
    # plt.ylim((np.min(dp_),np.max(dp_)))
    plt.ylim((-0.0002,0.0002))
    plt.xlabel(f'$t/2 \pi$')

    plt.subplot(122)
    w_, dpw_ = fft(t_, dp_)
    # plt.axvline(d_eq[0]*2, c='gray', lw=1)
    # plt.axvline(d_eq[1]*2, c='gray', lw=1)

    lw = np.copy(g)
    lw[lw>0.01] = 2
    lw[lw<0.01] = 0.8
    plt.axvline(d_eq[0]*1, c='gray', lw=lw[0], ls='-.')
    plt.axvline(d_eq[0]*2, c='gray', lw=lw[0], ls='-.')
    plt.axvline(d_eq[0]*3, c='gray', lw=lw[0], ls='-.')
    if len(d_eq)>1:
        plt.axvline(d_eq[1]*1, c='gray', lw=lw[1], ls='-')
        plt.axvline(d_eq[1]*2, c='gray', lw=lw[1])
        plt.axvline(d_eq[1]*3, c='gray', lw=lw[1], ls='-')
    plt.plot(w_, np.abs(dpw_))
    plt.xlim((0,4))
    plt.xlabel(f'$\omega$')
    plt.tight_layout()


def plotE(r):
    plt.figure()
    d = r['d_2'].real
    t = r['t']
    tau = r['tau']
    e = r['efield']
    w = r['w']
    d_eq = r['d_eq'][:,0,0]

    plt.subplot(121)
    plt.plot(t/pi2,e)
    plt.ylabel(f'$E$')
    plt.xlabel(f'$t/2 \pi$')
    plt.subplot(122)
    tw, ew = fft(t,e)
    plt.plot(tw, np.abs(ew))
    plt.xlim((0,w + 1/tau*5))
    plt.xlabel(f'$\omega$')
    plt.axvline(d_eq[0], c='gray', lw=1)
    if len(d_eq)>1: plt.axvline(d_eq[1], c='gray', lw=1)
    plt.tight_layout()

delays = np.sort(values('t_delay'))
w_pr = values('w_pr')
Nes = np.sort(values('Ne'))

r0 = sel(T=0.06,Ne=Nes[-1])[0]
d2_ref = r0['d_2'].real
t = r0['t']

xNe = []
yTc = []
yD = []

save = True
save = False
for Ne in Nes[:-1]:
    for r in sel(Ne=Ne,T=0.06):
        # plt.figure()
        diff = (r['d_2'].real - d2_ref)[:,0]
        # plt.plot(t,diff)
        ttt = t[diff>0.000001]
        if len(ttt) == 0:
            tc = max(t)
        else:
            tc = ttt[0]
        # plt.axvline(tc,c='red')
        # plt.ylim(0,0.0001)
        xNe.append(r['Ne'])
        yTc.append(tc)
        yD.append(r['duration'])

plt.figure()
plt.plot(xNe,yTc,'.')
plt.xlabel('Ne')
plt.ylabel('$t_c$ Time after which solution is not converged')
# plt.figure()
# plt.plot(xNe,yD)


"""
for r in sel(t_delay=delays[3],w_pr=w_pr[1],A0_pr=0.005, A0=0):
    # plotE(r)
    if save: plt.savefig(next('efield'))
    plotHiggs(r)
    if save: plt.savefig(next('higgs'))
    plotLegett(r)
    if save: plt.savefig(next('legett'))
    # print(r['duration']/60)
    plotJ(r, order=1)
    if save: plt.savefig(next('j1'))

    plotJ(r, order=3)
    if save: plt.savefig(next('j3'))
"""
