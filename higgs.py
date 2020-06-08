import numpy as np
from numpy import sqrt
# import matplotlib as mpl
# mpl.use('Agg')
import matplotlib.pyplot as plt
from scipy import integrate
import scipy
import time
import pickle
import sys

# this is an index that can be used for batch submission to the cluster
# idx = int(sys.argv[1])
idx = -1

Ne = 100

#Constants
hbar=1
kb=1
e_charge=1

#time constraints
tmin=-4*(2*np.pi)
tmax=8*(2*np.pi)

Nt=200

if Nt != None:
    t_points = np.linspace(tmin, tmax, Nt)
else:
    t_points = None


#superconductor parameters
T = 0.02
wd =  10
s = np.array([1,-1])
m=np.array([1.0,1.2])
ef=np.array([500,300])
gscale=np.array([10,10])
g = np.array([10,0.000000001])
U = np.array([[0.08,0.05],
              [0.05,0.18]])

nb = len(s) #nb=number of bands

#pump parameters
A0 = np.exp(0.5)
A0 = 1
tau = 0.2*np.pi
tau = 2*np.pi
w = 0
w = 1

pre_d0  = np.array([0.47921433, 1.06231377])

# Second pulse
A0_pr = 0
tau_pr = 0.1
w_pr = 0
t_delay = 0

# Calculate remaining dependent parameters from these parameters
try:
    B = 1/(kb*T)
except:
    B = np.inf

kf = np.sqrt(2*m*ef)/hbar
vf = kf/m
n = kf**3/(3*np.pi**2)
N0 = m*kf/(2*np.pi**2)
UN0 = U*N0[:, np.newaxis]

# if U is set, calculate d0 like this:
def d0_integrand(x, d):
    # This is an auxiliary function used in find_d0 to calculate an integral
    return 0.5*1/np.sqrt(x**2+d**2)*np.tanh(B/2*np.sqrt(x**2+d**2))
def find_d0(UN0):
    # this function finds the initial gap(s) given U*N(0). Works in the single or multiband case
    if nb == 2:
        d = np.array([1, 1])
        integral = np.zeros(2)
        for j in [0, 1]:
            integral[j] = integrate.quad(d0_integrand, -wd, wd, (d[j],))[0]
        d_new = np.sum(UN0*d*integral, axis=1)

        while (d != d_new).all(): # seems like convergence of gap is based on machine precision
            d = d_new
            integral = np.zeros(2)
            for j in [0, 1]:
                integral[j] = integrate.quad(d0_integrand, -wd, wd, (d[j],))[0]
            d_new = np.sum(UN0*d*integral, axis=1)
        return d_new
    elif nb == 1:
        d = 1
        d_new = UN0*d*integrate.quad(d0_integrand, -wd, wd, (d,))[0]
        while d != d_new:
            d = d_new
            d_new = UN0*d*integrate.quad(d0_integrand, -wd, wd, (d,))[0]
        return d_new
d_eq0 = find_d0(UN0)

# if d0 is set, calculate U like this:
#U = find_U(pre_d0, 0.05/0.08, integrate.simps(0.5*1/np.sqrt(ep**2+d0.reshape(2,1)**2)*np.tanh(B/2*np.sqrt(ep**2+d0.reshape(2,1)**2)),ep))
#UN0 = U*N0[:,np.newaxis]

if gscale.any() != None:
    g = gscale*(2*d_eq0)
# ep array


ne_ = 3*Ne + 4*Ne**2 # number of equations per band
ne = nb * ne_ # number of equations

# Variables used by ode solvers. They are computed once on script execution and then used many times

# We choose the following naming convention:

# variable name ending on 1_ indicates an 1D array of dimention [Ne] corresponding to energies epsilon in [-wd,wd]
# name ending on 1 indicates an array of dimension [nb,Ne] that also includes the band degree of freedom
# no ending indicates array of dimension [nb,Ne,Ne] where first index is for band, second is for epsilon,
# and the last index is epsilon'

# we will also make use of python broadcasting. If a variable does not depend on an index, we may insert a 
# dimension with only one element using ax. I.e. the 3D array e only depends on the epsilon index, but to make it 
# 3D we insert placeholder axes in the following way e = e1_[ax,:,ax]

ax = np.newaxis

d_eq1 = d_eq0[:,ax]
d_eq = d_eq0[:,ax, ax]
s1 = s[:,ax]
m1 = m[:,ax]
vf1 = vf[:,ax]

e1_ = np.linspace(-wd, wd, Ne)

de = e1_[1] - e1_[0]
de2 = de**2

e1 = e1_[ax,:]
e_ = e1_[:,ax]
ep_ = e1_[ax,:]

e = e1_[ax,:,ax]
ep = e1_[ax,ax,:]

E1 = np.sqrt(e1**2 + d_eq1**2)
E = np.sqrt(e**2 + d_eq**2)
Ep = np.sqrt(ep**2 + d_eq**2)

u1 = np.sqrt(0.5 * (1 + e1/E1))
v1 = np.sqrt(0.5 * (1 - e1/E1))
u = np.sqrt(0.5 * (1 + e/E))
v = np.sqrt(0.5 * (1 - e/E))
up = np.sqrt(0.5 * (1 + ep/Ep))
vp = np.sqrt(0.5 * (1 - ep/Ep))

l = u*up + v*vp
p = v*up - u*vp

W = 1/np.pi*g[:,ax,ax]/((e-ep)**2+g[:,ax,ax]**2)

f1 = 1.0/(np.exp(B*E1) + 1)
f = 1.0/(np.exp(B*E) + 1)
fp = 1.0/(np.exp(B*Ep) + 1)

ff1_ = 1.0/(np.exp(B*e1_) + 1)
ff1 = 1.0/(np.exp(B*e1) + 1)
ff = 1.0/(np.exp(B*e) + 1)
ffp = 1.0/(np.exp(B*ep) + 1)

def integ(x, axis):
    """ Integrate the function 'x' over the axis 'axis'. The integration can be performed over one or two dimensions """
    if hasattr(axis, "__len__"):
        return integrate.simps(integrate.simps(x, dx=de, axis=axis[1]), dx=de, axis=axis[0])
    else:
        return integrate.simps(x, dx=de, axis=axis)
    """
    dx = de
    if hasattr(axis, "__len__"):
        dx = de**(len(axis))
    return np.sum(x, axis=axis) * dx
    """


def A(t):
    """ Returns the vector potential at time t """
    return A0*np.exp(-t**2/(2*tau**2))*np.cos(w*t) #+ A0_pr*np.exp(-(t-t_delay)**2/(2*tau_pr**2))*np.cos(w_pr*(t-t_delay))


def sprime(t, s):
    """ function given to the integrator used to calculate the time evolution of the state vector in the second order calculation """
    ds = np.copy(s).reshape(nb, ne_)
    # unpack and reshape all variables from s
    r11 = ds[:, :Ne]
    r21o = ds[:, Ne:2*Ne]
    r21e = ds[:, 2*Ne:3*Ne]
    R11 = ds[:, 3*Ne:(3*Ne+Ne**2)].reshape(nb,Ne,Ne)
    R21 = ds[:, (3*Ne+Ne**2):(3*Ne+2*Ne**2)].reshape(nb,Ne,Ne)
    f11 = ds[:, (3*Ne+2*Ne**2):(3*Ne+3*Ne**2)].reshape(nb,Ne,Ne)
    f21 = ds[:, (3*Ne+3*Ne**2):(3*Ne+4*Ne**2)].reshape(nb,Ne,Ne)

    # symmetrize
    # r11 = (r11 + np.flip(r11, axis=1))/2
    # r21e = (r21e + np.flip(r21e, axis=1))/2
    # r21o = (r21o - np.flip(r21o, axis=1))/2

    r21 = r21o+r21e

    # self-consitent recalculation of the gap
    d = U @ ( N0*integ(-d_eq1/E1*r11 + u1**2*r21 - v1**2*np.conj(r21), axis=1) )

    # this is for python broadcasting convenience
    d1 = d[:,ax]
    d2 = d[:,ax, ax]

    # first order
    d_f11 = -1j * (Ep - E) * f11 - 1j * (fp - f) * A(t)
    d_f21 = -1j * (Ep + E) * f21 + 1j * (1 - fp - f) * A(t)

    # second order
    d_r11 =  -2*A(t)*(e_charge*vf1)**2/3.0 * integ((l**2*f11.imag - p**2 * f21.imag) * W, axis=2)

    d_r21o = -1j*2*E1 * r21o + 1j*(1-2*f1)*(e1/E1) * d1.real \
        + 1j*2*A(t)*(e_charge*vf1)**2/3.0 * integ(W * l * p * ( f21 - np.conj(f11)), axis=2)

    d_r21e = -2*1j*E1 * r21e - 1j*(1-2*f1)/E1 * e_charge**2 * A(t)**2 / 2 * d_eq1 * s1/m1 - (1-2*f1)*d1.imag

    # third order
    d_R11 = -1j*(Ep-E)*R11 - 1j*np.transpose(l, axes=(0,2,1)) * (A(t) * (r11[:,ax,:] - r11[:,:,ax]) + d2.real * (d_eq/Ep - d_eq/E) * f11 ) \
        - 1j * np.transpose(p, axes=(0,2,1)) * (A(t) * (r21o[:,ax,:] + r21o[:,:,ax].conj()) + d2.real * (-ep/Ep * f21.conj() - e/E * f21) )

    d_R21 = -1j*(Ep+E)*R21 + 1j*np.transpose(p, axes=(0,2,1)) * (A(t) * (r11[:,ax,:] + r11[:,:,ax]) + d2.real * (d_eq/Ep + d_eq/E) * f21 ) \
        - 1j * np.transpose(l, axes=(0,2,1)) * (A(t) * (r21o[:,ax,:] - r21o[:,:,ax])       + d2.real * ( e/E * f11 - ep/Ep * f11.conj()) )


    # reshape and repack into ds
    ds[:, :Ne] = d_r11
    ds[:, Ne:2*Ne] = d_r21o
    ds[:, 2*Ne:3*Ne] = d_r21e
    ds[:, 3*Ne:(3*Ne+Ne**2)] = d_R11.reshape(nb,Ne**2)
    ds[:, (3*Ne+Ne**2):(3*Ne+2*Ne**2)] = d_R21.reshape(nb,Ne**2)
    ds[:, (3*Ne+2*Ne**2):(3*Ne+3*Ne**2)] = d_f11.reshape(nb,Ne**2)
    ds[:, (3*Ne+3*Ne**2):(3*Ne+4*Ne**2)] = d_f21.reshape(nb,Ne**2)

    return ds.reshape((ne,))


start = time.time()

# initial conditions - everything is zero
s0 = np.zeros(ne, dtype=complex)

# the built in integrator solves for the r values numerically:
sols = integrate.solve_ivp( sprime, (tmin, tmax), s0, t_eval=t_points)

# extracting the solutions from the solver output:
Y = sols.y.reshape(nb, ne_, len(sols.t)).swapaxes(0,2).swapaxes(1,2)
t = sols.t
r11 = Y[:, :, :Ne]
r21o = Y[:, :, Ne:2*Ne]
r21e = Y[:, :, 2*Ne:3*Ne]
r21 = r21o + r21e
R11 = Y[:, :, 3*Ne:(3*Ne+Ne**2)].reshape(len(t),nb,Ne,Ne)
R21 = Y[:, :, (3*Ne+Ne**2):(3*Ne+2*Ne**2)].reshape(len(t),nb,Ne,Ne)
F11 = Y[:, :, (3*Ne+2*Ne**2):(3*Ne+3*Ne**2)].reshape(len(t),nb,Ne,Ne)
F21 = Y[:, :, (3*Ne+3*Ne**2):(3*Ne+4*Ne**2)].reshape(len(t),nb,Ne,Ne)

# compute first order currents
jp_1 = e_charge**2*np.sum((n/m) * integ(W*(l**2 * F11.real + p**2 * F21.real), axis=(2,3)), axis=1)

eep = (e - ep) * np.ones((nb,1,1))
for i in np.arange(nb):
    np.fill_diagonal(eep[i], 1)
jd_1 = A(t) * e_charge**2 * np.sum((n/m) * integ( (ff-ffp)  / eep * W, axis=(1,2)))
j_1 = jd_1 + jp_1

# compute e field
efield = A0*np.exp(-t**2/(2*tau**2))\
    *(t/tau**2*np.cos(w*t) + w*np.sin(w*t))\
        +A0_pr*np.exp(-(t-t_delay)**2/(2*tau**2)) \
    * ((t-t_delay)/tau**2*np.cos(w*(t-t_delay)) + w*np.sin(w*(t-t_delay)))


# compute the second order gap from above solutions for all times
d_2 = np.einsum('ij,tj->ti', U, N0*integ(-d_eq1/E1*r11 + u1**2*r21 - v1**2*np.conj(r21), axis=2))


# compute third order current
jp_3 = e_charge**2/2 * np.sum(n/m * integ(W * (l * (R11 + R11.conj()) + p * (R21 + R21.conj())), axis=(2,3)), axis=1)

jd_3 = - e_charge**2 * A(t) * np.sum(s*N0/m * integ((u1**2 - v1**2) * (r11 + r11.conj() ) + 2*u1*v1 * (r21 + r21.conj()), axis=2) )

j_3 = jp_3 + jd_3

end = time.time()
duration = end - start
# pr.disable()
print(f'finished in {duration}s')
# pr.print_stats(sort='time')

tp = t/(2*np.pi)

# plot
if idx == -1:
    plt.figure()
    plt.subplot(121)
    plt.plot(tp,d_2.real)
    plt.title(f'Re$[\Delta]$')
    plt.subplot(122)
    plt.plot(tp, d_2.imag)
    plt.title(f'Im$[\Delta]$')
    plt.tight_layout()


    plt.figure()
    plt.subplot(121)
    plt.plot(tp,jd_1, label=f'$j_D|_1$')
    plt.plot(tp,jp_1, label=f'$j_P|_1$')
    plt.legend()
    plt.subplot(122)
    plt.plot(tp,efield/np.max(efield)*np.max(j_1), '--', label=f'$E$')
    plt.plot(tp,j_1, 'k', label=f'$j|_1$')
    plt.legend()
    plt.tight_layout()


    plt.figure()
    plt.subplot(121)
    plt.plot(tp,jd_3, label=f'$j_D|_3$')
    plt.plot(tp,jp_3, label=f'$j_P|_3$')
    plt.legend()
    plt.subplot(122)
    plt.plot(tp,efield/np.max(efield)*np.max(j_3), '--', label=f'$E$')
    plt.plot(tp,j_3, 'k', label=f'$j|_3$')
    plt.legend()
    plt.tight_layout()


# save as dictionary using pickle
res = {'Ne': Ne,
        'T': T,
        'wd': wd,
        's': s,
        'm': m,
        'ef': ef,
        'gscale': gscale,
        'g': g,
        'U': U,
        'nb': nb,
        'A0': A0,
        'tau': tau,
        'w': w,
        't': t,
        'jp_1': jp_1,
        'jd_1': jd_1,
        'd_2': d_2,
        'jp_3': jp_3,
        'jd_3': jd_3,
        'efield': efield,
        'duration': duration,
        'hbar': hbar,
        'kb': kb,
        'e_charge': e_charge}


f1 = open(str(idx) + '.pickle', 'ab')
pickle.dump(res, f1)
f1.close()