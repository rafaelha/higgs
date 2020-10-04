
    ps = np.loadtxt('..\\material_parameters_Matsunaga\\pulse_shape.csv', delimiter=',')
    rm = 20
    t_pulse = running_mean(ps[:,0],rm)

    b = 0.3
    smooth = (np.tanh((t_pulse-1.5)/b)+1)/2
    smooth2 = (-np.tanh((t_pulse-7.5)/b)+1)/2

    e_pulse = running_mean(ps[:,1],rm)*smooth*smooth2
    e_pulse /= np.max(np.abs(e_pulse))

    # 1 time unit = 0.438856 ps
    t_pulse /= 0.438856

    # plt.figure('poly')
    # plt.clf()
    # plt.plot(t_pulse, e_pulse, '.-')
    # xx = np.linspace(min(t_pulse), max(t_pulse), 2000)
    # plt.plot(xx, poly(xx))


    a_pulse = -scipy.integrate.cumtrapz(e_pulse,t_pulse, initial=0)
    shift = 0.9
    A_pump = scipy.interpolate.interp1d(np.concatenate([[-1e4],t_pulse-shift*2*np.pi,[1e4]]),A0*np.concatenate([[0],a_pulse,[a_pulse[-1]]]),kind='cubic')
    A_probe = scipy.interpolate.interp1d(np.concatenate([[-1e4],t_pulse-shift*2*np.pi+t_delay,[1e4]]),A0_pr*np.concatenate([[0],a_pulse,[a_pulse[-1]]]),kind='cubic')

    # plt.plot(xx,A_pump(xx))
    # efield = -np.diff(A_pump(xx), prepend=0)/(xx[1]-xx[0])
    # plt.plot(xx,efield, '.-')

    def A(t):
        return A_pump(t) + A_probe(t)
        """ Returns the vector potential at time t """
        # return A0*np.exp(-t**2/(2*tau**2))*np.cos(w*t) \
            # +  A0_pr*np.exp(-(t-t_delay)**2/(2*tau_pr**2))*np.cos(w_pr*(t-t_delay))

    if False:
        plt.figure('A')
        plt.clf()
        plt.subplot(121)
        # plt.plot(tp,A(t))

        efield = -np.diff(A(t))/(t[1]-t[0])
        efield = np.concatenate([[efield[0]],efield])
        # efield = np.concatenate([[0],efield[:-1]])
        plt.plot(tp,efield)
        plt.xlim((-3,3))
        plt.ylabel(f'$A$')
        plt.xlabel(f'$t/2 \pi$')

        def nm(x):
            return x / np.max(np.abs(x))
        plt.subplot(122)
        tw, ew = pfft(t,efield)
        plt.plot(tw, np.real(nm(ew)))
        plt.plot(tw, np.imag(nm(ew)))
        plt.plot(tw, np.abs(nm(ew)),'.-')
        plt.xlim((0,5*d_eq0[0]))
        plt.xlabel(f'$\omega$')
        plt.axvline(d_eq[0], c='gray', lw=1)
        if len(d_eq)>1: plt.axvline(d_eq[1], c='gray', lw=1)