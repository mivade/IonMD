"""
runsim.py

Command line frontend to ion molecular dynamics simulations.

This file is part of IonMD.

IonMD is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your
option) any later version.
  
IonMD is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with IonMD.  If not, see <http://www.gnu.org/licenses/>.
"""

import time, datetime, shutil
import ctypes
from ctypes import c_int, c_double, POINTER
import numpy as np
from numpy import pi, sqrt
from numpy.random import random, uniform, normal, shuffle
from numpy.fft import fft, fftshift
from scipy.linalg import norm
import scipy.optimize
import scipy.constants as consts
import matplotlib.pyplot as plt
import ionvis
from params import Params

###############
## CONSTANTS ##
###############

# Physical constants
amu = consts.u
q_e = consts.e
kB = consts.k

# ctypes constants
int_p = POINTER(c_int)
double_p = POINTER(c_double)
double_pp = POINTER(POINTER(c_double))

###############
## FUNCTIONS ##
###############

def load_library():
    """
    Loads the C++ library and defines the functions callable through
    Python.

    """
    dll = ctypes.CDLL("./ionmd.so")
    dll.simulate.restype = c_int
    dll.simulate.argtypes = [double_p, double_p, POINTER(Params)]
    dll.printIonStatistics.restype = None
    dll.printIonStatistics.argtypes = [POINTER(Params)]
    dll.printParams.restype = None
    dll.printParams.argtypes = [POINTER(Params)]
    return dll

def init_ions(N, m_lc, Z_lc, N_sc=None, m_sc=None, Z_sc=None):
    """
    Return length N arrays of masses, charges, and boolean values for
    laser cooling to be given to the simulation. Masses are assumed
    given in amu, charges in units of e. To specify multiple
    sympathetically-cooled ion species, give N_sc, m_sc and Z_sc as
    array-like objects.

    Examples
    --------

    To create 10 ions, 7 of which are laser cooled, and 3 of which are
    the same sympathetically cooled mass::

         m, Z, lc = init_ions(10, 138, 1, 3, 136, 1)

    To create 10 ions, 7 laser cooled, 3 of different masses that are
    sympathetically cooled::

        m, Z, lc = init_ions(10, 138, 1, [1,1,1], [135,136,137], [1,1,1])

    """
    m, Z, lc = np.zeros(N), np.zeros(N), np.zeros(N, dtype=c_int)
    if not N_sc:
        m += m_lc*amu
        Z += Z_lc*q_e
        lc += 1
    elif not isinstance(N_sc, (list, tuple)):
        if N_sc >= N:
            raise ValueError("N_sc must be < N.")
        N_lc = N - N_sc
        m[:N_lc] = m_lc*amu
        m[N_lc:] = m_sc*amu
        Z[:N_lc] = Z_lc*q_e
        Z[N_lc:] = Z_sc*q_e
        lc[:N_lc] = 1
    else:
        N_lc = N - sum(N_sc)
        m[:N_lc] = m_lc*amu
        Z[:N_lc] = Z_lc*q_e
        lc[:N_lc] = 1
        N_i = N_lc
        for i in range(len(N_sc)):
            try:
                m[N_i:(N_i + N_sc[i])] = m_sc[i]*amu
                Z[N_i:(N_i + N_sc[i])] = Z_sc[i]*q_e
                N_i += N_sc[i]
            except IndexError:
                raise IndexError("N_sc must have the same dimensions as m_sc and Z_sc!")
    return m, Z, lc

def initial_conditions(N, pos_fname=None,
                       randomize=False, rlim=None, zlim=None, vlim=None):
    """
    Generate intial positions and velocities for N ions. rlim, zlim,
    and vlim are maximum values for initial r displacement, z
    displacement, and velocity, respectively when using
    randomization. If pos_fname and vel_fname are given, read initial
    conditions from those files.

    Loading from a file assumes xyz format for the positions.

    If no filename is given, and randomize is False, arrange ions in a
    chain intially with randomized velocities (limit set by vlim).

    """
    x0, v0 = np.zeros((N, 3)), np.zeros((N, 3))
    if pos_fname:
        x0 = np.loadtxt(pos_fname, skiprows=2, usecols=(1,2,3))*1e-6
        v0 = np.zeros((N,3))
        if x0.shape != (N,3) or v0.shape != (N,3):
            raise ValueError("Initial conditions do not have the right number of ions.")
        shuffle(x0)
    elif randomize:
        x0[:,:2] = uniform(-rlim, rlim, (N,2))
        x0[:,2] = uniform(-zlim, zlim, (N,))
        v0 = normal(0., vlim, (N,3))
        #v0 = zeros((N,3))
    else:
        x0[:,:2] = uniform(-rlim, rlim, (N,2))
        x0[:,2] = np.linspace(-zlim, zlim, N)
        v0 = uniform(-vlim, vlim, (N,3))
        shuffle(x0)
    return x0, v0

def langevinRate(m_ion, m_gas, P, T, alpha):
    """Compute the Langevin collision rate between an ion of mass
    m_ion and neutrals of mass m_gas at pressure P, temperature T, and
    polarizability alpha. Masses are given in amu, pressure in torr,
    temperature in K and alpha in angstrom**3. This assumes a charge
    of +e on the ion for now."""
    mu = amu*1000*m_ion*m_gas/(m_ion + m_gas)
    P = P*1333.2239
    alpha = alpha*1e-24
    Ze = 4.8032043e-10
    return 2*pi*Ze*sqrt(alpha/mu)*P/(kB*1e7*T)

def plot_trajectory(dt, t_max, N, traj_file="traj.dat",
                    outfile=None, start=0, end=1000):
    """Plot the saved ion trajectory."""
    traj = np.fromfile(traj_file, dtype=c_float)
    t = (np.arange(0, t_max, dt)/1e-6)[start:]
    traj.shape = (len(traj)/3,3)
    x, y, z = traj[:,0], traj[:,1], traj[:,2]
    plt.figure()
    plt.subplot(2,1,1)
    plt.hold(True)
    plt.plot(t[:end], x[:end], '-', label='$x$')
    plt.plot(t[:end], y[:end], '-', label='$y$')
    plt.ylabel('$x, y$ [mm]')
    plt.legend()
    plt.hold(False)
    plt.subplot(2,1,2)
    plt.plot(t[:end], z[:end], '-')
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'$z$ [mm]')
    if not outfile:
        plt.show()

def plot_fourier(dt, t_max, N, traj_file="traj.dat",
                 outfile=None, start=0, end=1000):
    """
    Plot the Fourier transform of the trajectory data to extract
    motional frequencies.

    """
    traj = np.fromfile(traj_file, dtype=c_float)
    t = (np.arange(0, t_max, dt)/1e-6)[start:]
    traj.shape = (len(traj)/3,3)
    x, y, z = traj[:,0], traj[:,1], traj[:,2]
    plt.figure()
    plt.subplot(2,1,1)
    plt.hold(True)
    plt.plot(fft(x[:end]), '-', label='$x$')
    plt.plot(fft(y[:end]), '-', label='$y$')
    plt.hold(False)
    plt.legend()
    plt.hold(False)
    plt.subplot(2,1,2)
    plt.plot(fft(z[:end]), '-')
    if not outfile:
        plt.show()

def plot_temperature(N, m=138*amu, temp_file="temperature.txt", outfile=None):
    """Plot the temperature over time."""
    t, v = np.loadtxt(temp_file, unpack=True)
    t /= 1e-6
    T = m*v**2/(3*N*kB)
    plt.figure()
    plt.plot(t, T)
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'T [???]')
    if not outfile:
        plt.show()

def main(dll, dt, t_max,
         params_file=None, save_final_params=False, N=None, **kwargs):
    """
    Set up initial conditions and run the simulation. Returns the
    input parameters used (useful for plotting routines after running
    the simulation).

    """
    # Load initial parameters from a file
    if params_file:
        pass # TODO

    # Defaults
    else:
        # Ion parameters
        if not N:
            N = 250
        if kwargs['all_lc']:
            m, Z, lc = init_ions(N, 138, 1)
            masses = np.array([138*amu])
        else:
            masses = np.array([138, 136])*amu
            m, Z, lc = init_ions(N, masses[0]/amu, 1, int(N*.1), masses[1]/amu, 1)
        try: # use passed ion parameters if given
            masses = kwargs['masses']
            m = kwargs['m']
            Z = kwargs['Z']
            lc = kwargs['lc']
        except KeyError:
            print "You didn't give masses, m, Z, or lc as a kwarg!"
        N_masses = len(masses)
        m_p = m.ctypes.data_as(double_p)
        Z_p = Z.ctypes.data_as(double_p)
        masses_p = masses.ctypes.data_as(double_p)
        lc_p = lc.ctypes.data_as(int_p)

        # Laser parameters
        khat = np.array([0,0,-1], dtype=np.float64)
        khat /= norm(khat)
        khat_p = khat.ctypes.data_as(double_p)
        lmda = 493.5e-9
        r_l = 0.5e-3
        s = 3.
        Gamma = 2*pi*15e6
        delta = -2*Gamma
        beta = kwargs.get('beta', 2e-22)
        F0 = kwargs.get('F0', 1.3e-19)

        # Trap parameters
        r0, z0 = 3.18e-3, (25e-3)/2.
        kappa = kwargs.get('kappa', 0.006) #.008
        Omega = kwargs.get('Omega', 2*pi*2.7e6)
        V, U, UEC = kwargs.get('V', 125.), 0., kwargs.get('UEC', 300.)
        Vsec, wsec = 1, 2*pi*90e3

        # Background gas parameters
        #gamma_col = langevinRate(138, 28, 5e-9, 298, 1.71)
        gamma_col = kwargs.get('gamma_col', 1.) #11.95 # K/s
        m_gas = 28*amu
        T_gas = 298.

        # CCD settings
        sim_ccd = kwargs.get('sim_ccd', 1)
        ccd_bins = kwargs.get('ccd_bins', 500)
        ccd_extent = kwargs.get('ccd_extent', 1000)

        # Simulation parameters
        dt, t_max, traj_start = dt, t_max, kwargs.get('traj_start', 0.)
        min_time = kwargs.get('min_time', 2e-3)
        use_rfmm = kwargs.get('use_rfmm', 0)
        use_coulomb = kwargs.get('use_coulomb', 1)
        use_laser = kwargs.get('use_laser', 1)
        use_secular = kwargs.get('use_secular', 0)
        use_stochastic = kwargs.get('use_stochastic', 1)
        abort_bounds = r0
        use_abort = 1

        # Data recording parameters
        traj_fname = kwargs.get("traj_fname", "traj.dat")
        com_fname = kwargs.get("com_fname", "com_traj.dat")
        fpos_fname = "fpos.xyz"
        fvel_fname = "fvel.txt"
        ccd_fname = "ccd"
        temp_fname = "temperature.txt"
        T_steps = kwargs.get('T_steps', 1000)
        record_traj = 1

        # Create Params object to pass to the compiled library
        p = Params(N=N, N_masses=len(masses),
                   m=m_p, Z=Z_p, masses=masses_p, lc=lc_p,
                   khat=khat_p, lmbda=lmda, r_l=r_l,
                   delta=delta, s=s, Gamma=Gamma, beta=beta, F0=F0,
                   r0=r0, z0=z0, Omega=Omega,
                   V=V, U=U, UEC=UEC, kappa=kappa,
                   Vsec=Vsec, w=wsec,
                   gamma_col=gamma_col, m_gas=m_gas, T_gas=T_gas,
                   sim_ccd=sim_ccd,
                   ccd_bins=ccd_bins, ccd_extent=ccd_extent,
                   dt=dt, t_max=t_max, min_time=min_time,
                   abort_bounds=abort_bounds,
                   t_steps=len(np.arange(0, t_max, dt)),
                   use_rfmm=use_rfmm,
                   use_coulomb=use_coulomb,
                   use_laser=use_laser,
                   use_secular=use_secular,
                   use_stochastic=use_stochastic,
                   use_abort=use_abort,
                   num_threads=4,
                   quiet=0,
                   com_fname=com_fname,
                   traj_fname=traj_fname,
                   fpos_fname=fpos_fname,
                   ccd_fname=ccd_fname,
                   temp_fname=temp_fname,
                   record_traj=record_traj,
                   traj_start=traj_start,
                   T_steps=T_steps)

    # Initial conditions
    T0 = 10e-3
    if kwargs.has_key('ipos_fname'):
        ipos = kwargs['ipos_fname']
        x0, v0 = initial_conditions(N, pos_fname=ipos)
    else:
        x0, v0 = initial_conditions(N, randomize=True, rlim=r0/2.,
                                    zlim=z0, vlim=sqrt(3*kB*T0/masses[0]/2.))
    x0 = x0.flatten()
    v0 = v0.flatten()
    x0_p = x0.ctypes.data_as(double_p)
    v0_p = v0.ctypes.data_as(double_p)

    # Run the simulation
    print ":::", time.ctime(), ":::"
    if kwargs.get('print_params', False):
        print "=== Ion statistics ==="
        if N <= 10:
            dll.printIonStatistics(p)
        else:
            print "Lots of ions."
        print "\n=== Simulation parameters ==="
        dll.printParams(p)
    print "\n=== Begin simulation ==="
    t0 = time.time()
    dll.simulate(x0_p, v0_p, p)
    print "Finish on", time.ctime()
    t_total = time.time() - t0
    print "Total run time: %s" % str(datetime.timedelta(seconds=t_total))
    return p

##########
## MAIN ##
##########

if __name__ == "__main__":
    dll = load_library()
    dt, t_max = 20e-9, 5e-3
    min_time = 1.5e-3
    traj_start = dt*100
    ccd_bins, ccd_extent = 768, 1024
    kappa = 9e-3
    all_lc = False
    gcol, beta, F0 = 8., 4e-22, 5.0e-20 # F0 corresponds to s ~ 4
    V = 120.
    N = 50
    masses = np.array([138, 136])*amu
    N_ccd = len(masses)
    if False:
        t0 = time.time()
        Nsc = int(N*0.1)
        m, Z, lc = init_ions(N, masses[0]/amu, 1., 
                             Nsc, masses[1]/amu, 1.)
        for UEC in np.arange(10, 101, 10):
            prefix = "UEC%.0f" % UEC
            p = main(dll, dt, t_max, min_time=min_time,
                     N=N, all_lc=all_lc, print_params=False,
                     masses=masses, m=m, Z=Z, lc=lc,
                     ccd_bins=ccd_bins, ccd_extent=ccd_extent,
                     V=V, UEC=UEC,
                     kappa=kappa,
                     gamma_col=gcol, beta=beta, F0=F0,
                     use_stochastic=1,
                     T_steps=1200,
                     traj_start=traj_start)
            ## for i in range(1,N_ccd+1):
            ##     fname = "ccd_%d.dat" % i
            ##     shutil.copy(fname, "data/"+prefix+fname)
            ## shutil.copy('fpos.xyz', "data/"+prefix+".xyz")
            ## shutil.copy('traj.dat', "data/"+prefix+"_traj.dat")
            shutil.copy('com_traj.dat', "data/"+prefix+"_com_traj.dat")
        t = time.time() - t0
        print "Finished all in: %s" % str(datetime.timedelta(seconds=t))
    else:
        if True:
            p = main(dll, dt, t_max, min_time=min_time,
                     N=N, all_lc=all_lc, print_params=True,
                     ccd_bins=ccd_bins, ccd_extent=ccd_extent,
                     V=V,
                     kappa=kappa,
                     gamma_col=gcol, beta=beta, F0=F0,
                     use_stochastic=1,
                     T_steps=1200,
                     traj_start=traj_start)
        #plot_trajectory(dt, t_max, N, start=traj_start/dt, end=-1)
        #plot_fourier(dt, t_max, N, start=traj_start/dt, end=-1)
        #ionvis.display(fpos_fname="ipos.xyz")
        ionvis.display()#outfile='images/test.png')
        #plot_temperature(N, 138*amu)
        if False:
            #for N_ccd in range(5,0,-1):
            if all_lc:
                N_ccd = 1
            else:
                N_ccd = 2
            ionvis.simCCD("ccd", N_ccd, ccd_bins, ccd_extent,
                          outfile="images/CCD_latest.png",
                          show=True, brightness=2)
            
