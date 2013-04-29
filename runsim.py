"""
runsim.py

Command line frontend to ion molecular dynamics simulations.
"""

import time, datetime, shutil
#import cPickle as pickle
from ctypes import *
from numpy import *
from numpy.random import random, uniform, normal, shuffle
from scipy.linalg import norm
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

def loadLibrary():
    """Loads the C++ library and defines the functions callable
    through Python."""
    dll = CDLL("./ionmd.so")
    dll.simulate.restype = c_int
    dll.simulate.argtypes = [double_p, double_p, POINTER(Params)]
    dll.printIonStatistics.restype = None
    dll.printIonStatistics.argtypes = [POINTER(Params)]
    dll.printParams.restype = None
    dll.printParams.argtypes = [POINTER(Params)]
    return dll

def saveParams(filename, p):
    """Save parameters p to file filename."""
    raise NotImplementedError("Saving parameters to a file is not yet implemented.")

def loadParams(filename):
    """Load and return pickled parameters object."""
    raise NotImplementedError("Loading parameters from a file is not yet implemented.")

def initIons(N, m_lc, Z_lc, N_sc=None, m_sc=None, Z_sc=None):
    """Return length N arrays of masses, charges, and boolean values
    for laser cooling to be given to the simulation. Masses are
    assumed given in amu, charges in units of e. To specify multiple
    sympathetically-cooled ion species, give N_sc, m_sc and Z_sc as
    array-like objects.

    Examples:

    To create 10 ions, 7 of which are laser cooled, and 3 of which are
    the same sympathetically cooled mass:

       m, Z, lc = initIons(10, 138, 1, 3, 136, 1)

    To create 10 ions, 7 laser cooled, 3 of different masses that are
    sympathetically cooled (this doesn't actually work yet):

       m, Z, lc = initIons(10, 138, 1, [1,1,1], [135,136,137], [1,1,1])"""
    m, Z, lc = zeros(N), zeros(N), zeros(N, dtype=c_int)
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
        raise NotImplementedError("This doesn't work yet.")
        N_lc = N - sum(N_sc)
        m[:N_lc] = m_lc*amu
        Z[:N_lc] = Z_lc*q_e
        lc[:N_lc] = 1
        N_i = N_lc
        for i in range(len(N_sc)):
            try:
                m[N_i:N_sc[i]] = m_sc[i]*amu
                Z[N_i:N_sc[i]] = Z_sc[i]*q_e
                N_i += N[i]
            except TypeError:
                raise TypeError("If N_sc is array-like, so must be m_sc and Z_sc!")
            except IndexError:
                raise IndexError("N_sc must have the same dimensions as m_sc and Z_sc!")
    return m, Z, lc

def initialConditions(N, pos_fname=None,
                      randomize=False, rlim=None, zlim=None, vlim=None):
    """Generate intial positions and velocities for N ions. rlim,
    zlim, and vlim are maximum values for initial r displacement, z
    displacement, and velocity, respectively when using
    randomization. If pos_fname and vel_fname are given, read initial
    conditions from those files.

    Loading from a file assumes xyz format for the positions.

    If no filename is given, and randomize is False, arrange ions in
    a chain intially with randomized velocities (limit set by vlim)."""
    x0, v0 = zeros((N,3)), zeros((N,3))
    if pos_fname:
        x0 = loadtxt(pos_fname, skiprows=2, usecols=(1,2,3))*1e-6
        v0 = zeros((N,3))
        if x0.shape != (N,3) or v0.shape != (N,3):
            raise ValueError("Initial conditions do not have the right number of ions.")
        shuffle(x0)
    elif randomize:
        nhat = norm(random((N,3)))
        #x0[:,:2] = uniform(-rlim, rlim, (N,2))
        #x0[:,2] = uniform(-zlim, zlim, (N,))
        x0[:,:2] = normal(0, rlim, size=(N,2))
        x0[:,2] = normal(0, zlim, size=(N,))
        v0 = normal(0, vlim/3, (N,3))
    else:
        x0[:,:2] = uniform(-rlim, rlim, (N,2))
        x0[:,2] = linspace(-zlim, zlim, N)
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

def plotTrajectory(p, dt, t_max, N, traj_file="traj.dat", outfile=None):
    """Plot the saved ion trajectory."""
    if not isinstance(p, Params):
        raise TypeError("p must be a Params instance.")
    traj = fromfile(traj_file, dtype=c_float)
    t = arange(0, t_max, dt)
    traj.shape = ((3,len(t)))
    x, y, z = traj[0,:], traj[1,:], traj[2,:]
    plt.figure()
    plt.subplot(211)
    plt.hold(True)
    plt.plot(t, x, label='$x$')
    plt.plot(t, y, label='$y$')
    plt.ylabel('$x, y$ [mm]')
    plt.legend()
    plt.hold(False)
    plt.subplot(212)
    plt.plot(t, z)
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'$z$ [mm]')
    if not outfile:
        plt.show()

def main(dll, dt, t_max,
         params_file=None, save_final_params=False, N=None, **kwargs):
    """Set up initial conditions and run the simulation. Parameters
    can either be set here or loaded from a file. See saveParams for
    file specifications. Returns the input parameters used (useful for
    plotting routines after running the simulation)."""
    # Load initial parameters from a file
    if params_file:
        p = loadParams(params_file)

    # Defaults
    else:
        # Ion parameters
        if not N:
            N = 250
        if kwargs['all_lc']:
            m, Z, lc = initIons(N, 138, 1)
        else:
            m, Z, lc = initIons(N, 138, 1, int(N*.3), 136, 1)
        m_p = m.ctypes.data_as(double_p)
        Z_p = Z.ctypes.data_as(double_p)
        lc_p = lc.ctypes.data_as(int_p)

        # Simulation parameters
        dt, t_max, traj_start = dt, t_max, 0.
        use_rfmm = 0
        use_coulomb = 1
        use_laser = 1
        use_secular = 0
        use_stochastic = kwargs.get('use_stochastic', 1)

        # Laser parameters
        khat = array([0,0,1])
        khat_p = khat.ctypes.data_as(double_p)
        lmda = 493.5e-9
        r_l = 0.5e-3
        s0 = 10.
        Gamma = 2*pi*15e6
        delta = -Gamma

        # Trap parameters
        r0, z0, kappa = 3.18e-3, 25.4e-3/2., 0.008
        Omega = 2*pi*3.0e6
        V, U, UEC = 235., 0., 300.
        Vsec, wsec = 1, 2*pi*90e3

        # Background gas parameters
        gamma_col = langevinRate(138, 28, 5e-10, 298, 1.71)

        # Data recording parameters
        traj_fname = "traj.dat"
        fpos_fname = "fpos.xyz"
        fvel_fname = "fvel.txt"
        ccd_fname = "ccd.dat"
        record_traj = 1

        # Create Params object to pass to the compiled library
        p = Params(N=N, m=m_p, Z=Z_p, lc=lc_p,
                   khat=khat_p, lmbda=lmda, r_l=r_l,
                   delta=delta, s0=s0, Gamma=Gamma,
                   r0=r0, z0=z0, Omega=Omega,
                   V=V, U=U, UEC=UEC, kappa=kappa,
                   Vsec=Vsec, w=wsec,
                   gamma_col=gamma_col,
                   dt=dt, t_max=t_max,
                   t_steps=len(arange(0, t_max, dt)),
                   use_rfmm=use_rfmm,
                   use_coulomb=use_coulomb,
                   use_laser=use_laser,
                   use_secular=use_secular,
                   use_stochastic=use_stochastic,
                   num_threads=4,
                   quiet=0,
                   traj_fname=traj_fname,
                   fpos_fname=fpos_fname,
                   ccd_fname=ccd_fname,
                   record_traj=record_traj,
                   traj_start=traj_start)
        #saveParams("default.par", p)

    # Initial conditions
    T0, mass = 10e-3, 138*amu
    vlim = sqrt(3*kB*T0/mass)
    if kwargs.has_key('ipos_fname'):
        ipos = kwargs['ipos_fname']
        x0, v0 = initialConditions(N, pos_fname=ipos)
    elif True:
        x0, v0 = initialConditions(N, randomize=True, rlim=r0/3.,
                                   zlim=z0/4., vlim=vlim)
    else:
        x0, v0 = initialConditions(N, rlim=r0/4., zlim=z0 %N, vlim=vlim)
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
    dt, t_max = .05e-6, 2e-3
    dll = loadLibrary()
    if False:
        # 2 ms is plenty sufficient for generating initial conditions
        for N in range(600,601):
            print "N = %d" % N
            p = main(dll, dt, 2e-3, N=N, all_lc=True, print_params=False)
            shutil.copyfile("fpos.xyz", "init/fpos%i.xyz" % N)
    else:
        N = 25
        print "\nMinimizing..."
        main(dll, 0.08e-6, 2e-3, N=N, all_lc=True, print_params=False,
             use_stochastic=0)
        print "\nSimulating..."
        p = main(dll, dt, t_max, N=N, all_lc=False, print_params=True,
                 ipos_fname="fpos.xyz")
                 #ipos_fname="init/fpos%i.xyz" % N)
        #plotTrajectory(p, dt, t_max, N)
        extents = linspace(-600,600,512)
        ionvis.simCCD(bins=extents, clim=[0,5e3])
        ionvis.display()
