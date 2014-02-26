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

from __future__ import print_function
from __future__ import division
import time, datetime
import ctypes
import numpy as np
from numpy import pi, sqrt
from numpy.random import uniform, normal, shuffle
import scipy.constants as consts
import ionvis
from params import Params as CxxParams
from settings import SimParams
from paultrap import PaulTrap

###############
## CONSTANTS ##
###############

# Physical constants
amu = consts.u
q_e = consts.e
kB = consts.k

# ctypes constants
_DOUBLE_P = ctypes.POINTER(ctypes.c_double)

###############
## FUNCTIONS ##
###############

def load_library():
    """
    Loads the C++ library and defines the functions callable through
    Python.

    """
    dll = ctypes.CDLL("./libionmd.so")
    dll.simulate.restype = ctypes.c_int
    dll.simulate.argtypes = [_DOUBLE_P, _DOUBLE_P, ctypes.POINTER(CxxParams)]
    dll.print_ion_statistics.restype = None
    dll.print_ion_statistics.argtypes = [ctypes.POINTER(CxxParams)]
    dll.print_params.restype = None
    dll.print_params.argtypes = [ctypes.POINTER(CxxParams)]
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

        m, Z, lc = init_ions(10, 138, 1, [1, 1, 1],
                             [135, 136, 137], [1, 1, 1])

    """
    m, Z, lc = np.zeros(N), np.zeros(N), np.zeros(N, dtype=ctypes.c_int)
    if not N_sc:
        m += m_lc
        Z += Z_lc
        lc += 1
    elif not isinstance(N_sc, (list, tuple)):
        if N_sc >= N:
            raise ValueError("N_sc must be < N.")
        N_lc = N - N_sc
        m[:N_lc] = m_lc
        m[N_lc:] = m_sc
        Z[:N_lc] = Z_lc
        Z[N_lc:] = Z_sc
        lc[:N_lc] = 1
    else:
        N_lc = N - sum(N_sc)
        m[:N_lc] = m_lc
        Z[:N_lc] = Z_lc
        lc[:N_lc] = 1
        N_i = N_lc
        for i in range(len(N_sc)):
            try:
                m[N_i:(N_i + N_sc[i])] = m_sc[i]
                Z[N_i:(N_i + N_sc[i])] = Z_sc[i]
                N_i += N_sc[i]
            except IndexError:
                raise IndexError("N_sc must have the same dimensions as " + \
                                 "m_sc and Z_sc!")
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
        x0 = np.loadtxt(pos_fname, skiprows=2, usecols=(1, 2, 3))*1e-6
        v0 = np.zeros((N, 3))
        if x0.shape != (N, 3) or v0.shape != (N, 3):
            raise ValueError("Initial conditions do not have the right " + \
                             "number of ions.")
        shuffle(x0)
    elif randomize:
        x0[:,:2] = uniform(-rlim, rlim, (N, 2))
        x0[:,2] = uniform(-zlim, zlim, (N,))
        v0 = normal(0., vlim, (N, 3))
        #v0 = zeros((N,3))
    else:
        x0[:,:2] = uniform(-rlim, rlim, (N, 2))
        x0[:,2] = np.linspace(-zlim, zlim, N)
        v0 = uniform(-vlim, vlim, (N, 3))
        shuffle(x0)
    return x0, v0

def langevin_rate(m_ion, m_gas, P, T, alpha):
    """
    Compute the Langevin collision rate between an ion of mass m_ion
    and neutrals of mass m_gas at pressure P, temperature T, and
    polarizability alpha. Masses are given in amu, pressure in torr,
    temperature in K and alpha in angstrom**3. This assumes a charge
    of +e on the ion for now.
    """
    mu = amu*1000*m_ion*m_gas/(m_ion + m_gas)
    P = P*1333.2239
    alpha = alpha*1e-24
    Ze = 4.8032043e-10
    return 2*pi*Ze*sqrt(alpha/mu)*P/(kB*1e7*T)

def run_simulation(sim_settings, **kwargs):
    """
    Run the simulation.

    Keyword Arguments
    -----------------
    quit_after_minimizing : bool
        If True, quit execution after running the minimization
        procedure. Mostly useful for debugging purposes to verify
        things like time step size are sufficiently small.

    """
    dll = load_library()
    
    # Check that the parameters were correctly given.
    if not isinstance(sim_settings, SimParams):
        raise TypeError("sim_settings must be an instance of " + \
                        "settings.SimParams.")

    # C++ settings struct
    p = sim_settings.to_cxx_format()
    try:
        p.quit_after_minimizing = ctypes.c_int(kwargs['quit_after_minimizing'])
    except KeyError:
        p.quit_after_minimizing = 0

    # Initial conditions
    # TODO: clean this up, maybe allow for settings
    T0 = 10e-3
    N = p.N
    r0 = p.r0
    z0 = p.z0
    m = sim_settings.ions['m'][0]*amu
    if 'ipos_fname' in kwargs:
        ipos = kwargs['ipos_fname']
        x0, v0 = initial_conditions(N, pos_fname=ipos)
    else:
        x0, v0 = initial_conditions(N, randomize=True, rlim=r0/2.,
                                    zlim=z0, vlim=sqrt(3*kB*T0/m/2.))
    x0 = x0.flatten()
    v0 = v0.flatten()
    print(x0, v0)
    x0_p = x0.ctypes.data_as(_DOUBLE_P)
    v0_p = v0.ctypes.data_as(_DOUBLE_P)

    # Run the simulation
    print(":::", time.ctime(), ":::")
    print("=== Ion statistics ===")
    if N <= 10:
        dll.print_ion_statistics(p)
    else:
        print("Lots of ions.")
    print("\n=== Simulation parameters ===")
    dll.print_params(p)
    print("\n=== Begin simulation ===")
    t0 = time.time()
    dll.simulate(x0_p, v0_p, p)
    print("Finish on", time.ctime())
    t_total = time.time() - t0
    print("Total run time: %s" % str(datetime.timedelta(seconds=t_total)))

##########
## MAIN ##
##########

if __name__ == "__main__":
    parameters = SimParams("default.json")
    #m, Z, lc = init_ions(20, 40, 1, 2, 42, 1)
    N = 20
    m = N*[40]
    Z = N*[1]
    lc = N*[1]
    parameters.set_ions(m, Z, lc)
    
    trap = PaulTrap(parameters.trap['r0'], parameters.trap['f_RF'])
    q, a = trap.get_mathieu_qa(parameters.ions['m'][0]*amu,
                               parameters.trap['V'],
                               parameters.trap['U'],
                               parameters.ions['Z'][0])
    print("q =", q, ", a =", a)

    parameters.control['display'] = False
    run_simulation(parameters, quit_after_minimizing=False)

    dt = parameters.control['dt']
    t_max = parameters.control['t_max']
    N = len(parameters.ions['m'])

    if parameters.control['sim_ccd']:
        ionvis.sim_ccd(parameters.files['ccd_fname'],
                       len(np.unique(parameters.ions['m'])),
                       parameters.ccd['ccd_bins'],
                       parameters.ccd['ccd_extent'],
                       show=True, brightness=2)
    if parameters.control['plot_trajectory']:
        ionvis.plot_trajectory(dt, t_max, N)
    if parameters.control['plot_fourier']:
        pass
        #ionvis.plot_fourier(dt, t_max, N, start=traj_start/dt, end=-1)
    if parameters.control['plot_temperature']:
        pass
        #plot_temperature(N, 138*amu)
    if parameters.control['display']:
        #ionvis.display("ipos.xyz")
        ionvis.display("fpos.xyz", scale=1)
