"""
settings.py

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
import json
import numpy as np
from numpy import pi
from params import Params
import ctypes
import scipy.constants as consts

_DOUBLE_P = ctypes.POINTER(ctypes.c_double)
_INT_P = ctypes.POINTER(ctypes.c_int)
_AMU = consts.u
_ECHARGE = consts.e

class SimParams(object):
    """
    Pythonic representation of all simulation parameters. This allows
    for reading and writing of simulation setup to JSON files for easy
    storage of parameters used and helps clean up the Python front end
    to the MD simulations.

    Attributes
    ----------
    control : dict
        Simulation control parameters.
    ions : dict
        Ion simulation parameters.
    trap : dict
        Trap simulation parameters.
    stochastic : dict
        Simulation parameters for stochastic forces.
    laser : dict
        Laser simulation parameters.
    ccd : dict
        Parameters for CCD image smiulation.
    files : dict
        File names for data files.
    cxx_params : Params
        ctypes Struct object for simulation parameters to pass to the
        C++ library.

    """

    def __init__(self, init_file="default.json"):
        """
        Create the SimParams object using a default settings file if
        none is specified. Using a default file is helpful so that you
        can vary some parameters, but keep the rest the same by
        defining them only once.

        """
        self.control, self.ions, self.trap = {}, {}, {}
        self.stochastic, self.laser, self.ccd, self.files = {}, {}, {}, {}
        self.load(init_file)
        self.cxx_params = None

    def set_control(self, **kwargs):
        """
        Set simulation control parameters.

        Keyword arguments
        -----------------
        use_rfmm : bool
            Use RF micromotion in the simulation (not yet
            implemented) if True.
        use_coulomb : bool
            Include the Coulomb interaction if True.
        use_laser : bool
            Include laser cooling if True.
        use_stochastic : bool
            Include stochastic processes if True.
        sim_ccd : bool
            Simulate CCD images if True.
        dt : float
            Time step size in s.
        t_max : float
            Maximum simulation time (in sim time units).
        min_time : float
            Simulation time to run the minimization procedure.
        plot_trajectory : bool
            Plot the COM trajectory after the simulation is done if
            True.
        plot_fourier : bool
            Plot the Fourier transform of the COM when the simulation
            is done if True.
        display : bool
            Use MayaVI to display the final ion positions if True.
        num_threads : int
            Number of threads to use for multithreaded calculations.

        """
        self.control['use_rfmm'] = kwargs.get('use_rfmm', False)
        self.control['use_coulomb'] = kwargs.get('use_colomb', True)
        self.control['use_laser'] = kwargs.get('use_laser', True)
        self.control['sim_ccd'] = kwargs.get('sim_ccd', True)
        self.control['dt'] = kwargs.get('dt', 20e-9)
        self.control['t_max'] = kwargs.get(5e-3)
        self.control['min_time'] = kwargs.get(1.5e-3)
        self.control['plot_trajectory'] = kwargs.get('plot_trajectory', False)
        self.control['plot_fourier'] = kwargs.get('plot_fourier', False)
        self.control['display'] = kwargs.get('display', True)
        self.control['num_threads'] = kwargs.get('num_threads', 1)

    def set_files(self, traj_fname="traj.dat", com_fname="com_traj.dat",
                  fpos_fname="fpos.xyz", ccd_fname="ccd",
                  temp_fname="temperature.txt"):
        """
        Set the filenames (or prefix for filenames in the case of
        simulated CCD files.

        Parameters
        ----------
        traj_fname : str
            Filename for trajectory data of the ion being tracked.
        fpos_fname : str
            Filename for xyz storage of final ion positions.
        com_fname : str
            Filename for the center of mass trajectory data.
        ccd_fname : str
            Prefix for simulated CCD filenames.
        temp_fname : str
            Filename for ion temperature data.
        
        """
        self.files['traj_fname'] = traj_fname
        self.files['com_fname'] = com_fname
        self.files['fpos_fname'] = fpos_fname
        self.files['ccd_fname'] = ccd_fname
        self.files['temp_fname'] = temp_fname

    def set_ions(self, m, Z, lc):
        """
        Set ion parameters.

        Parameters
        ----------
        m : list
            List of masses of every ion in amu.
        Z : list
            List of charges in units of e for each ion.
        lc : list
            List of booleans or ints. For each ion, True indicates
            laser cooling and False indicates no laser cooling.

        """
        if len(m) is len(Z) and len(Z) is len(lc):
            self.ions['m'] = m
            self.ions['Z'] = Z
            self.ions['lc'] = lc
        else:
            raise ValueError("m, Z, and lc must have the same length.")

    def set_trap(self, **kwargs):
        """
        Set trap parameters.

        Keyword arguments
        -----------------
        r0 : float
            Trap radius.
        z0 : float
            Half the distance between the end cap electrodes.
        kappa : float
            Trap geometrical factor.
        V : float
            RF voltage.
        U : float
            RF DC offset.
        UEC : float
            End cap voltage.

        """
        self.trap['r0'] = kwargs.get('r0', 3e-3)
        self.trap['z0'] = kwargs.get('z0', 25e-3)
        self.trap['kappa'] = kwargs.get('kappa', 1e-4)
        self.trap['f_RF'] = kwargs.get('f_RF', 4.7e6)
        self.trap['V'] = kwargs.get('V', 100)
        self.trap['U'] = kwargs.get('U', 0)
        self.trap['UEC'] = kwargs.get('UEC', 30)

    def set_laser(self, **kwargs):
        """
        Set laser parameters.

        Keyword arguments
        -----------------
        khat : length 3 list
            Unit vector pointing in the direction of laser
            propagation.
        lmbda : float
            Cooling laser wavelength.
        Gamma_Hz : float
            Cooling transition linewidth in Hz.
        delta_Gamma : float
            Detuning in units of Gamma_Hz.
        s : float
            Saturation parameter I/I_s.
        beta : float
            Damping parameter for laser cooling.
        F0 : float
            Constant laser radiation pressure force.
        
        """
        self.laser['khat'] = kwargs.get('khat', [0., 0., 1.])
        self.laser['lmbda'] = kwargs.get('lmbda', 397e-9)
        self.laser['Gamma_Hz'] = kwargs.get('Gamma_Hz', 20e6)
        self.laser['delta_Gamma'] = kwargs.get('delta_gamma', -10.)
        self.laser['s'] = kwargs.get('s', 5.0)
        self.laser['beta'] = kwargs.get('beta', 2e-22)
        self.laser['F0'] = kwargs.get('F0', 1.3e-19)

    def set_stochastic(self, gamma_col=1.0, m_gas_amu=28., T_gas=298.):
        """
        Set stochastic force parameters.

        Parameters
        ----------
        gamma_col : float
            Background gas collision rate.
        m_gas_amu : float
            Background gas mass in amu.
        T_gas : float
            Temperature of the background gas.

        """
        self.stochastic['gamma_col'] = gamma_col
        self.stochastic['m_gas_amu'] = m_gas_amu
        self.stochastic['T_gas'] = T_gas

    def set_ccd(self, ccd_bins=768, ccd_extent=1024):
        """
        Set simulated CCD parameters.

        Parameters
        ----------
        ccd_bins : int
            Number of horizontal and vertical bins for the CCD.
        ccd_extent : numeric
            Extent of the CCD in microns.

        """
        self.ccd['ccd_bins'] = ccd_bins
        self.ccd['ccd_extent'] = ccd_extent

    def save(self, fname):
        """Write parameters to JSON file fname."""
        obj = {"control": self.control,
               "files": self.files,
               "params": {"ions": self.ions,
                          "trap": self.trap,
                          "laser": self.laser,
                          "stochastic": self.stochastic,
                          "ccd": self.ccd}}
        with open(fname, 'w') as json_file:
            json.dump(obj, json_file, indent=4)

    def load(self, fname):
        """Read parameters from JSON file fname."""
        with open(fname, 'r') as json_file:
            params = json.load(json_file)
        self.control = params['control']
        self.files = params['files']
        self.ions = params['params']['ions']
        self.trap = params['params']['trap']
        self.laser = params['params']['laser']
        self.stochastic = params['params']['stochastic']
        self.ccd = params['params']['ccd']

    def to_cxx_format(self):
        """
        Return a Params representation of the simulation settings
        which can be passed to the C++ module.

        Possible bug: Might need to explicitly tell np.array the
        dtypes when creating pointers.

        """
        m = np.array(self.ions['m'])*_AMU
        Z = np.array(self.ions['Z'])*_ECHARGE
        lc = np.array(self.ions['lc'], dtype=ctypes.c_int)
        m_p = m.ctypes.data_as(_DOUBLE_P)
        Z_p = Z.ctypes.data_as(_DOUBLE_P)
        masses = np.unique(self.ions['m'])
        masses_p = masses.ctypes.data_as(_DOUBLE_P)
        lc_p = lc.ctypes.data_as(_INT_P)
        khat_p = np.array(self.laser['khat']).ctypes.data_as(_DOUBLE_P)
        
        self.cxx_params = Params(
            N=len(self.ions['m']),
            N_masses=len(masses),
            m=m_p, Z=Z_p, masses=masses_p, lc=lc_p,
            khat=khat_p,
            lmbda=self.laser['lmbda'], r_l=0.,
            delta=2*pi*self.laser['delta_Gamma']*self.laser['Gamma_Hz'],
            s=self.laser['s'],
            Gamma=2*pi*self.laser['Gamma_Hz'],
            beta=self.laser['beta'],
            F0=self.laser['F0'],
            r0=self.trap['r0'], z0=self.trap['z0'],
            Omega=2*pi*self.trap['f_RF'],
            V=self.trap['V'], U=self.trap['U'], UEC=self.trap['UEC'],
            kappa=self.trap['kappa'],
            Vsec=0., w=0., # TODO
            gamma_col=self.stochastic['gamma_col'],
            m_gas=self.stochastic['m_gas_amu']*_AMU,
            T_gas=self.stochastic['T_gas'],
            sim_ccd=self.control['sim_ccd'],
            ccd_bins=self.ccd['ccd_bins'],
            ccd_extent=self.ccd['ccd_extent'],
            dt=self.control['dt'], t_max=self.control['t_max'],
            min_time=self.control['min_time'],
            abort_bounds=self.trap['r0'],
            t_steps=len(
                np.arange(0, self.control['t_max'], self.control['dt'])),
            use_rfmm=self.control['use_rfmm'],
            use_coulomb=self.control['use_coulomb'],
            use_laser=self.control['use_laser'],
            use_secular=False,
            use_stochastic=self.control['use_stochastic'],
            use_abort=True,
            num_threads=self.control['num_threads'],
            quiet=0, # TODO
            com_fname=self.files['com_fname'],
            traj_fname=self.files['traj_fname'],
            fpos_fname=self.files['fpos_fname'],
            ccd_fname=self.files['ccd_fname'],
            temp_fname=self.files['temp_fname'],
            record_traj=1, # TODO
            traj_start=self.control['dt']*100, # TODO
            T_steps=1200) # TODO
        return self.cxx_params

if __name__ == "__main__":
    params = SimParams()
    params.set_ions([40, 42], [1, 1], [1, 0])
    params.set_ccd(768, 1024)
    params.save("test.json")
    params.load("default.json")
    print(params.to_cxx_format())
