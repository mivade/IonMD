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
from numpy import pi
from params import Params

class SimParams(object):
    """
    Pythonic representation of all simulation parameters. This allows
    for reading and writing of simulation setup to JSON files for easy
    storage of parameters used and helps clean up the Python front end
    to the MD simulations.

    """

    def __init__(self, init_file="default.json"):
        """
        Create the SimParams object using a default settings file if
        none is specified. Using a default file is helpful so that you
        can vary some parameters, but keep the rest the same by
        defining them only once.

        """
        self.load(init_file)

    def set_control(self, use_rfmm=False, use_coulomb=True,
                    use_laser=True, use_stochastic=False,
                    sim_ccd=True, plot_trajectory=False,
                    plot_fourier=False, display=True):
        """
        Set simulation control parameters.

        Parameters
        ----------
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
        plot_trajectory : bool
            Plot the COM trajectory after the simulation is done if
            True.
        plot_fourier : bool
            Plot the Fourier transform of the COM when the simulation
            is done if True.
        display : bool
            Use MayaVI to display the final ion positions if True.

        """

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
        self.ions['m'] = m
        self.ions['Z'] = Z
        self.ions['lc'] = lc

    def set_trap(self, r0, z0, kappa, f_RF,
                 V, U, UEC):
        """
        Set trap parameters.

        Parameters
        ----------
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
        self.trap['r0'] = r0
        self.trap['z0'] = z0
        self.trap['kappa'] = kappa
        self.trap['f_RF'] = f_RF
        self.trap['V'] = V
        self.trap['U'] = U
        self.trap['UEC'] = UEC

    def set_laser(self, khat, lmbda, Gamma_Hz, delta_Gamma,
                  s, beta, F0):
        """
        Set laser parameters.

        Parameters
        ----------
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
        self.laser['khat'] = khat
        self.laser['lmbda'] = lmbda
        self.laser['Gamma_Hz'] = Gamma_Hz
        self.laser['delta_Gamma'] = delta_gamma
        self.laser['s'] = s
        self.laser['beta'] = beta
        self.laser['F0'] = F0

    def set_stochastic(self, gamma_col, m_gas_amu, T_gas):
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

    def set_ccd(self, ccd_bins, ccd_extent):
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
        self.ions = params['params']['ions']
        self.trap = params['params']['trap']
        self.laser = params['params']['laser']
        self.stochastic = params['params']['stochastic']
        self.ccd = params['params']['ccd']

    def to_cxx_format(self):
        """
        Return a Params representation of the simulation settings
        which can be passed to the C++ module.

        """

if __name__ == "__main__":
    params = SimParams()
    params.set_ions([40, 42], [1, 1], [1, 0])
    params.set_ccd(768, 1024)
    params.save("test.json")
    params.load("default.json")
