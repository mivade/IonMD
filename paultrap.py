"""
paultrap.py

Defines a utility class for calculating Mathieu stability parameters.

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
import numpy as np
from numpy import pi
import scipy.constants as consts

class PaulTrap:
    """Utility class for calculating Mathieu stability parameters."""
    def __init__(self, radius, frequency, z0=None, kappa=None):
        """
        Initialize PaulTrap object.

        Parameters
        ----------
        radius : float
            Trap radius in m.
        frequency : float
            Trap RF frequency in Hz.
        z0 : float, optional
            Half the distance between the trap end caps in m. Only
            necessary for calculating secular frequencies.
        kappa : float, optional
            Trap geometrical factor. Only necessary for calculating
            secular frequencies.

        """
        self.r0 = radius
        self.Omega = 2*np.pi*frequency
        self.z0 = z0
        self.kappa = kappa

    def get_mathieu_q(self, m, V, Z=1):
        """
        Calculate and return the Mathieu stability parameter :math:`q`
        for a given ion and RF voltage.

        Parameters
        ----------
        m : float
            Ion mass in SI units.
        V : float
            RF voltage V as defined in P. K. Ghosh, *Ion Traps*
            (Oxford University Press, Oxford, 1995).
        Z : int, optional
            Ionization state (i.e., ``Z = 1`` means an ion with charge
            ``Z*e``).

        Returns
        -------
        q : float
            The calculated stability parameter :math:`q`.

        """
        q = 2*Z*consts.e*V/(m*self.r0**2*self.Omega**2)
        return q

    def get_mathieu_a(self, m, U, Z=1):
        """
        Calculate and return the Mathieu stability parameter :math:`a`
        for a given ion and RF voltage.

        Parameters
        ----------
        m : float
            Ion mass in SI units.
        U : float
            RF DC offset voltage U as defined in P. K. Ghosh, *Ion
            Traps* (Oxford University Press, Oxford, 1995).
        Z : int, optional
            Ionization state (i.e., ``Z = 1`` means an ion with charge
            ``Z*e``).

        Returns
        -------
        a : float
            The calculated stability parameter :math:`a`.

        """
        a = 4*Z*consts.e*U/(m*self.r0**2*self.Omega**2)
        return a

    def get_mathieu_qa(self, m, V, U, Z=1):
        """
        Calculate and return both the Mathieu stability parameters
        :math:`q` and :math:`a`.

        Parameters
        ----------
        m : float
            Ion mass in SI units.
        V : float
            RF voltage V as defined in P. K. Ghosh, *Ion Traps*
            (Oxford University Press, Oxford, 1995).
        U : float
            RF DC offset voltage U as defined in P. K. Ghosh, *Ion
            Traps* (Oxford University Press, Oxford, 1995).
        Z : int, optional
            Ionization state (i.e., ``Z = 1`` means an ion with charge
            ``Z*e``).

        Returns
        -------
        q : float
            The calculated stability parameter :math:`q`.
        a : float
            The calculated stability parameter :math:`a`.

        """
        q = self.get_mathieu_q(m, V, Z)
        a = self.get_mathieu_a(m, U, Z)
        return q, a

    def _get_w0(self, m, V, U, Z=1):
        """
        Calculate the nominal radial secular frequency for an ion of
        mass m and charge Z*e. That is, calculate the ideal radial
        frequency in the case of an infinitely long ion guide.

        Parameters
        ----------
        m : float
            Ion mass in SI units.
        V : float
            RF voltage V as defined in P. K. Ghosh, *Ion Traps*
            (Oxford University Press, Oxford, 1995).
        U : float
            RF DC offset voltage U as defined in P. K. Ghosh, *Ion
            Traps* (Oxford University Press, Oxford, 1995).
        Z : int, optional
            Ionization state (i.e., ``Z = 1`` means an ion with charge
            ``Z*e``).

        Returns
        -------
        float
            Nominal radial secular (angular) frequency.

        """
        q, a = self.get_mathieu_qa(m, V, U, Z)
        return np.sqrt(a + q**2/2.0)*self.Omega/2.0

    def get_axial_freq(self, m, UEC, Z=1):
        """
        Calculate the axial secular frequency for an ion of mass m and
        charge Z*e.

        Parameters
        ----------

        Returns
        -------
        float
            Axial secular (angular) frequency.

        """
        try:
            return np.sqrt(2*Z*consts.e*self.kappa*UEC/(m*self.z0**2))
        except TypeError:
            raise ValueError("z0 and kappa must be defined to " + \
                             "calculate secular frequencies.")
        

    def get_radial_freq(self, m, V, U, UEC, Z=1):
        """
        Calculate the actual radial secular frequency.

        Parameters
        ----------
        m : float
            Ion mass in SI units.
        V : float
            RF voltage V as defined in P. K. Ghosh, *Ion Traps*
            (Oxford University Press, Oxford, 1995).
        U : float
            RF DC offset voltage U as defined in P. K. Ghosh, *Ion
            Traps* (Oxford University Press, Oxford, 1995).
        UEC : float
            End cap electrode voltage.
        Z : int, optional
            Ionization state (i.e., ``Z = 1`` means an ion with charge
            ``Z*e``).

        Returns
        -------
        float
            Radial secular (angular) frequency.

        """
        w_0 = self._get_w0(m, V, U, Z)
        w_z = self.get_axial_freq(m, UEC, Z)
        return np.sqrt(w_0**2 - 0.5*w_z**2)

if __name__ == "__main__":
    trap = PaulTrap(3.5e-3, 4.7e6, 2.70e-3, 0.248)
    m = 40*consts.u
    V = 500.
    UEC = 0.27
    print(trap.get_mathieu_qa(m, V, 0, 1))
    print("w_r =", trap.get_radial_freq(m, V, 0, UEC, 1)/2/pi/1000, "kHz")
    print("w_z =", trap.get_axial_freq(m, UEC, 1)/2/pi/1000, "kHz")
