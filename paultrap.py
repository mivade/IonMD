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
import scipy.constants as consts

class PaulTrap:
    """Utility class for calculating Mathieu stability parameters."""
    def __init__(self, radius, frequency):
        """
        Initialize PaulTrap object with radius :math:`r_0` and
        frequency :math:`f` in SI units.

        """
        self.r0 = radius
        self.Omega = 2*np.pi*frequency

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
        Z : float, optional
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
        Z : float, optional
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
        Z : float, optional
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

if __name__ == "__main__":
    trap = PaulTrap(3.5e-3, 4.7e6)
    print(trap.get_mathieu_qa(40*consts.u, 1000, 0, 1))
