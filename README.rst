=====
IonMD
=====

Overview
========

IonMD is a molecular dynamics (MD) simulation of ions in a linear Paul
trap. The implmentation uses the leapfrog integration technique and
Newtonian equations of motion to clasically simulate the forces on
individual ions. Methods are largely based on the MD simulation
described in [1]_ and [2]_.

Model Description
-----------------

IonMD simulates ions in a linear Paul trap by solving each ion's
classical equation of motion :math:`\ddot{\vec{x}}_i = \vec{F}_i/m_i`
where the force on each ion is the sum of forces due to the trapping
potential, cooling lasers, Coulomb repulsion, and stochastic processes
such as background gas collisions. In other words,

.. math::
   
   \vec{F}_i = \vec{F}_{T,i} + \vec{F}_{L,i} + \vec{F}_{C,i} + \vec{F}_{S,i}

Minimization
------------

For obtaining good simulation results, it is best to start with good
initial conditions, otherwise the simulation would need to run for an
exorbitant amount of time. One method is to minimize the potential
energy :math:`U` of the system. For the :math:`i`\th particle, the
potential energy is given by 

.. math:: U_i = U_{T,i} + U_{L,i} + U_{C,i},

where we ignore contributions of stochastic effects. Specifically,
there is the trapping potential energy

.. math::

   U_{T,i} &= Ze\psi \\
           &= Ze \left[ (C-D)x_i^2 + (C-D)y_i^2 + 2Dz_i^2 \right]

where :math:`C = ZeV^2/m\Omega^2r_0^4` and :math:`D = \kappa
U_{EC}/2z_0^2`;

.. math::

   U_{L,i} = -|\vec{F}_0|(\vec{r}_i \cdot \hat{r}_L) + \frac{1}{2}
   \beta |\vec{r}_i|^2

is the potential energy contribution of the laser; and the Coulomb
potential energy is

.. math::

   U_{C,i} = \frac{Ze}{4\pi\epsilon_0} \sum_{i\neq j}
   \frac{1}{|\vec{r_i}-\vec{r_j}|}.

Finding the positions for minimization is a complex problem that
requires finding the minimum of a scalar function of :math:`3N`
variables. Presently, this is *not* the method used for
minimization. This is mainly because of difficulties in getting NLopt
to give sensible results (which might be due to an error in the
implementation).

An alternative (and the currently implemented minimization method) is
to treat *all* ions as laser cooled in the beginning and let the
system evolve for a short period of time. This has the elegance of not
requiring any extra code, but seems to suffer from not allowing
isotope sorting.

License
=======

IonMD is freely distributable under the terms of the GNU GPL version 3
(see LICENSE for details).

System Requirements
===================

IonMD uses Python as a frontend for accessing compiled C++
routines. IonMD is built and tested on Linux (Debian Wheezy), but
should work on any platform that can meet the following dependencies:

* Python requirements

  * Python 2.6 or 2.7
  * Numpy
  * Scipy
  * MayaVI (optional; for displaying ions in 3D)
  * Matplotlib (optional; for plotting/simulating CCD images)
  * PIL (optional; for simulating CCD images)

* C++ requirements

  * `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_ (for
    random number generation and 2D binning for simulated CCD images)
  * `NLopt <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_ (for
    finding the potential energy minimum for good initial conditions)

Usage
=====

The simulation at present is run through the runsim.py
file. Simulation parameters may be entered in the main() function, or
any other valid Python code can be written by the user to interface
with the C++ library.

If modifying the C++ code, the included Makefile should be sufficient
if using a standard Linux environment with g++.

Data Output
===========

Large data files containing, e.g., trajectories are written in low
level binary format. This means that they may not be portable to
another computer, but within a machine, there should be no
issues. Data can be easily read with numpy.fromfile. See the simCCD
function in runsim.py for an example.

Simulated CCD data gives mass in amu and position data in
microns. Trajectory data gives time in microseconds and position in
mm. In both cases, floats are used rather than the internal doubles in
order to reduce file size (precision should not be an issue in this
case). A "final positions" file gives x y z positions of all ions at
the end of the simulation in the plaintext xyz format (a common
chemical file format).

Future Features
===============

Some things that are planned (or possible) for the future:

* Better stability -> less random initial conditions could help this
* Converting from raw double arrays to either std::vector or possibly
  `Armadillo <http://arma.sourceforge.net/>`_ (mostly for code
  readability more than anything else)

Known Bugs
==========

* The params.py generation has an error for 64 bit systems which
  causes problems with CCD simulation. Workaround: either turn off CCD
  simulation or edit the file to make it read "Params._pack_ = 8".
* The Python function for setting up the number of ions of various
  masses doesn't currently work for more than 2 distinct masses.

Authors
=======

IonMD is principally written by Michael V. DePalatis <mvd@gatech.edu>
with some optimization enhancements by Ben Land.

References
==========

.. [1] C.B. Zhang *et al.*, Phys. Rev. A **76**, 012719 (2007).
.. [2] C.B. Zhang, *Production and Sympathetic Cooling of Complex
       Molecular Ions*, PhD thesis, Heinrich-Heine-Universität
       Düsseldorf (2008).

.. |Ba+| replace:: Ba\ :sup:`+`\ 
