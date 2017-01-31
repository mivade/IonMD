=====
IonMD
=====

.. note::

   The current master branch is undergoing a large overhaul to modernize the C++
   code and doesn't actually do much of anything yet. See the ``thesis`` tag for
   a version that nominally worked before.


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
exorbitant amount of time. Essentially, this means minimizing the
potential energy of the system. This is currently implemented by first
treating *all* ions as laser cooled for a period of time in order to
obtain good starting positions.


System Requirements
===================

FIXME: Installation of requirements with conda

FIXME: Python requirements

* C++ requirements

  * Armadillo_ and its dependencies (on Debian wheezy, this requires
    manually adding ``libboost-math-dev`` since it is not listed as a
    dependency).

  * Boost_

.. _Armadillo: http://arma.sourceforge.net/
.. _Boost: http://www.boost.org/


Building
========

The C++ components are built with CMake::

  $ cd build
  $ cmake ..
  $ cmake --build .


Usage
=====

TODO


Data Output
===========

The following is old:

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


Known Bugs
==========

* Temperature calculation is... wrong.


Authors
=======

IonMD is principally written by Michael V. DePalatis <depalatis@phys.au.dk> with
some optimization enhancements by Ben Land.


License
=======

IonMD is freely distributable under the terms of the GNU GPL version 3
(see COPYING for details).


References
==========

.. [1] C.B. Zhang *et al.*, Phys. Rev. A **76**, 012719 (2007).
.. [2] C.B. Zhang, *Production and Sympathetic Cooling of Complex
       Molecular Ions*, PhD thesis, Heinrich-Heine-Universität
       Düsseldorf (2008).

.. |Ba+| replace:: Ba\ :sup:`+`\
