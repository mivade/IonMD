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
exorbitant amount of time. Essentially, this means minimizing the
potential energy of the system. This is currently implemented by first
treating *all* ions as laser cooled for a period of time in order to
obtain good starting positions.

License
=======

IonMD is freely distributable under the terms of the GNU GPL version 3
(see COPYING for details).

System Requirements
===================

IonMD uses Python as a frontend for accessing compiled C++
routines. IonMD is built and tested on Linux (Debian Wheezy), but
should work on any platform that can meet the following dependencies:

* Python requirements

  * Python 2.6 or 2.7
  * ctypeslib_
  * Numpy_
  * Scipy_
  * MayaVI_ (optional; for displaying ions in 3D)
  * Matplotlib_ (optional; for plotting/simulating CCD images; tested
    successfully with version 1.2.1, fails with Debian wheezy's
    version)
  * PIL_ (optional; for simulating CCD images)

.. _ctypeslib: https://pypi.python.org/pypi/ctypeslib/
.. _Numpy: http://www.numpy.org/
.. _Scipy: http://www.scipy.org/
.. _MayaVI: http://code.enthought.com/projects/mayavi/
.. _Matplotlib: http://matplotlib.org/
.. _PIL: http://www.pythonware.com/products/pil/

* C++ requirements

  * `GNU Scientific Library <https://www.gnu.org/software/gsl/>`_ (for
    random number generation and 2D binning for simulated CCD images)
  * Armadillo_ and its dependencies (on Debian wheezy, this requires
    manually adding ``libboost-math-dev`` since it is not listed as a
    dependency).

..  * `NLopt <http://ab-initio.mit.edu/wiki/index.php/NLopt>`_ (for
    finding the potential energy minimum for good initial
    conditions). Note that in the current implementation, this is not
    actually used, but the code still exists for it. I plan to later
    try using this along with nearest neighbor Coulomb potential
    approximation which is why it's staying in.

.. _Armadillo: http://arma.sourceforge.net/

Windows
-------

In principle, this should work under Cygwin_. That said, I have been
unable to get it to work due to some linking errors. An alternative is
to use a virtual machine such as VirtualBox_ with Linux installed on
it. Since modern CPUs support hardware virtualization, this should
work sufficiently well (especially if you enable multiple CPUs for the
VM).

.. _Cygwin: http://cygwin.com/
.. _VirtualBox: https://www.virtualbox.org/

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
* Python frontend cleanup

In general, see `issues
<https://github.com/mivade/IonMD/issues?state=open>`_.

Known Bugs
==========

* Temperature calculation is... wrong.

Authors
=======

IonMD is principally written by Michael V. DePalatis
<depalatis@phys.au.dk> with some optimization enhancements by Ben
Land.

References
==========

.. [1] C.B. Zhang *et al.*, Phys. Rev. A **76**, 012719 (2007).
.. [2] C.B. Zhang, *Production and Sympathetic Cooling of Complex
       Molecular Ions*, PhD thesis, Heinrich-Heine-Universität
       Düsseldorf (2008).

.. |Ba+| replace:: Ba\ :sup:`+`\ 
