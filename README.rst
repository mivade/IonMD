IonMD
=====

.. note::

   The current master branch is undergoing a large overhaul to modernize the C++
   code and doesn't actually do much of anything yet. See the ``thesis`` tag for
   a version that nominally worked before.


Overview
--------

IonMD is a molecular dynamics (MD) simulation of ions in a linear Paul
trap. The implmentation uses the leapfrog integration technique and
Newtonian equations of motion to clasically simulate the forces on
individual ions. Methods are largely based on the MD simulation
described in [1]_ and [2]_.


Model Description
^^^^^^^^^^^^^^^^^

IonMD simulates ions in a linear Paul trap by solving each ion's
classical equation of motion :math:`\ddot{\vec{x}}_i = \vec{F}_i/m_i`
where the force on each ion is the sum of forces due to the trapping
potential, cooling lasers, Coulomb repulsion, and stochastic processes
such as background gas collisions. In other words,

.. math::

   \vec{F}_i = \vec{F}_{T,i} + \vec{F}_{L,i} + \vec{F}_{C,i} + \vec{F}_{S,i}


Building
--------

.. warning:: This section is a work in progress.

To build, Armadillo_ and CMake_ are required. If using Anaconda or Miniconda,
these can be installed with::

    $ conda install cmake
    $ conda install -c conda-forge armadillo

.. note:: The Armadillo package on conda-forge does not include Windows
          binaries. Instead, download and install binaries from Armadillo
          website.

Building the Python bindings requires pybind11_ and scikit-build_. The former
is included here as a git submodule::

    $ git submodule init && git submodule update

The latter is installed with pip::

    $ pip install scikit-built

To build without Python bindings::

    $ mkdir -p build && cd build && cmake .. && cmake --build .

.. _Armadillo: http://arma.sourceforge.net/
.. _CMake: https://cmake.org/
.. _pybind11: https://pybind11.readthedocs.io/en/master/
.. _scikit-build: https://github.com/scikit-build/scikit-build


Usage
-----

See ``demo/demo.cpp``.


Authors
-------

IonMD is principally written by Michael V. DePalatis with some optimization
enhancements by Ben Land.


License
-------

IonMD is freely distributable under the terms of the GNU GPL version 3
(see COPYING for details).


References
----------

.. [1] C.B. Zhang *et al.*, Phys. Rev. A **76**, 012719 (2007).
.. [2] C.B. Zhang, *Production and Sympathetic Cooling of Complex
       Molecular Ions*, PhD thesis, Heinrich-Heine-Universität
       Düsseldorf (2008).

.. |Ba+| replace:: Ba\ :sup:`+`\
