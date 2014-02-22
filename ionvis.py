"""
ionvis.py

IonMD ion visualization functions.

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

import tarfile
import os.path
import numpy as np
from numpy.fft import fft
import scipy.constants as consts
import ctypes
import matplotlib.pyplot as plt
import Image, ImageEnhance
from mayavi import mlab

# TODO: PEP8 renaming of functions
# TODO: Conversion to use settings.SimParams instead of passing
#       arguments one by one

# Physical constants
amu = consts.u
q_e = consts.e
kB = consts.k

# color order = red, green, blue, cyan, gold
colors = [(1,0,0), (0,1,0), (0,0,1), (0,1,1), (1,.84,0)]

# Utility functions
# -----------------

def toRGB(data, channels):
    """Convert 2D grayscale data to 2D 3 channel 'RGB' data."""
    new_data = np.zeros((data.shape[0],data.shape[1],3))
    for i in range(3):
        new_data[:,:,i] = data*channels[i]
    return new_data

# Basic plotting
# --------------

def plot_trajectory(dt, t_max, N, traj_file="traj.dat",
                    outfile=None, start=0, end=1000):
    """Plot the saved ion trajectory."""
    traj = np.fromfile(traj_file, dtype=ctypes.c_float)
    t = (np.arange(0, t_max, dt)/1e-6)[start:]
    traj.shape = (len(traj)/3, 3)
    x, y, z = traj[:,0], traj[:,1], traj[:,2]
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.hold(True)
    plt.plot(t[:end], x[:end], '-', label='$x$')
    plt.plot(t[:end], y[:end], '-', label='$y$')
    plt.ylabel('$x, y$ [mm]')
    plt.legend()
    plt.hold(False)
    plt.subplot(2, 1, 2)
    plt.plot(t[:end], z[:end], '-')
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'$z$ [mm]')
    if not outfile:
        plt.show()

def plot_fourier(dt, t_max, N, traj_file="traj.dat",
                 outfile=None, start=0, end=1000):
    """
    Plot the Fourier transform of the trajectory data to extract
    motional frequencies.

    """
    traj = np.fromfile(traj_file, dtype=ctypes.c_float)
    #t = (np.arange(0, t_max, dt)/1e-6)[start:]
    traj.shape = (len(traj)/3, 3)
    x, y, z = traj[:,0], traj[:,1], traj[:,2]
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.hold(True)
    plt.plot(fft(x[:end]), '-', label='$x$')
    plt.plot(fft(y[:end]), '-', label='$y$')
    plt.hold(False)
    plt.legend()
    plt.hold(False)
    plt.subplot(2, 1, 2)
    plt.plot(fft(z[:end]), '-')
    if not outfile:
        plt.show()

def plot_temperature(N, m=138*amu, temp_file="temperature.txt", outfile=None):
    """Plot the temperature over time."""
    t, v = np.loadtxt(temp_file, unpack=True)
    t /= 1e-6
    T = m*v**2/(3*N*kB)
    plt.figure()
    plt.plot(t, T)
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'T [???]')
    if not outfile:
        plt.show()

# 3D display functions
# --------------------

def display(fpos_fname='fpos.xyz', scale=25, outfile=None):
    scale = scale
    ions, x, y, z = np.loadtxt(fpos_fname, skiprows=2, unpack=True)
    xlist = []
    m_last, ci = ions[0], 0
    fig = mlab.figure(size=(640,480), bgcolor=(0,0,0))
    for i, ion in enumerate(ions):
        if ion != m_last or i == len(ions)-1:
            xlist = np.array(xlist)
            mlab.points3d(xlist[:,0], xlist[:,1], xlist[:,2],
                          color=colors[ci], scale_factor=scale)
            xlist = []
            ci += 1
            m_last = ions[i]
        xlist.append([x[i], y[i], z[i]])
    mlab.view(azimuth=45, elevation=90)
    mlab.roll(180)
    mlab.orientation_axes()
    if outfile:
        mlab.savefig(outfile, figure=fig)
        mlab.close(all=True)
    else:
        mlab.show()

# CCD simulation
# --------------

def simCCD(ccd_file, N_ccd, bins, extents,
           outfile=None, show=False, brightness=1, imgcmd="eog"):
    """Simulate a CCD image. If ccd_file is a string, opens N_ccd
    files named ccd_file + '_<integer>.dat'. If a length N_ccd list of
    file objects, it reads those file objects instead."""
    tmpimg = "images/tmp.png"
    ccd = np.zeros((bins,bins,3))
    for i in range(N_ccd):
        if isinstance(ccd_file, str):
            hist = np.fromfile(ccd_file + "_" + str(i+1) + ".dat")
        elif isinstance(ccd_file, list):
            hist = np.fromstring(ccd_file[i])
        else:
            raise TypeError("ccd_file must be a string or a list of file objects.")
        ranges = hist[:(bins*2+2)]
        data = hist[(bins*2+2):]
        data /= data.max()
        data.shape = (bins,bins)
        ccd += toRGB(data, colors[i%len(colors)])
    plt.imsave(fname=tmpimg, arr=ccd)
    imgA = Image.open(tmpimg)
    imgB = ImageEnhance.Brightness(imgA)
    imgC = imgB.enhance(brightness)
    if outfile:
        imgC.save(outfile)
    if show:
        imgC.show(command=imgcmd)
    return ccd

def simCCDArchive(tar_fname, N_ccd, bins, extents, **kwargs):
    """
    Open a tar archive containing a bunch of simulated CCD data
    files. The data files in the archive must have a name of the
    form::

      <prefix>_<integer>.dat

    where `<prefix>` is some number or otherwise distinguising bit of
    the file name in order to differentiate between simulation runs
    and the integer following the underscore is (inclusively) between
    1 and N_ccd and indicates the CCD for a particular ion
    species. For example, an archive might contain the following
    filenames::
    
      ['fraction0.3ccd_1.dat',
       'fraction0.3ccd_2.dat',
       'fraction0.3ccd_3.dat',
       'fraction0.4ccd_1.dat',
       'fraction0.4ccd_2.dat',
       'fraction0.4ccd_3.dat]

    In this case, there are two `<prefix>` values: `fraction0.3ccd`
    and `fraction0.4ccd` which correspond to different simulation
    runs, and there being three of each indicates there were 3
    simulated CCDs for 3 different ion species.

    TODO: add kwargs for things to go to simCCD
    """
    tar = tarfile.open(tar_fname, 'r')
    names, members = tar.getnames(), tar.getmembers()
    prefixes = sorted(set([x.split('_')[0] for x in names]))
    for prefix in prefixes:
        ccd_file = []
        for i in range(N_ccd):
            data = tar.extractfile(prefix + "_%i.dat" % (i+1)).read()
            ccd_file.append(data)
        simCCD(ccd_file, N_ccd, bins, extents,
               outfile=(os.path.dirname(tar_fname)+"/"+prefix+'.png'))
    tar.close()

if __name__ == "__main__":
    #ccd_bins, ccd_extent = 512, 600
    #display()
    #simCCD("ccd", 2, ccd_bins, ccd_extent, brightness=1.5,
    #       outfile="images/CCD_latest.png", show=True)
    bins, extents = 768, 1024
    simCCDArchive("data/N8/ccd/ccddat.tar.gz", 2, bins, extents)
