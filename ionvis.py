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
import numpy.fft as npfft
import scipy.misc
import scipy.constants as consts
import ctypes
import matplotlib.pyplot as plt
from mayavi import mlab
try:
    from PIL import Image, ImageEnhance
except ImportError:
    import Image
    import ImageEnhance

# TODO: Conversion to use settings.SimParams instead of passing
#       arguments one by one

# color order = red, green, blue, cyan, gold
colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (1, .84, 0)]

# Utility functions
# -----------------

def to_rgb(data, channels):
    """Convert 2D grayscale data to 2D 3 channel 'RGB' data."""
    new_data = np.zeros((data.shape[0], data.shape[1], 3))
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
    plt.plot(npfft.fft(x[:end]), '-', label='$x$')
    plt.plot(npfft.fft(y[:end]), '-', label='$y$')
    plt.hold(False)
    plt.legend()
    plt.hold(False)
    plt.subplot(2, 1, 2)
    plt.plot(npfft.fft(z[:end]), '-')
    if not outfile:
        plt.show()

def plot_temperature(N, m, temp_file="temperature.txt", outfile=None):
    """
    Plot the temperature over time.

    Parameters
    ----------
    N : int
        Number of ions.
    m : float
        Average mass of ions in SI units.
    temp_file : str, optional
        Path to file storing temperature data.
    outfile : str or None
        If not None, a string containing the filename to save the plot
        to.

    """
    t, v = np.loadtxt(temp_file, unpack=True)
    t /= 1e-6
    T = m*v**2/(3*N*consts.k)
    plt.figure()
    plt.plot(t, T)
    plt.xlabel(r'$t$ [$\mu$s]')
    plt.ylabel(r'T [???]')
    if not outfile:
        plt.show()

# 3D display functions
# --------------------

def display(pos_fname='fpos.xyz', scale=25, outfile=None):
    """
    Display the 3D positions of ions given by a file in xyz format.

    Parameters
    ----------
    pos_fname : string, optional
        Path to xyz file to use. Defaults to fpos.xyz.
    scale : float, optional
        Scaling factor for MayaVI glyphs.
    outfile : string, optional
        If specified, save the resulting visualization to this
        filename instead of displaying the results.

    """
    scale = scale
    ions, x, y, z = np.loadtxt(pos_fname, skiprows=2, unpack=True)
    xlist = []
    m_last, ci = ions[0], 0
    fig = mlab.figure(size=(640, 480), bgcolor=(0, 0, 0))
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

def sim_ccd(ccd_file, N_ccd, bins, extents,
            outfile=None, show=False, brightness=1):
    """
    Simulate a CCD image. If ccd_file is a string, opens N_ccd files
    named ccd_file + '_<integer>.dat'. If a length N_ccd list of file
    objects, it reads those files instead.

    Parameters
    ----------
    ccd_file : string or list
        If a string, the prefix for the CCD file names. If a list or
        tuple, the list of all absolute file names for CCD images.
    N_ccd : int
        Number of CCD files to load and process.
    bins : int
        Bins contained in the CCD histograms.
    extents : float
        Physical extents of the CCD data.
    outfile : string, optional
        Filename to save the output to. If None, the image will not be
        saved.
    show : bool, optional
        If True, display the image.
    brightness : float, optional
        Factor to enhance the image brightness by.

    Returns
    -------
    ccd : numpy.ndarray
        A numpy array representation of the CCD data.

    """
    tmp_img = "images/tmp.png"
    ccd = np.zeros((bins, bins, 3))
    for i in range(N_ccd):
        if isinstance(ccd_file, (str, unicode)):
            hist = np.fromfile(ccd_file + "_" + str(i + 1) + ".dat")
        elif isinstance(ccd_file, list):
            hist = np.fromstring(ccd_file[i])
        else:
            raise TypeError("ccd_file must be a string or a list of " + \
                            "file objects.")
        ranges = hist[:(bins*2 + 2)]
        data = hist[(bins*2 + 2):]
        data /= data.max()
        data.shape = (bins, bins)
        ccd += to_rgb(data, colors[i%len(colors)])
    scipy.misc.imsave(tmp_img, ccd)
    img_a = Image.open(tmp_img)
    img_b = ImageEnhance.Brightness(img_a)
    img_c = img_b.enhance(brightness)
    if outfile is not None:
        img_c.save(outfile)
    if show:
        img_c.show()
    return ccd

def sim_ccd_archive(tar_fname, N_ccd, bins, extents, **kwargs):
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

    TODO: add kwargs for things to go to sim_ccd

    """
    tar = tarfile.open(tar_fname, 'r')
    names = tar.getnames()
    prefixes = sorted(set([x.split('_')[0] for x in names]))
    for prefix in prefixes:
        ccd_file = []
        for i in range(N_ccd):
            data = tar.extractfile(prefix + "_%i.dat" % (i+1)).read()
            ccd_file.append(data)
        sim_ccd(ccd_file, N_ccd, bins, extents,
               outfile=(os.path.dirname(tar_fname)+"/"+prefix+'.png'))
    tar.close()

if __name__ == "__main__":
    #ccd_bins, ccd_extent = 512, 600
    #display()
    #sim_ccd("ccd", 2, ccd_bins, ccd_extent, brightness=1.5,
    #       outfile="images/CCD_latest.png", show=True)
    #bins, extents = 768, 1024
    #sim_ccd_archive("data/N8/ccd/ccddat.tar.gz", 2, bins, extents)
    import argparse
    parser = argparse.ArgumentParser()
