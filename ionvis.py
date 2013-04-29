import numpy as np
import ctypes
import matplotlib.pyplot as plt
from mayavi import mlab

def display(fpos_fname='fpos.xyz', m_lc=138):
    scale = 25.
    ions, x, y, z = np.loadtxt(fpos_fname, skiprows=2, unpack=True)
    c_lc = (0,1,1)
    c_sc = (1,0,0)
    x_lc = []
    x_sc = []
    for i, ion in enumerate(ions):
        if ion == m_lc:
            x_lc.append([x[i], y[i], z[i]])
        else:
            x_sc.append([x[i], y[i], z[i]])
    x_lc = np.array(x_lc)
    x_sc = np.array(x_sc)
    mlab.points3d(x_lc[:,0], x_lc[:,1], x_lc[:,2], color=c_lc, scale_factor=scale)
    try:
        mlab.points3d(x_sc[:,0], x_sc[:,1], x_sc[:,2], color=c_sc, scale_factor=scale)
    except:
        pass
    mlab.show()

def simCCD(ccd_fname="ccd.dat", bins=500, clim=None, outfile=None):
    """Simulate a CCD image."""
    tmp = np.fromfile(ccd_fname, dtype=ctypes.c_float)
    tmp.shape = (tmp.shape[0]/3, 3)
    m = tmp[:,0]
    x = tmp[:,1]
    y = tmp[:,2]
    tmp = None
    h, xe, ye = np.histogram2d(x, y, bins=bins)
    #print xe, ye
    plt.figure()
    plt.hot()
    plt.imshow(h)
    plt.clim(clim)
    #plt.xticks([])
    #plt.yticks([])
    plt.xlabel(r'$x$ [$\mu$m]')
    plt.ylabel(r'$y$ [$\mu$m]')
    if not outfile:
        plt.show()
    else:
        plt.savefig(outfile)

if __name__ == "__main__":
    display()
    
