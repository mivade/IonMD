import numpy as np
import ctypes
import matplotlib.pyplot as plt
from mayavi import mlab

colors = [(1,0,0), (0,1,0), (0,0,1)]

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

def toRGB(data, channels):
    """Convert 2D grayscale data to 2D 3 channel 'RGB' data."""
    new_data = np.zeros((data.shape[0],data.shape[1],3))
    for i in range(3):
        new_data[:,:,i] = data*channels[i]
    return new_data

def simCCD(ccd_fname_prefix, N_ccd, bins, extents,
           clim=None, xlim=None, ylim=None, outfile=None):
    """Simulate a CCD image."""
    ccd = np.zeros((bins,bins,3))
    for i in range(N_ccd):
        hist = np.fromfile(ccd_fname_prefix + "_" + str(i+1) + ".dat")
        ranges = hist[:(bins*2+2)]
        data = hist[(bins*2+2):]
        data /= data.max()
        data.shape = (bins,bins)
        ccd += toRGB(data, colors[i%len(colors)])
    plt.figure()
    plt.imshow(ccd)
    plt.clim(clim)
    plt.xticks([])
    plt.yticks([])
    #plt.xlabel(r'$x$ [$\mu$m]')
    #plt.ylabel(r'$y$ [$\mu$m]')
    plt.xlim(xlim)
    plt.ylim(ylim)
    if not outfile:
        plt.show()
    else:
        plt.savefig(outfile, bbox_inches='tight')

if __name__ == "__main__":
    display()
    
