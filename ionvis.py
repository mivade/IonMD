import numpy as np
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

if __name__ == "__main__":
    display()
    
