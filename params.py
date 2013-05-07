from ctypes import *

STRING = c_char_p


class Params(Structure):
    pass
Params._pack_ = 4
Params._fields_ = [
    ('N', c_int),
    ('N_masses', c_int),
    ('m', POINTER(c_double)),
    ('Z', POINTER(c_double)),
    ('masses', POINTER(c_double)),
    ('lc', POINTER(c_int)),
    ('khat', POINTER(c_double)),
    ('lmbda', c_double),
    ('r_l', c_double),
    ('delta', c_double),
    ('s0', c_double),
    ('Gamma', c_double),
    ('r0', c_double),
    ('z0', c_double),
    ('kappa', c_double),
    ('Omega', c_double),
    ('V', c_double),
    ('U', c_double),
    ('UEC', c_double),
    ('Vsec', c_double),
    ('w', c_double),
    ('gamma_col', c_double),
    ('sim_ccd', c_int),
    ('ccd_bins', c_double),
    ('ccd_extent', c_double),
    ('dt', c_double),
    ('t_max', c_double),
    ('abort_bounds', c_double),
    ('minimizing', c_int),
    ('t_steps', c_int),
    ('use_rfmm', c_int),
    ('use_coulomb', c_int),
    ('use_laser', c_int),
    ('use_secular', c_int),
    ('use_stochastic', c_int),
    ('use_abort', c_int),
    ('num_threads', c_int),
    ('quiet', c_int),
    ('traj_fname', STRING),
    ('fpos_fname', STRING),
    ('ccd_fname', STRING),
    ('temp_fname', STRING),
    ('record_traj', c_int),
    ('traj_start', c_double),
    ('T_steps', c_int),
]
__all__ = ['Params']
