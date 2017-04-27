import os
import time
import numpy as np
import matplotlib.pyplot as plt
from ionmd import Params, Trap, Simulation, Status

params = Params()
trap = Trap()
sim = Simulation()

filename = "trajectories.bin"
try:
    os.remove(filename)
except OSError:
    pass
params.filename = filename

sim.set_params(params)
sim.set_trap(trap)

print("Adding ions...")

n_ions = 2

z = np.linspace(-20, 20, n_ions)

for n in range(n_ions):
    sim.add_ion(40, 1, [0, 0, z[n]])

print("Running simulation...")
sim.start()
if sim.status != Status.FINISHED:
    print("waiting...")
    time.sleep(0.5)
print(sim.status)

shape = (n_ions*3, params.num_steps)
# shape = (params.num_steps, n_ions*3)
data = np.fromfile(params.filename, dtype=np.double).reshape(shape).T

fig, ax = plt.subplots(3, n_ions)
t = np.linspace(0, params.num_steps * params.dt, params.num_steps)
lim = 200
for n in range(n_ions):
    for k in range(3):
        ax[k, n].plot(t[:lim], data[:lim, k])
plt.show()
