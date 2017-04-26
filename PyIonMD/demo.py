import numpy as np
from pyionmd import SimParams, Trap, Simulation

params = SimParams()
trap = Trap()
sim = Simulation()

sim.set_params(params)
sim.set_trap(trap)

print("Adding ions...")

n_ions = 20

z = np.linspace(-200, 200, n_ions)

for n in range(n_ions):
    sim.add_ion(40, 1, [0, 0, z[n]])

print("Running simulation...")
sim.run()
print(sim.status)
