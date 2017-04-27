import numpy as np
import matplotlib.pyplot as plt
from ionmd import Params, Trap, Simulation

params = Params()
trap = Trap()
sim = Simulation()

# params.filename = "trajectories.bin"

sim.set_params(params)
sim.set_trap(trap)

print("Adding ions...")

n_ions = 2

z = np.linspace(-20, 20, n_ions)

for n in range(n_ions):
    sim.add_ion(40, 1, [0, 0, z[n]])

print("Running simulation...")
sim.run()
print(sim.status)

shape = (n_ions*3, params.num_steps)
data = np.fromfile(params.filename, dtype=np.double).reshape(shape).T

plt.figure()
for n in range(n_ions):
    # plt.plot(data[n+2,:100], label=str(n + 1))
    plt.plot(data[:100,n+2], label=str(n + 1))
plt.legend()
plt.show()
