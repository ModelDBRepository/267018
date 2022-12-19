from neuron import h, rxd
from neuron.units import ms, mV, µm, mM
import plotnine as p9
import pandas as pd
import math
from functools import lru_cache
import time
import random

h.load_file("stdrun.hoc")

initial_concentration = 1 * mM
start_x = -2 * µm
stop_x = 2 * µm
start_y = -2 * µm
stop_y = 2 * µm
start_z = -2 * µm
stop_z = 2 * µm
tstop = 20 * ms
D = 1 * µm ** 2 / ms

def fundamental_solution(t, D, distance):
    return 1 / (4 * math.pi * D * t) ** (1.5) * math.exp(-(distance ** 2) / (4 * D * t))


@lru_cache(maxsize=None)
def solution_from_domain(
    x, y, z, t, D, interval_steps=100
):
    """solution of diffusion on an infinite space from a finite rectangular prism source"""
    assert start_x < stop_x
    dx = (stop_x - start_x) / interval_steps
    dy = (stop_y - start_y) / interval_steps
    dz = (stop_z - start_z) / interval_steps
    return initial_concentration * dx * dy * dz * sum(
        fundamental_solution(t, D, 
            math.sqrt(
                (x - (start_x + (i + 0.5) * dx)) ** 2
                + (y - (start_y + (j + 0.5) * dy)) ** 2
                + (z - (start_z + (k + 0.5) * dz)) ** 2
        )) for i in range(interval_steps) for j in range(interval_steps) for k in range(interval_steps)
    )

domain = h.Section(name="domain")
domain.pt3dadd(-20, 0, 0, 40)
domain.pt3dadd(20, 0, 0, 40)


# initialization rule
def init_concentration(node):
    return initial_concentration if start_x < node.x3d < stop_x and start_y < node.y3d < stop_y and start_z < node.z3d < stop_z else 0

rxd.set_solve_type(domain, dimension=3)
cyt = rxd.Region([domain], name="cyt", nrn_region="i", dx=0.25)
ca = rxd.Species(cyt, name="ca", d=D, charge=2, initial=init_concentration)

start = time.time()
h.finitialize(-65 * mV)
print(f"initializing took {time.time() - start} seconds")

start = time.time()
h.continuerun(tstop)
print(f"simulating took {time.time() - start} seconds")

pts_checked = set()
nodes = ca.nodes
r = []
errors = []
while len(r) < 100:
    i = random.randint(0, len(nodes) - 1)
    if i not in pts_checked:
        pts_checked.add(i)
        node = nodes[i]
        x = node.x3d
        y = node.y3d
        z = node.z3d
        if x ** 2 + y ** 2 + z ** 2 < 100:
            solution = solution_from_domain(x, y, z, tstop, D)
            r.append(math.sqrt(x ** 2 + y ** 2 + z ** 2))
            errors.append(100 * abs(solution - node.concentration) / solution)
            print(f"len(r) = {len(r)}")

p9.options.figure_size = (4, 3)
g = (
    p9.ggplot(pd.DataFrame({"r": r, "error": errors}), p9.aes(x="r", y="error")) + p9.geom_point() + 
    p9.labs(x='Distance from origin (µm)', y='Relative error (%)')
)
g.save("diffusion-3d-comparison.pdf")
