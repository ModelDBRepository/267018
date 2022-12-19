# NEURON Methods paper, Figure 1A generation code. Plot wave over time in 3D.
from neuron import h, rxd
import time
import numpy as np
from matplotlib import pyplot as plt
from neuron.rxd.node import Node3D

h.load_file('stdrun.hoc')
h.load_file('import3d.hoc')


class Cell:
    def __init__(self,filename):
        """Read geometry from a given SWC file and create a cell with a K+ source"""
        cell = h.Import3d_Neurolucida3()
        cell.input(filename)
        h.Import3d_GUI(cell, 0)
        i3d = h.Import3d_GUI(cell, 0)
        i3d.instantiate(self)
        for sec in self.all:
            sec.nseg = 1 + 10 * int(sec.L / 5)
            sec.insert('steady_k')


mycell = Cell('070314F_11.ASC')    # load cell 070314F_11.ASC from local directory

min_diam = min(min(sec.diam3d(i) for i in range(sec.n3d())) for sec in mycell.all)

secs3d = [mycell.apic[0], mycell.apic[1]] + [dend for dend in h.allsec() if h.distance(dend(0.5), mycell.soma[0](0.5)) < 70]
min_diam_3d = min(min(sec.diam3d(i) for i in range(sec.n3d())) for sec in secs3d)

print(f"The minimum diameter in the whole cell is {min_diam} µm")
print(f"The minimum diameter in the 3D part is {min_diam_3d} µm")

rxd.set_solve_type(secs3d, dimension=3)
rxd.nthread(4)
# Set nseg for our 1D sections
secs1d = [sec for sec in h.allsec() if sec not in secs3d]
for sec in secs1d:
    sec.nseg = 11


def plot_contours(species, i, perspective=1):
    r = species.nodes[0].region
    # perspective1 = xz axes
    # perspective2 = xy axes

    def replace_nans(a, b):
        if np.isnan(a):
            return b
        return max(a, b)

    if perspective == 1:
        flat = np.empty((max(r._xs)+1, max(r._zs)+1))
        flat.fill(np.nan)

    elif perspective == 2:
        flat = np.empty((max(r._xs)+1, max(r._ys)+1))
        flat.fill(np.nan)

    for node in ca.nodes:
        if isinstance(node, Node3D):
            if h.distance(node, mycell.soma[0](0.5)) >= h.distance(mycell.apic[3](0), mycell.soma[0](0.5)):
                continue
                # apic[3] is the section cutoff for this particular cell. change section by choice
            if perspective==1:
                flat[node._i, node._k] = replace_nans(flat[node._i, node._k], node.value)
            elif perspective==2:
                flat[node._i, node._j] = replace_nans(flat[node._i, node._j], node.value)

    xs, ys = np.meshgrid(range(flat.shape[1]), range(flat.shape[0]))

    plt.contour(xs, ys, np.nan_to_num(flat), [0.5], colors='k', linewidths=0.5)
    plt.axis('equal')
    plt.axis('off')


dx=0.17
r = rxd.Region(h.allsec(), nrn_region='i', dx=dx)
ca = rxd.Species(r, d= 0.25, name='ca', charge=2, initial= lambda node: 1 if node.sec in [mycell.apic[8]] else 0)
bistable_reaction = rxd.Rate(ca, -ca * (1 - ca) * (0.01 - ca))
h.dt = .115     # We choose dt = 0.1 here because the ratio of d * dt / dx**2 must be less than 1
print(f"starting initialization at {time.perf_counter()}")
h.finitialize(-65)
print(f"finished initialization at {time.perf_counter()}")

rng = 190   # number of timesteps
run = 3     # time-step length in ms
perspective = 1     # get both perspectives
for i in range(rng):
    start = time.perf_counter()
    print(f"started {i} at: {start}")
    h.continuerun(i*run)
    plt.figure(2, figsize=(15,27.6))    # this choice of size is arbitrary
    if max(ca.nodes(mycell.apic[1]).concentration) > 0.5:   # changed from 0.5 to 0 to finally get SOME results
        print("Plotting contours...")
        plt.figure(1)
        plot_contours(ca, i, perspective=perspective)   # get both perspectives
        plt.figure(2)
        plot_contours(ca, i, perspective=2)

    print(f"time for {i}: {time.perf_counter()-start}")


for i in [1,2]:
    plt.figure(i)
    plt.savefig(f"fig1a/p_{i}_Figure1A_hybrid_3d_dx_{dx}_run_{run}ms_rng{rng}.svg")
    plt.savefig(f"fig1a/p_{i}Figure1A_hybrid_3d_dx_{dx}_run_{run}ms_rng{rng}.pdf")

plt.show()
