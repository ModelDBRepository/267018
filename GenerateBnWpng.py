# import numpy as np
import png
import pickle
import math
import json
# NEURON Methods paper, Figure 1A generation code. Plot wave over time in 3D.
from neuron import h, rxd
import neuron
import time
import numpy as np
from neuron.rxd.node import Node3D
import matplotlib.pyplot as plt

neuron.rxd.options.ics_partial_volume_resolution = 1

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

mycell = Cell('asc/070314F_11.ASC')    # load cell 070314F_11.ASC from local directory

secs3d = [mycell.apic[0], mycell.apic[1]] + [dend for dend in h.allsec() if h.distance(dend(0.5), mycell.soma[0](0.5)) < 70]
rxd.set_solve_type(secs3d, dimension=3)
rxd.nthread(4)
# Set nseg for our 1D sections
secs1d = [sec for sec in h.allsec() if sec not in secs3d]
for sec in secs1d:
    sec.nseg = 11

def get_png(species, i, timestep, perspective=1):
    r = species.nodes[0].region
    # perspective1 = xz axes
    # perspective2 = xy axes

    def replace_nans(a, b):
        if np.isnan(a):
            return b
        return max(a, b)

    if perspective == 1:
        flat = np.empty((max(r._xs)+1, max(r._zs)+1, max(r._ys)+1), dtype=np.uint8)

    elif perspective == 2:
        flat = np.empty((max(r._xs)+1, max(r._ys)+1, max(r._zs)+1), dtype=np.uint8)

    for node in ca.nodes:
        if isinstance(node, Node3D):
            # if h.distance(node, mycell.soma[0](0.5)) >= h.distance(mycell.apic[3](0), mycell.soma[0](0.5)):
            #     continue
            # apic[3] is the section cutoff for this particular cell. change section by choice
            if perspective==1:
                flat[node._i, node._k, node._j] = 255
            elif perspective==2:
                flat[node._i, node._j, node._k] = 255

    with open(f"bwpng{timestep:04d}.pk", 'wb') as f:
        f.write(pickle.dumps(flat))

    print(f"Meshgrid xlo = {r._mesh_grid['xlo']}, ylo = {r._mesh_grid['ylo']}, zlo = {r._mesh_grid['zlo']}")
    # xs, ys = np.meshgrid(range(flat.shape[1]), range(flat.shape[0]))
    ct = 0
    for k in range(flat.shape[2]):
        # arr = flat[]
        png.from_array(flat[:,:,k].copy(), mode='L').save(f"bwpngno150/cell_{ct:04d}.png")
        # plt.imshow(flat[:,:,k], vmin=0, vmax=1)
        # plt.savefig(f"bwpng/Sliced_time_{timestep:04d}_z_{k:04d}.png")
        # plt.close()
        print(f"Printed png #{k}")
        ct+=1

def get_png_pk(file):
    f = pickle.load(file)
    ct=0
    for k in range(150, f.shape[2]):

        # arr = flat[]
        png.from_array(file[:,:,k].copy(), mode='L').save(f"bwpng/cell_{ct:04d}.png")
        # plt.imshow(flat[:,:,k], vmin=0, vmax=1)
        # plt.savefig(f"bwpng/Sliced_time_{timestep:04d}_z_{k:04d}.png")
        # plt.close()
        print(f"Printed png #{k}")
        ct+=1


dx=0.25
r = rxd.Region(h.allsec(), nrn_region='i', dx=dx)
# ca = rxd.Species(r, d= 0.25, name='ca', charge=2, initial= lambda node: 1 if node.sec in [mycell.apic[8]] else 0)
ca = rxd.Species(r, d=0.25, name='ca', charge=2, initial=lambda node: 1 if node.x3d > 50 else 0)

nodes3d = [node for node in ca.nodes if any(node in sec for sec in secs3d)]
print(f"#nodes is: {len(nodes3d)}")
print(f"vcell projected volume: {len(nodes3d) * (dx**3)}")
neuronvolume = [node.volume for node in nodes3d]
print(f"neuron volume with dx={dx} is: {sum(neuronvolume)}")
plt.hist(neuronvolume, bins=200)
plt.show()
# bistable_reaction = rxd.Rate(ca, -ca * (1 - ca) * (0.01 - ca))
# h.dt = .115     # We choose dt = 0.1 here because the ratio of d * dt / dx**2 must be less than 1
# print(f"starting initialization at {time.perf_counter()}")
#
# h.finitialize(-65)
#
# print(f"Meshgrid xlo = {r._mesh_grid['xlo']}, ylo = {r._mesh_grid['ylo']}, zlo = {r._mesh_grid['zlo']}")
# print(f"Meshgrid xhi = {r._mesh_grid['xhi']}, yhi = {r._mesh_grid['yhi']}, zhi = {r._mesh_grid['zhi']}")
# print("maxs: ", max(r._xs)+1, max(r._ys)+1, max(r._zs)+1)
#
# h.continuerun(225)
# print(f"finished initialization at {time.perf_counter()}")
# get_png(ca, 0, timestep=225, perspective=2)
# with open(f"bwpng{225:04d}.pk", 'rb') as f:
#     # f2list = pickle.load(f2)
#     get_png_pk(f)
    # print(f'zs: {len(sorted(set(f2list)))}')


# rng = 19   # number of timesteps
# run = 30     # time-step length in ms
# for i in range(rng):
#     start = time.perf_counter()
#     print(f"started {i} at: {start}")
#     h.continuerun(i*run)
#
#     get_png(ca, i, timestep=i*run, perspective=2)
#     print(f"time for {i}: {time.perf_counter()-start}")
