import numpy
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
from neuron import h, rxd

h.load_file("stdrun.hoc")

soma = h.Section(name="soma")
soma.L = 5
soma.diam = 5
soma.nseg = 9

soma.pt3dclear()
soma.pt3dadd(0, 0, 0, 5)
soma.pt3dadd(5, 0, 0, 5)

dend = h.Section(name="dend")
dend.L = 5
dend.diam = 1
dend.nseg = 5


dend.connect(soma(1))
rxd.set_solve_type(dimension=3)
allsecs = soma.wholetree()
r = rxd.Region(allsecs, dx=0.05)
ca = rxd.Species(r)
h.finitialize(-65)

my_nodes = [
    node for node in ca.nodes if node.surface_area and node.y3d > 0 and node.z3d > 0 and node.x3d > 0
]
xs = [node.x3d for node in my_nodes]
ys = [node.y3d for node in my_nodes]
zs = [node.z3d for node in my_nodes]

segs = {node.segment for node in ca.nodes if node.surface_area}
print("len(segments with surface area):", len(segs))

fig = plt.figure(figsize=(5, 10))
ax = fig.add_subplot(211, projection="3d")
ax2d = fig.add_subplot(212)

for seg in itertools.chain.from_iterable(allsecs):
    seg_nodes = [node for node in my_nodes if node.segment == seg]
    xs = [node.x3d for node in seg_nodes]
    ys = [node.y3d for node in seg_nodes]
    zs = [node.z3d for node in seg_nodes]

    ax.scatter(xs, ys, zs, s=1)
    ax2d.scatter(xs, ys, s=1)
    print("%-20r %20g" % (seg, sum(node.volume for node in ca.nodes(seg))))


ax.set_xlim(0, 10)
ax.set_ylim(-5, 5)
ax.set_zlim(-5, 5)

ax2d.vlines(numpy.arange(0, 5, 5 / soma.nseg), -5, 5)
ax2d.vlines(numpy.arange(5, 10, 5 / dend.nseg), -5, 5)
ax2d.plot([2.5, 6], [12.5, -1.5], color="gray")
ax2d.set_xlim(0, 10)
ax2d.set_ylim(-5, 5)


plt.show()
