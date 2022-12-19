from neuron import h, rxd
from neuron.units import um, ms, mV, uM
import plotly.graph_objects as go
import numpy as np
import plotnine as p9
import pandas as pd

h.load_file("stdrun.hoc")

rxd.set_solve_type(dimension=3)
dend = h.Section(name="dend")
dend.L = 5 * um
dend.diam = 2 * um

cyt = rxd.Region([dend], dx=0.05 * um)
c = rxd.Species(cyt, d=0.1 * um ** 2 / ms)

# pick a node and verify that it is on the surface
my_nodes = c.nodes((2.5, 0.975, 0.275))
assert len(my_nodes) == 1
my_node = my_nodes[0]
assert my_node.surface_area > 0

nearby_nodes = [c.nodes((my_node.x3d + dx, my_node.y3d, my_node.z3d))[0] for dx in [0.1, 0.5, 1, 2]]

r = h.RxDSyn(my_node.segment)
my_node.include_flux(r._ref_g)
ns = h.NetStim()
ns.number = 1
ns.start = 5 * ms

nc = h.NetCon(ns, r)
nc.weight[0] = 10000
nc.delay = 0

t = h.Vector().record(h._ref_t)
my_c = h.Vector().record(my_node._ref_concentration)
nearby_cs = [h.Vector().record(node._ref_concentration) for node in nearby_nodes]


def plot_it(title, species_to_plot, region, eye=None):
    return
    data = species_to_plot.nodes.value_to_grid() / uM
    grid = region._mesh_grid
    xs, ys, zs = np.meshgrid(
        [grid["xlo"] + i * grid["dx"] for i in range(data.shape[0])],
        [grid["ylo"] + i * grid["dy"] for i in range(data.shape[1])],
        [grid["zlo"] + i * grid["dz"] for i in range(data.shape[2])],
        indexing="ij",
    )
    fig = go.Figure(
        data=go.Volume(
            x=xs.flatten(),
            y=ys.flatten(),
            z=zs.flatten(),
            value=np.nan_to_num(data.flatten(), nan=-1),
            isomin=0,
            isomax=10,
            opacity=0.1,
            surface_count=7,
        )
    )
    fig.update_layout(title=title)
    if eye:
        fig.update_layout(scene_camera={"eye": eye})

    fig.show()


h.finitialize(-65 * mV)
h.continuerun(5.05 * ms)
plot_it(f"h.t = {h.t}", c, cyt)
print(c.nodes((2.575, 0.15, 0.25)).value)
h.continuerun(5.1 * ms)
plot_it(f"h.t = {h.t}", c, cyt)
h.continuerun(5.5 * ms)
plot_it(f"h.t = {h.t}", c, cyt)
plot_it(f"h.t = {h.t}", c, cyt, eye={"x": 2, "y": -0.75, "z": 0.5})
h.continuerun(6 * ms)
plot_it(f"h.t = {h.t}", c, cyt)

h.continuerun(10 * ms)

plot = p9.ggplot(
    pd.DataFrame({"t": t, "concentration": my_c / uM}), p9.aes(x="t", y="concentration")
) + p9.geom_line(size=1)

colors = ["red", "purple", "blue", "green"]
for timeseries, color in zip(nearby_cs, colors):
    plot += p9.geom_line(
        pd.DataFrame({"t": t, "concentration": timeseries / uM}), size=1, color=color
    )

print(
    plot + p9.xlim(4, 10) + p9.scale_y_log10(breaks=[10**i for i in range(-3, 4)])
)
