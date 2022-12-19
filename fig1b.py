# Insert (2, more) spines /into/ dendrite

from neuron import h, rxd
from neuron.units import um, ms, mV, uM, mM, mm
import plotly.graph_objects as go
import plotly.express as px
import numpy as np
import math
from matplotlib import pyplot
import os

try:
    os.makedirs("concentration_plots")
except:
    ...

h.load_file("stdrun.hoc")
# model1
""" Enter angles at which to insert spines """
angle_1 = math.pi / 6  # angle (relative to dendrite/axis) of 1st spine neck
angle_2 = math.pi / 3

# model2
angle_3 = math.pi  # these remain as is (as 180 def difference) - not to be changed atm
angle_4 = 0

""" enter information here """
""" note: following information is to be constant across both models. 1->3, 2->4"""
dend_l = 4  # dendrite length
dend_diam = 1.5
location = 0.5  # x value of where to locate spinesjobs

s1_l = 3  # length of spine1 neck - must be longer than <dendrite depth>
s2_l = 3
s1_diam = 0.12 * um  # spine neck diameters
s2_diam = 0.12 * um
s1_head_diam = 0.6 * um  # spine head diameters
s2_head_diam = 0.6 * um
insert_where_1 = 0.3  # i.e. dend(insert_where_1); (e.g. dend(0.3))
insert_where_2 = 0.3

d = 0.01 * (um ** 2 / ms)  # diffusion constant
dx = 0.05  # change this accordingly, smaller = better resolution.

""" start of plotting sim """
h.load_file('stdrun.hoc')
x = location

#   model 1
# set up cylinder section
dend = h.Section(name="dend")
dend.L = dend_l * um
dend.diam = dend_diam * um
dend.nseg = 11

sp1 = h.Section(name="spine1")
sp1.L = s1_l
sp1.diam = s1_diam

sp2 = h.Section(name="spine2")
sp2.L = s2_l
sp2.diam = s2_diam

#   model 2
dend2 = h.Section(name="dend2")
dend2.L = dend_l * um
dend2.diam = dend_diam * um
dend2.nseg = 11

sp3 = h.Section(name="spine3")
sp3.L = s1_l
sp3.diam = s1_diam

sp4 = h.Section(name="spine4")
sp4.L = s2_l
sp4.diam = s2_diam

h.define_shape()


def plot_it(title, species_to_plot, region, eye=None):
    data = species_to_plot.nodes.value_to_grid()

    grid = region._mesh_grid
    xx = [grid["xlo"] + i * grid["dx"] for i in range(data.shape[0])]

    yy = [grid["ylo"] + i * grid["dy"] for i in range(data.shape[1])]

    zz = [grid["zlo"] + i * grid["dz"] for i in range(data.shape[2])]

    xs, ys, zs = np.meshgrid(
        xx,
        yy,
        zz,
        indexing="ij",
    )
    fig = go.Figure(
        data=go.Volume(
            x=xs.flatten(),
            y=ys.flatten(),
            z=zs.flatten(),
            value=np.nan_to_num(data.flatten(), nan=-1),
            isomin=-0.5,
            isomax=2,
            opacity=0.1,
            surface_count=20,
            colorscale=px.colors.sequential.Viridis
        )
    )
    fig.update_layout(title=title)
    if eye:
        fig.update_layout(scene_camera={"eye": eye})

    fig.write_image("fig1B-3D-dendrite-w-spines.svg")
    fig.write_image("fig1B-3D-dendrite-w-spines.pdf")
    fig.write_image("fig1B-3D-dendrite-w-spines.png")
    fig.show()


def find_min_depth(dend_radius, spine_radius):
    d = dend_radius - np.sqrt(dend_radius ** 2 - spine_radius ** 2)
    return d


def get_A(x, radius, depth, theta):
    """ Get first end of spine neck """
    return [x, (radius - depth) * np.cos(theta), (radius - depth) * np.sin(theta)]


def get_B(x, radius, theta):
    """ Get tangent point (on dendrite surface) of connection with spine neck """
    return [x, radius * np.cos(theta), radius * np.sin(theta)]


def get_C(x, spine_length, theta):
    """ Get end point of spine neck (outside dendrite) """
    return [x, spine_length * np.cos(theta), spine_length * np.sin(theta)]


def insert_spine_neck(spine, angle, dendrite, where):
    depth = find_min_depth(dendrite.diam / 2, spine.diam / 2)
    spine_a = get_A(x, spine.diam / 2, depth, angle)
    spine_b = get_B(x, spine.diam / 2, angle)
    spine_c = get_C(x, spine.L, angle)

    spine.pt3dclear()
    h.pt3dstyle(1, spine_b[0], spine_b[1], spine_b[2], sec=spine)
    spine.pt3dadd(spine_a[0], spine_a[1], spine_a[2], spine.diam)
    spine.pt3dadd(spine_c[0], spine_c[1], spine_c[2], spine.diam)
    spine.connect(dendrite(where))

    return spine


def insert_spine_head(spine, head_diam, angle):
    """ ATTN:   NOT finalized """
    head = h.Section(name=f"head_{str(spine)}")
    head.diam = head_diam
    head.L = 0.5
    head.pt3dclear()
    h.pt3dstyle(1, spine.x3d(1), spine.y3d(1), spine.z3d(1), sec=head)
    head_a = get_C(spine.x3d(1), spine.L, angle)
    head_c = get_C(spine.x3d(1), head.L + spine.L, angle)
    head.pt3dadd(head_a[0], head_a[1], head_a[2], head_diam)
    head.pt3dadd(head_c[0], head_c[1], head_c[2], head_diam)
    head.connect(spine(1))
    return head


def initialize(node):
    return 0 * mM if node in dend else 2 * mM


def initialize2(node):
    return 0 * mM if node in dend2 else 2 * mM


# model1
sp1 = insert_spine_neck(sp1, angle_1, dend, insert_where_1)
head1 = insert_spine_head(sp1, s1_head_diam, angle_1)

sp2 = insert_spine_neck(sp2, angle_2, dend, insert_where_2)
head2 = insert_spine_head(sp2, s2_head_diam, angle_2)

# model2
sp3 = insert_spine_neck(sp3, angle_3, dend2, insert_where_1)
head3 = insert_spine_head(sp3, s1_head_diam, angle_3)

sp4 = insert_spine_neck(sp4, angle_4, dend2, insert_where_2)
head4 = insert_spine_head(sp4, s2_head_diam, angle_4)

rxd.set_solve_type(dimension=3)

cyt = rxd.Region(dend.wholetree(), dx=dx * um)
c = rxd.Species(cyt, d=d)
c.initial = initialize

cyt2 = rxd.Region(dend2.wholetree(), dx=dx * um)
c2 = rxd.Species(cyt2, d=d)
c2.initial = initialize2

h.finitialize(-65 * mV)

# if just plotting 3d figure: uncomment this
print("Plotting 3D figure...")
plot_it("3D rendition of dendrite & spines", c, cyt)


import plotly.express as px
b_time_above_threshold = 0
r_time_above_threshold = 0


def get_time_above_threshold(data, threshold, times):
    threshold_crossed_at = np.diff(data > threshold, prepend=False)
    t = [times[x] for x in np.argwhere(threshold_crossed_at)[:,0]]
    if len(t) == 2:
        return t[1] - t[0]
    else:
        raise TypeError


def plot_max_concs(species, label, t=50 * ms, color='blue'):
    node_vectors = [h.Vector().record(node._ref_concentration) for node in species]
    times = h.Vector().record(h._ref_t)
    h.finitialize(-65 * mV)
    h.continuerun(t)
    ts = [element for element in times]
    blueplot = []
    redplot = []
    bp = []
    rp =[]
    maxes = [max(values) for values in zip(*node_vectors)]
    data = np.array(maxes)
    threshold = 0.15
    pyplot.plot(ts, maxes, color=color, label=label)
    pyplot.axhline(y=0.15, color='black', linestyle='--')

    pyplot.title(f"Maximum concentration in dendrite with 2 spines, d={d}")
    pyplot.xlabel("t [ms]")
    pyplot.ylabel("concentration in dendrite [mM]")
    pyplot.legend()
    # filename ='300'
    pyplot.savefig(
        f'concentration_plots/max_conc_plot_t_{t}_d_{d}_angles_{round(math.degrees(angle_1))}_{round(math.degrees(angle_2))}.svg')
    pyplot.savefig(
        f'concentration_plots/max_conc_plot_t_{t}_d_{d}_angles_{round(math.degrees(angle_1))}_{round(math.degrees(angle_2))}.pdf')

    if color == 'blue':
        b_time_above_threshold = get_time_above_threshold(data, threshold, ts)
        print(f"blue time above threshold: {b_time_above_threshold}")
        print(f"overall max blue: {max(maxes)}")
        print(f"overall min blue: {min(maxes[1:])}")
        print(f"final value of blue: {maxes[-1]}")
        print(f"Rise from minimum of blue: {(maxes[-1] / min(maxes[1:])) * 100}%")
        return b_time_above_threshold
    else:
        r_time_above_threshold = get_time_above_threshold(data, threshold, ts)
        print(f"red time above threshold: {r_time_above_threshold}")
        print(f"maxes: {maxes[:30]}")
        print(f"Overall max red: {max(maxes)}")
        print(f"overall min red: {min(maxes[1:])}")
        print(f"final value of red: {maxes[-1]}")

        print(f"Rise from minimum of red: {(maxes[-1] / min(maxes[1:])) * 100}%")
        return r_time_above_threshold


# if diffusion:
# print("Calculating max values...")
# r_time_above_threshold = plot_max_concs(t=2500 * ms, color='red', species=c.nodes(dend),
#                label=f'angle difference = {round(abs(math.degrees(angle_1) - math.degrees(angle_2)))}')
# b_time_above_threshold = plot_max_concs(t=2500 * ms, color='blue', species=c2.nodes(dend2),
#                label=f'angle difference = {round(abs(math.degrees(angle_3) - math.degrees(angle_4)))}')
#
# print(f"post function: red: {r_time_above_threshold}, blue: {b_time_above_threshold}")
# print(f"Difference in time above threshold: {abs(r_time_above_threshold - b_time_above_threshold)}")
# pyplot.show()
#
# print(f"percent difference in max : red is {(0.4743593 / 0.23633265) * 100} % of blue")
#