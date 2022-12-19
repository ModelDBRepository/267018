from neuron import h, rxd
from neuron.units import µm, ms, mV, nM, µM
import plotnine as p9
import pandas as pd
import numpy as np

h.load_file("stdrun.hoc")

# create a Y-shaped geometry
# we're letting NEURON pick the join angle, but this can be viewed
# using an h.PlotShape
sec = [h.Section(name=f"sec[{i}]") for i in range(3)]
for child_sec in sec[1:]:
    child_sec.connect(sec[0])

# select a small enough scale to run quickly
# all sections have length 10 and diameter 2
sec[0].pt3dadd(0, 0, 0, 2)
for my_sec in sec:
    my_sec.pt3dadd(10, 0, 0, 2)
for mul, my_sec in [[1, sec[1]], [-1, sec[2]]]:
    my_sec.pt3dadd(10 * (1 + h.cos(mul * h.PI / 6)), 10 * h.sin(mul * h.PI / 6), 0, 2)


# initial conditions
def ca_init(node):
    if node in sec[0]:
        return 1 * µM
    else:
        return 100 * nM


# setup the pure diffusion test
cyt = rxd.Region(sec[0].wholetree(), name="cyt", nrn_region="i")
ca = rxd.Species(
    cyt, name="ca", initial=ca_init, d=1 * µm ** 2 / ms, charge=2, atolscale=1e-6
)

# indicate 3D simulation
rxd.set_solve_type(domain=sec[0].wholetree(), dimension=3)

# measure the total calcium in mM * µm ** 3 (is there a better unit?)
def total_ca():
    return sum(node.concentration * node.volume for node in ca.nodes)


# analysis routine, compare totals at different time points
def do_analysis(method, tstops=np.logspace(0, 5) * ms):

    if method == "fixed":
        cvode_active = False
    elif method == "variable":
        cvode_active = True
    else:
        raise Exception("unsupported integration method")

    h.CVode().active(cvode_active)
    print(f"*** {method} step integration ***")

    # initial membrane potential doesn't matter for this simulation
    # but initialization is required
    h.finitialize(-65 * mV)

    initial_total = total_ca()
    print(f"After initialization, total calcium = {initial_total} mM * µm ** 3")

    # advance until tstop (event ensures variable step does not go past tstop)
    percent_changes = []
    for tstop in tstops:
        h.CVode().event(tstop)
        h.continuerun(tstop)

        ending_total = total_ca()
        percent_change = 100 * (initial_total - ending_total) / initial_total
        percent_changes.append(abs(percent_change))
        print(f"At t = {h.t}, total calcium = {ending_total} mM * µm ** 3")
        print(f"    Change: {percent_change}%")
        print(
            f"    Maximum variation: {max(ca.nodes.concentration) - min(ca.nodes.concentration)} mM"
        )
        print()

    return p9.geom_line(
        data=pd.DataFrame(
            {"t": tstops, "percent error": percent_changes, "method": method}
        )
    )


# do the analysis and generate the graph
fixed_step_analysis = do_analysis("fixed")
variable_step_analysis = do_analysis("variable")

print(f"Total number of voxels: {len(ca.nodes)}")

print(
    p9.ggplot(p9.aes(x="t", y="percent error", color="method"))
    + fixed_step_analysis
    + variable_step_analysis
    + p9.scale_x_continuous(trans="log10")
    + p9.scale_y_continuous(trans="log10")
)
