from neuron import h, rxd
from neuron.units import ms, mV, µm, mM
import plotnine as p9
import pandas as pd
import multiprocessing as mp
import math
from functools import lru_cache

h.load_file("stdrun.hoc")

# centered at 76.5
start_x = 70 * µm
stop_x = 83 * µm
initial_concentration = 1 * mM
D = 1
tstop = 50 * ms


def fundamental_solution(t, D, distance):
    return 1 / math.sqrt(4 * math.pi * D * t) * math.exp(-(distance ** 2) / (4 * D * t))


@lru_cache(maxsize=None)
def solution_from_interval(
    x, t, D, start_x, stop_x, initial_concentration, interval_steps=100
):
    """solution of diffusion on an infinite line from a finite interval source"""
    assert start_x < stop_x
    dx = (stop_x - start_x) / interval_steps
    return initial_concentration * sum(
        dx * fundamental_solution(t, D, x - (start_x + (i + 0.5) * dx))
        for i in range(interval_steps)
    )


def solution_from_interval_with_boundaries(
    x, t, D, start_x, stop_x, initial_concentration, interval_steps=100
):
    # 153 µm is the total length of the line
    # in principle, the mathematical solution involves an infinite sum
    return (
        solution_from_interval(
            x, t, D, start_x, stop_x, initial_concentration, interval_steps
        )
        + solution_from_interval(x, t, D, -stop_x, -start_x, 0, interval_steps)
        + solution_from_interval(
            x, t, D, 153 - stop_x, 153 - start_x, 0, interval_steps
        )
    )


def main(dx=0.25, ics_partial_volume_resolution=2):
    # NOTE: NEURON's default is to use an ics_partial_volume_resolution of 2.
    print("starting setup...")
    # set accuracy of volumes
    rxd.options.ics_partial_volume_resolution = ics_partial_volume_resolution
    # small time step (half of NEURON's default)... needed for numerical stability with small dx
    h.dt = 0.0125 * ms
    # construct 4 axons
    # axon1: just 1D
    # axon2: 3D in the middle, 1D on edges
    # axon3: 3D everywhere
    # axon4: 3D on edges, 1D in the middle
    axon1 = h.Section(name="axon1")
    axon2a = h.Section(name="axon2a")
    axon2b = h.Section(name="axon2b")
    axon2c = h.Section(name="axon2c")
    axon3 = h.Section(name="axon3")
    axon4a = h.Section(name="axon4a")
    axon4b = h.Section(name="axon4b")
    axon4c = h.Section(name="axon4c")

    axon2b.connect(axon2a)
    axon2c.connect(axon2b)
    axon4b.connect(axon4a)
    axon4c.connect(axon4b)

    axon3.L = axon1.L = 153 * µm
    axon2a.L = axon2b.L = axon2c.L = axon4a.L = axon4b.L = axon4c.L = 51 * µm
    for sec in [axon1, axon2a, axon2b, axon2c, axon3, axon4a, axon4b, axon4c]:
        sec.diam = 2 * µm
        sec.nseg = int(sec.L) * 2

    rxd.set_solve_type([axon2b, axon3, axon4a, axon4c], dimension=3)

    # initialization rule
    def init_concentration(node):
        return initial_concentration if start_x < node.x3d < stop_x else 0

    # set up the model
    cytosol = rxd.Region(
        [axon1, axon2a, axon2b, axon2c, axon3, axon4a, axon4b, axon4c],
        name="cytosol",
        nrn_region="i",
        dx=dx,
    )
    ca = rxd.Species(cytosol, name="ca", d=D, charge=2, initial=init_concentration)

    # run the simulation
    print("initializing...")
    h.finitialize(-65 * mV)
    initial_mass = sum(
        node.concentration * node.volume
        for node in ca.nodes
        if node.sec in [axon2a, axon2b, axon2c]
    )
    print(
        "volume in 1D part:",
        sum(node.volume for node in ca.nodes if node.sec in [axon2a, axon2c]),
    )
    print(
        "volume in 3D part:",
        sum(node.volume for node in ca.nodes if node.sec in [axon2b]),
    )
    print("initial mass:", initial_mass)
    print("running...")
    h.continuerun(tstop)
    ending_mass = sum(
        node.concentration * node.volume
        for node in ca.nodes
        if node.sec in [axon2a, axon2b, axon2c]
    )
    print("ending mass:", ending_mass)
    print(
        f"true concentration at x=76.5: {solution_from_interval_with_boundaries(76.5 * µm, tstop, D, start_x, stop_x, initial_concentration)}"
    )
    print(
        f"true concentration 1/3 of the way in: {solution_from_interval_with_boundaries(51 * µm, tstop, D, start_x, stop_x, initial_concentration)}"
    )
    print(
        f"change over 1 timestep in true concentration 1/3 of the way in: {solution_from_interval_with_boundaries(51 * µm, tstop + h.dt, D, start_x, stop_x, initial_concentration) - solution_from_interval_with_boundaries(51 * µm, 50 * ms, D, start_x, stop_x, initial_concentration)}"
    )


    # prepare the data
    # note: axon1, axon2, axon3, and axon4 all have same number of segments
    sec2id = {
        axon1: "1D",
        axon2a: "3D middle",
        axon2b: "3D middle",
        axon2c: "3D middle",
        axon3: "full 3D",
        axon4a: "3D edges",
        axon4b: "3D edges",
        axon4c: "3D edges",
    }
    all_nodes = [
        node
        for sec in [axon1, axon2a, axon2b, axon2c, axon3, axon4a, axon4b, axon4c]
        for node in ca.nodes(sec)
    ]
    all_ids = [sec2id[node.sec] for node in all_nodes]
    data = pd.DataFrame(
        {
            "x": [node.x3d for node in all_nodes],
            "vol": [node.volume for node in all_nodes],
            "cai": [node.concentration for node in all_nodes],
            "id": all_ids,
            "true_cai": [
                solution_from_interval_with_boundaries(
                    node.x3d, h.t, D, start_x, stop_x, initial_concentration
                )
                for node in all_nodes
            ],
        }
    )
    data["id"] = data["id"].astype("category")
    data["dx"] = dx
    data["ics_partial_volume_resolution"] = ics_partial_volume_resolution
    data["mass"] = data["vol"] * data["cai"]
    data["average_cai"] = data.groupby(["x", "id"]).mass.transform("sum") / data.groupby(["x", "id"]).vol.transform("sum")

    data["abs_error"] = data["average_cai"] - data["true_cai"]
    data["unaveraged_abs_error"] = data["cai"] - data["true_cai"]
    data["rel_error"] = data["abs_error"] / data["true_cai"]

    # Plot
    p9.options.figure_size = (4, 1)
    g = (
        p9.ggplot(data, p9.aes(x="x", y="cai", color="id"))
        + p9.geom_line()
        + p9.labs(x="Position (µm)", y="Ca$^{2+}$ concentration (mM)")
    )
    g.save(f"comparison-to-truth-dx-{dx}-ics-{ics_partial_volume_resolution}.pdf")

    # Plot
    p9.options.figure_size = (4, 3)
    g = (
        p9.ggplot(data, p9.aes(x="x", y="abs_error", color="id"))
        + p9.geom_line()
        + p9.labs(x="Position (µm)", y="Absolute error (mM)", title=f"dx={dx}")
    )
    g.save(f"comparison-to-truth-abs-error-dx-{dx}-ics-{ics_partial_volume_resolution}.pdf")

    # print maximum absolute error

    g = (
        p9.ggplot(data, p9.aes(x="x", y="rel_error", color="id"))
        + p9.geom_line()
        + p9.labs(x="Position (µm)", y="relative error", title=f"dx={dx}")
    )
    g.save(f"comparison-to-truth-rel-error-dx-{dx}-ics-{ics_partial_volume_resolution}.pdf")

    return data


if __name__ == "__main__":
    with mp.Pool() as pool:
        data = pd.concat(pool.starmap(main, [(0.25, 2),  (0.25, 6), (0.125, 2)]))

    data["original_dx"] = data["dx"]
    data["dx"] = data["dx"].astype("category")

    p9.options.figure_size = (4, 3)
    g = (
        p9.ggplot(data[data["ics_partial_volume_resolution"] == 2], p9.aes(x="x", y="abs_error", color="id", linetype="dx"))
        + p9.geom_line()
        + p9.labs(x="Position (µm)", y="Signed absolute error (mM)")
    )
    g.save(f"comparison-to-truth-abs-error-by-dx.pdf")

    data["ics_partial_volume_resolution"] = data["ics_partial_volume_resolution"].astype("category")
    g = (
        p9.ggplot(data[data["original_dx"] == 0.25], p9.aes(x="x", y="abs_error", color="id", linetype="ics_partial_volume_resolution"))
        + p9.geom_line()
        + p9.labs(x="Position (µm)", y="Signed absolute error (mM)")
    )
    g.save(f"comparison-to-truth-abs-error-by-ics.pdf")

    pd.set_option('precision', 12)

    print("Averaged concentration error:")
    print(data.groupby(["dx", "ics_partial_volume_resolution", "id"]).abs_error.agg(["max", "min"]))

    print("Unaveraged:")
    print(data.groupby(["dx", "ics_partial_volume_resolution", "id"]).unaveraged_abs_error.agg(["max", "min"]))
