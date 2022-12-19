from neuron import h, rxd
from neuron.units import mV, ms
import matplotlib.pyplot as plt

h.load_file("stdrun.hoc")
import time
import multiprocessing


def run_sim(tstop, dim, d):
    start = time.perf_counter()

    soma = h.Section(name="soma")
    soma.L = 10
    soma.diam = 10

    soma.insert(h.hh)

    rxd.set_solve_type(dimension=dim)
    cyt = rxd.Region([soma], nrn_region="i", name="cyt")
    na = rxd.Species(cyt, name="na", charge=1, initial=h.nai0_na_ion, d=d)

    ic = h.IClamp(soma(0.5))
    ic.amp = 0.1
    ic.delay = 0 * ms
    ic.dur = tstop

    t = h.Vector().record(h._ref_t)
    v = h.Vector().record(soma(0.5)._ref_v)
    sodium = h.Vector().record(soma(0.5)._ref_nai)

    h.finitialize(-65 * mV)
    h.continuerun(tstop)

    finished = time.perf_counter()
    return {
        "time": finished - start,
        "num_nodes": len(na.nodes),
        "total_na": sum(node.volume * node.concentration for node in na.nodes),
        "surface_area": sum(node.surface_area for node in na.nodes),
        "volume": sum(node.volume for node in na.nodes),
        "na_d": na.d,
        "t": t,
        "v": v,
        "sodium": sodium,
        "dimension": dim,
    }



def main():
    tstop = 100 * ms

    fig = plt.figure(figsize=(6, 6))
    voltage_axis = fig.add_subplot(2, 1, 1)
    voltage_axis_zoom = voltage_axis.inset_axes([0.01, 0.38, 0.3, 0.6])
    sodium_axis = fig.add_subplot(2, 1, 2)
    sodium_axis_zoom = sodium_axis.inset_axes([0.01, 0.38, 0.3, 0.6])

    pool = multiprocessing.Pool()

    for data in pool.starmap(
        run_sim,
        [
            [tstop, 3, 1e-4],
            [tstop, 3, 1e-3],
            [tstop, 3, 1e-2],
            [tstop, 3, 1e-1],
            [tstop, 3, 1],
            [tstop, 1, 0],
        ],
    ):
        if data["dimension"] == 3:
            for axis in [voltage_axis, voltage_axis_zoom]:
                axis.plot(data["t"], data["v"], label=f"D={data['na_d']}")
            for axis in [sodium_axis, sodium_axis_zoom]:
                axis.plot(data["t"], data["sodium"], label=f"D={data['na_d']}")
        else:
            for axis in [voltage_axis, voltage_axis_zoom]:
                axis.plot(data["t"], data["v"], "k--", label="1D")
            for axis in [sodium_axis, sodium_axis_zoom]:
                axis.plot(data["t"], data["sodium"], "k--", label="1D")


        print(
            f"""
dimension = {data["dimension"]}
na.d = {data["na_d"]}
len(na.nodes) = {data["num_nodes"]}
total(na) = {data["total_na"]}
surface area = {data["surface_area"]}
volume = {data["volume"]}
elapsed time = {data["time"]}
"""
        )

        if data["dimension"] == 1:
            print("(surface area interpreted differently for 1D (no edge faces), but total current flux is the same)")

    sodium_axis.legend(ncol=1, loc="upper right", facecolor="white", framealpha=1)

    voltage_axis.set_xlim(0, tstop)
    voltage_axis.set_xticklabels([])
    sodium_axis.set_xlim(0, tstop)
    sodium_axis.set_xlabel("Time (ms)")
    sodium_axis.set_ylabel("[Na$^+$] (mM)")
    sodium_axis.set_ylim(10, 15)
    voltage_axis.set_ylim(-80, 40)
    voltage_axis.set_ylabel("Membrane potential (mV)")
    voltage_axis_zoom.set_xlim(79, 84)
    voltage_axis_zoom.set_ylim(-75, 30)
    sodium_axis_zoom.set_xlim(68, 76)
    sodium_axis_zoom.set_ylim(10.2, 11.4)
    voltage_axis_zoom.set_xticks([])
    voltage_axis_zoom.set_yticks([])
    sodium_axis_zoom.set_xticks([])
    sodium_axis_zoom.set_yticks([])
    voltage_axis.indicate_inset_zoom(voltage_axis_zoom, edgecolor="black")
    sodium_axis.indicate_inset_zoom(sodium_axis_zoom, edgecolor="black")
    # subplot labeling based on https://stackoverflow.com/questions/18344939/matplotlib-panel-label-out-of-the-box-above-the-ylabel
    voltage_axis.text(-0.1, 1.15, "A", fontweight="bold", va="top", ha="right", transform=voltage_axis.transAxes, fontsize=16)
    sodium_axis.text(-0.1, 1.15, "B", fontweight="bold", va="top", ha="right", transform=sodium_axis.transAxes, fontsize=16)

    fig.savefig("response_to_currents.pdf")

    plt.show()


if __name__ == "__main__":
    main()
