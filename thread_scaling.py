import time
import multiprocessing
import sqlite3
import pandas as pd
from neuron import h, rxd
from neuron.units import mV, ms, um, mM

h.load_file("stdrun.hoc")
h.load_file("import3d.hoc")

DB_FILENAME = "thread_scaling.db"
SWC_FILENAME = "B4-CA1-L-D63x1zACR3_1.CNG.swc.txt"
NUM_RUNS = 3

class Cell:
    def __init__(self, dx):
        cell = h.Import3d_SWC_read()
        cell.input(SWC_FILENAME)
        i3d = h.Import3d_GUI(cell, False)
        i3d.instantiate(self)
        self.start = self.soma[0]
        self.name = "cell"
        self.dx = dx

class Cylinder:
    def __init__(self, dx):
        self.all = [h.Section(name=f"dend{i}") for i in range(2)]
        for dend in self.all:
            dend.L = 25 * um
            dend.diam = 1 * um
        self.all[1].connect(self.all[0])
        self.start = self.all[1]
        self.name = "cylinder"
        self.dx = dx

def diffusion_only(obj):
    return {
        "name": "diffusion",
        "regions": [cyt := rxd.Region(obj.all, name="cyt", dx=obj.dx)],
        "species": [rxd.Species(cyt, d=1*um**2/ms, initial=lambda node: 1 if node in obj.start else 0)]
    }

def bistable(obj):
    return {
        "name": "bistable-wave",
        "regions": [cyt := rxd.Region(obj.all, name="cyt", dx=obj.dx)],
        "species": [c := rxd.Species(cyt, d=1*um**2/ms, initial=lambda node: 1 * mM if node in obj.start else 0)],
        "reactions": [rxd.Rate(c, -c * (1 * mM - c) * (0.3 * mM - c))]
    }

def cawave(obj,
    caDiff = 0.08,
    ip3Diff = 1.41,
    cac_init = 1.e-4,
    ip3_init = 0.1,
    gip3r = 12040,
    gserca = 0.3913,
    gleak = 6.020,
    kserca = 0.1,
    kip3 = 0.15,
    kact = 0.4,
    ip3rtau = 2000,
    fc = 0.8,
    fe = 0.2,
    average_ca_inside = 0.0017
):
    cae_init = (average_ca_inside - cac_init * fc) / fe
    return {
        "name": "cawave",
        "regions": [
            cyt := rxd.Region(obj.all, nrn_region='i', geometry=rxd.FractionalVolume(fc, surface_fraction=1), dx=obj.dx),
            er := rxd.Region(obj.all, geometry=rxd.FractionalVolume(fe), dx=obj.dx),
            cyt_er_membrane := rxd.Region(obj.all, geometry=rxd.DistributedBoundary(1), dx=obj.dx)
        ],
        "species": [
            ca := rxd.Species([cyt, er], d=caDiff, name='ca', charge=2, initial=lambda node: cac_init if node in cyt else cae_init, atolscale=1e-6),
            ip3 := rxd.Species(cyt, d=ip3Diff, initial=lambda node: 2 * mM if node in obj.start else ip3_init),
            ip3r_gate_state := rxd.State(cyt_er_membrane, initial=0.8)
        ],
        "misc": [
            h_gate := ip3r_gate_state[cyt_er_membrane],
            minf := ip3[cyt] * 1000. * ca[cyt] / (ip3[cyt] + kip3) / (1000. * ca[cyt] + kact),
            k := gip3r * (minf * h_gate) ** 3
        ],
        "reactions": [
            serca := rxd.MultiCompartmentReaction(ca[cyt], ca[er], gserca / ((kserca / (1000. * ca[cyt])) ** 2 + 1), membrane=cyt_er_membrane, custom_dynamics=True),
            leak := rxd.MultiCompartmentReaction(ca[er], ca[cyt], gleak, gleak, membrane=cyt_er_membrane),
            ip3r := rxd.MultiCompartmentReaction(ca[er], ca[cyt], k, k, membrane=cyt_er_membrane),
            ip3rg := rxd.Rate(h_gate, (1. / (1 + 1000. * ca[cyt] / (0.3)) - h_gate) / ip3rtau)
        ]
    }


try:
    with sqlite3.connect(DB_FILENAME) as conn:
        old_data = pd.read_sql("SELECT * FROM data", conn)
except:
    old_data = pd.DataFrame({"nthread": [], "morphology": [], "kinetics": [], "dx": [], "runcount": [], "runtime": []})


def run_sim(nthread, morphology, kinetics, dx):
    # setup the model
    rxd.nthread(nthread)
    rxd.set_solve_type(dimension=3)
    morph = morphology(dx)
    my_kinetics = kinetics(morph)

    # skip if we've already done this
    if any((old_data["dx"] == dx) & (old_data["morphology"] == morph.name) & (old_data["kinetics"] == my_kinetics["name"]) & (old_data["nthread"] == nthread)):
        print(f"skipping: dx: {dx}, morph: {morph.name}, kinetics: {my_kinetics['name']}, nthread: {nthread}")
        return
    
    print(f"running: dx: {dx}, morph: {morph.name}, kinetics: {my_kinetics['name']}, nthread: {nthread}")

    # run the sim several times
    times = []
    for run in range(NUM_RUNS):
        print(f"  run #{run + 1}")
        initial_time = time.perf_counter()
        h.finitialize(-65 * mV)
        start_time = time.perf_counter()
        print(f"    initialization time: {start_time - initial_time}")
        h.continuerun(100 * ms)
        end_time = time.perf_counter()
        times.append(end_time - start_time)
        print(f"    elapsed: {end_time - start_time} s")
    
    # store the data in the database
    data = pd.DataFrame(
        {
            "nthread": nthread,
            "morphology": morph.name,
            "kinetics": my_kinetics['name'],
            "dx": dx,
            "runcount": range(NUM_RUNS),
            "runtime": times
        }
    )
    with sqlite3.connect(DB_FILENAME) as conn:
        data.to_sql("data", conn, if_exists="append", index=False)    


if __name__ == "__main__":
    for dx in [0.12, 0.06]:
        for morphology in [Cylinder, Cell]:
            for kinetics in [cawave, diffusion_only, bistable]:
                for nthread in [1, 2, 3, 4, 5, 6, 7, 8]:
                    p = multiprocessing.Process(
                        target=run_sim, args=(nthread, morphology, kinetics, dx)
                    )
                    p.start()
                    p.join()      

 