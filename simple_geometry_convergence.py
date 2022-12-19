import time
from neuron import h, rxd
import multiprocessing
import numpy as np
import sqlite3
import tqdm
import pandas as pd

DB_FILENAME = "simple_geometry_convergence.db"

try:
    with sqlite3.connect(DB_FILENAME) as conn:
        old_data = pd.read_sql("SELECT * FROM data", conn)
except:
    old_data = pd.DataFrame({"dx": [], "resolution": []})


def run_sim(dx, res=2, L=20, diam=2):
    start = time.perf_counter()

    rxd.options.ics_partial_volume_resolution = res

    dend = h.Section(name="dend")
    dend.L = L
    dend.diam = diam
    rxd.set_solve_type(dimension=3)

    cyt = rxd.Region([dend], name="cyt", dx=dx)
    ca = rxd.Species(cyt, name="ca", charge=2)

    true_volume = h.PI * dend.diam ** 2 * 0.25 * dend.L
    true_area = h.PI * dend.diam * dend.L + 0.5 * h.PI * dend.diam ** 2

    data = pd.DataFrame(
        {
            "dx": [dx],
            "L": L,
            "diam": diam,
            "surface_area": sum(ca.nodes.surface_area),
            "volume": sum(ca.nodes.volume),
            "surface_area_relative_error": 1 - sum(ca.nodes.surface_area) / true_area,
            "volume_relative_error": 1 - sum(ca.nodes.volume) / true_volume,
            "runtime": time.perf_counter() - start,
            "resolution": res
        }
    )
    with sqlite3.connect(DB_FILENAME) as conn:
        data.to_sql("data", conn, if_exists="append", index=False)



if __name__ == "__main__":
    for dx in tqdm.tqdm(np.logspace(np.log10(0.5), -2)):
        for res in [10, 8, 6, 4, 2]:
            if any((old_data["dx"] == dx) & (old_data["resolution"] == res)):
                continue
            else:
                p = multiprocessing.Process(target=run_sim, args=(dx, res))
                p.start()
                p.join()
    