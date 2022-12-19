import time
from neuron import h, rxd

h.load_file("import3d.hoc")


class Cell:
    def __init__(self, filename):
        cell = h.Import3d_SWC_read()
        cell.input(filename)
        i3d = h.Import3d_GUI(cell, False)
        i3d.instantiate(self)


if __name__ == "__main__":
    import pandas as pd
    import sqlite3
    import sys

    conn = sqlite3.connect("discretization.db")
    filename = sys.argv[1]
    dx = float(sys.argv[2])
    print(f"processing {filename} at dx={dx}")
    cell = Cell(filename)
    rxd.set_solve_type(cell.all, dimension=3)
    cyt = rxd.Region(cell.all, name="cyt", dx=dx)
    x = rxd.Species(cyt, name="x")
    start = time.perf_counter()
    rxd.re_init()
    elapsed = time.perf_counter() - start
    print(f"elapsed time: {elapsed} sec")
    surface_area = sum(x.nodes.surface_area)
    volume = sum(x.nodes.volume)
    data = pd.DataFrame(
        {
            "morphology": [filename],
            "dx": [dx],
            "volume": [volume],
            "surface_area": [surface_area],
            "num_voxels": [len(x.nodes)],
            "num_surface_voxels": [
                len([node for node in x.nodes if node.surface_area])
            ],
            "discretization_time": [elapsed],
            "num_sections": [len(cell.all)],
            "sum_lengths": [sum([sec.L for sec in cell.all])],
        }
    )
    data.to_sql("morphology", conn, if_exists="append", index=False)
