import multiprocessing
import random
import sqlite3
import pandas as pd
from neuron import h, rxd
import sys

h.load_file("stdrun.hoc")

NUM_ORIENTATIONS = 1000

try:
    with sqlite3.connect("cylinder_convergence.db") as conn:
        old_data = pd.read_sql("SELECT * FROM data", conn)
except:
    old_data = pd.DataFrame({"theta": [], "phi": [], "volume": [], "area": [], "dx": []})


def save_data(theta, phi, dx, volume, area):
    # connect to the database (or create it if it doesn't exist)
    with sqlite3.connect("cylinder_convergence.db") as conn:
        c = conn.cursor()

        # ensure the table exists
        c.execute(
            """
            CREATE TABLE IF NOT EXISTS data (
                theta REAL,
                phi REAL,
                dx REAL,
                volume REAL,
                area REAL
            )
            """
        )

        # store the data
        c.execute(
            "INSERT INTO data VALUES (?, ?, ?, ?, ?)",
            (theta, phi, dx, volume, area),
        )
        conn.commit()


def run_sim(theta, phi, dx, L=5, diam=2):
    # theta, phi are polar angle and azimuthal angle, respectively
    # per ISO 80000-2:2019... this is physics style not math convention
    # theta \in [0, \pi), phi \in [0, 2*pi)
    import numpy as np


    # setup the model geometry
    dend = h.Section(name="dend")
    dend.pt3dadd(0, 0, 0, diam)
    dend.pt3dadd(
        L * h.cos(phi) * h.sin(theta),
        L * h.sin(phi) * h.sin(theta),
        L * h.cos(theta),
        diam,
    )

    rxd.options.ics_partial_volume_resolution = 1
    rxd.set_solve_type(dimension=3)

    cyt = rxd.Region([dend], name="cyt", dx=dx, nrn_region="i")
    c = rxd.Species(
        cyt, name="c"
    )

    volume = sum(node.volume for node in c.nodes)
    area = sum(node.surface_area for node in c.nodes)

    save_data(
        theta, phi, dx, volume, area
    )


if __name__ == "__main__":
    # ensure deterministic randomness
    random.seed(1)

    # pick random orientations
    orientations = [
        (random.random() * h.PI, random.random() * 2 * h.PI)
        for _ in range(NUM_ORIENTATIONS)
    ]

    # do the parameter study
    start_i = int(sys.argv[1])
    stop_i = int(sys.argv[2])
    for dx in [2**-1, 2**-1.5, 2**-2, 2**-2.5, 2**-3, 2**-3.5, 2**-4]:
        for theta, phi in orientations[start_i:stop_i]:
            if any((old_data["dx"] == dx) & (old_data["theta"] == theta) & (old_data["phi"] == phi)):
                print(f"Skipping: dx={dx}, theta={theta}, phi={phi}")
            else:
                print(f"Running: dx={dx}, theta={theta}, phi={phi}")
                p = multiprocessing.Process(
                    target=run_sim, args=(theta, phi, dx)
                )
                p.start()
                p.join()
