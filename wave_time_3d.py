import multiprocessing
import random
import sqlite3
import pandas as pd
from neuron import h, rxd
from neuron.units import mV, ms

h.load_file("stdrun.hoc")

THRESHOLD_CONCENTRATION = 0.5
NUM_ORIENTATIONS = 100

try:
    with sqlite3.connect("wave_time_3d.db") as conn:
        old_data = pd.read_sql("SELECT * FROM data", conn)
except:
    old_data = pd.DataFrame({"theta": [], "phi": [], "alpha": [], "dx": []})

def on_stopevent():
    h.stoprun = True


def save_data(theta, phi, dx, alpha, length, diam, speed, error, sim_time):
    # connect to the database (or create it if it doesn't exist)
    conn = sqlite3.connect("wave_time_3d.db")
    c = conn.cursor()

    # ensure the table exists
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS data (
            theta REAL,
            phi REAL,
            dx REAL,
            alpha REAL,
            length REAL,
            diam REAL,
            speed REAL,
            relative_error REAL,
            sim_time REAL
        )
        """
    )

    # store the data
    c.execute(
        "INSERT INTO data VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (theta, phi, dx, alpha, length, diam, speed, error, sim_time),
    )
    conn.commit()
    conn.close()


def run_sim(theta, phi, dx, alpha=0.25, L=251, diam=2):
    # theta, phi are polar angle and azimuthal angle, respectively
    # per ISO 80000-2:2019... this is physics style not math convention
    # theta \in [0, \pi), phi \in [0, 2*pi)
    import time
    import numpy as np

    start = time.perf_counter()

    # setup the model geometry
    dend = h.Section(name="dend")
    dend.pt3dadd(0, 0, 0, diam)
    dend.pt3dadd(
        L * h.cos(phi) * h.sin(theta),
        L * h.sin(phi) * h.sin(theta),
        L * h.cos(theta),
        diam,
    )
    dend.nseg = 251

    # setup the model kinetics
    cyt = rxd.Region([dend], name="cyt", dx=dx, nrn_region="i")
    c = rxd.Species(
        cyt, d=1, name="c", initial=lambda node: 1 if node.x * dend.L < 50 else 0
    )
    wave_reaction = rxd.Rate(c, -c * (1 - c) * (alpha - c))

    # integration options
    rxd.nthread(4)
    rxd.set_solve_type(dimension=3)

    # the locations we'll monitor
    pt1 = dend(100 / dend.L)
    pt2 = dend(200 / dend.L)
    distance = h.distance(pt1, pt2)

    # monitor concentration timecourses at the points above
    threshold = h.ref(THRESHOLD_CONCENTRATION)
    c_pt1 = h.Vector().record(pt1._ref_ci)
    c_pt2 = h.Vector().record(pt2._ref_ci)
    t = h.Vector().record(h._ref_t)

    # stop simulation when pt2 crosses the threshold
    ste = h.StateTransitionEvent(1)
    ste.transition(0, 0, pt2._ref_ci, threshold, on_stopevent)

    def on_finitialize():
        ste.state(0)

    fih = h.FInitializeHandler(on_finitialize)

    # use variable step integration
    h.CVode().active(True)
    h.CVode().atol(1e-6)

    # actually run the simulation
    h.finitialize(-65 * mV)
    h.continuerun(3000 * ms)
    print(f"end time: {h.t}")

    # interpolate to estimate the crossing times
    pt1_crossing_time = np.interp(THRESHOLD_CONCENTRATION, c_pt1, t)
    pt2_crossing_time = np.interp(THRESHOLD_CONCENTRATION, c_pt2, t)

    measured_speed = distance / (pt2_crossing_time - pt1_crossing_time)
    expected_speed = 2 ** 0.5 * (0.5 - alpha)
    speed_error = abs(1 - measured_speed / expected_speed)

    print(f"pt1_crossing_time = {pt1_crossing_time}")
    print(f"pt2_crossing_time = {pt2_crossing_time}")
    print(f"speed = {measured_speed}")
    print(f"expected speed = {expected_speed}")
    print(f"relative error = {speed_error}")

    finished = time.perf_counter()
    print(f"elapsed time = {finished - start} s")
    save_data(
        theta, phi, dx, alpha, L, diam, measured_speed, speed_error, finished - start
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
    for dx in [2 ** -1, 2 ** -2, 2 ** -3, 2 ** -4, 1, 2 ** -5]:
        for alpha in [0.25, 0.15, 0.35]:
            for theta, phi in orientations:
                if any((old_data["dx"] == dx) & (old_data["alpha"] == alpha) & (old_data["theta"] == theta) & (old_data["phi"] == phi)):
                    print(f"Skipping: dx={dx}, alpha={alpha}, theta={theta}, phi={phi}")
                else:
                    print(f"Running: dx={dx}, alpha={alpha}, theta={theta}, phi={phi}")
                    p = multiprocessing.Process(
                        target=run_sim, args=(theta, phi, dx, alpha)
                    )
                    p.start()
                    p.join()
