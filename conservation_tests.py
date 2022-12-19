from neuron import h, rxd
from matplotlib import pyplot
from neuron.units import s, µm, nM, mV
import pandas as pd
import numpy as np
import sqlite3
import sys
h.load_file('stdrun.hoc')


def line(dx=0.25, dt=0.025, source='1d', hybrid=False):
    dend1 = h.Section(name='dend1')
    dend1.diam = 2
    dend1.nseg = 11
    dend1.L = 10
   
    dend2 = h.Section(name='dend2')
    dend2.diam = 2
    dend2.nseg = 11
    dend2.L = 10
    dend2.connect(dend1)

    diff_constant = 1

    r = rxd.Region(h.allsec(),dx=dx)
    if hybrid:
        rxd.set_solve_type([dend2], dimension=3)
    else:
        rxd.set_solve_type([dend1, dend2], dimension=3)

    source_sec = dend1 if source == '1d' else dend2

    ca = rxd.Species(r, d=diff_constant, atolscale=nM,
                        initial=lambda nd: 1 * µm if nd in source_sec and
                                                nd.x < 0.75 else 0)
    h.finitialize(-70 * mV)
    initial_amount = (np.array(ca.nodes.concentration) * np.array(ca.nodes.volume)).sum()
    h.continuerun(100 * s)
    final_amount = (np.array(ca.nodes.concentration) * np.array(ca.nodes.volume)).sum()

    return initial_amount, final_amount


def split(dx=0.25, dt=0.025, align=False, source='1d', hybrid=False):

    dend1 = h.Section(name='dend1')
    dend1.pt3dclear()
    dend1.pt3dadd(0,0,0,2)
    dend1.pt3dadd(10,0,0,2)
    dend1.nseg = 11

    dend2 = h.Section(name='dend2')
    dend2.pt3dclear()
    dend2.pt3dadd(10 - dx, 0, 0, 2)
    dend2.pt3dadd(10, 0, 0, 2)
    if align:
        dend2.pt3dadd(20, 0, 0, 2)
    else:
        dend2.pt3dadd(10 + 5 *(3**0.5), 5, 0, 2)
    dend2.nseg = 11
    dend2.connect(dend1)

    dend3 = h.Section(name='dend3')
    dend3.pt3dclear()
    dend3.pt3dadd(10 - dx, 0, 0, 2)
    dend3.pt3dadd(0, 0, 0, 2)
    if align:
        dend3.pt3dadd(20, 0, 0, 1)
    else:
        dend2.pt3dadd(10 + 5 *(3**0.5), -5, 0, 2)
    dend3.nseg = 11
    dend3.connect(dend1)
    diff_constant = 1

    r = rxd.Region([dend1, dend2, dend3], dx=dx)
    if hybrid:
        rxd.set_solve_type([dend1], dimension=3)
    else:
        rxd.set_solve_type([dend1, dend2, dend3], dimension=3)
    source_sec = [dend1] if source == '3d' else [dend2, dend3]
    ca = rxd.Species(r, d=diff_constant, atolscale=nM,
                        initial=lambda nd: 1 * µm if nd.sec in source_sec and 
                                                nd.x < 0.75 else 0)
    h.finitialize(-70 * mV)
    initial_amount = (np.array(ca.nodes.concentration) * np.array(ca.nodes.volume)).sum()
    h.continuerun(100 * s)

    final_amount = (np.array(ca.nodes.concentration) * np.array(ca.nodes.volume)).sum()

    return initial_amount, final_amount


if __name__ == "__main__":
    dx = float(sys.argv[1])
    if 'cvode' in sys.argv[2]:
        h.CVode().active(True)
        h.CVode().atol
        dt = 'cvode'
    else:
        dt = float(sys.argv[2])
        h.dt = dt
    source = sys.argv[3]
    model = sys.argv[4]
    hybrid = sys.argv[5] == 'hybrid'
    print(f"processing {model} in {hybrid} with dx={dx} dt={dt}")
    conn = sqlite3.connect("conservation_tests.db")
    if model == "split_align":
        initial_amount, final_amount = split(dx=dx, dt=dt, align=True,
                                                    source=source,
                                                    hybrid=hybrid)
    elif model == "split_y":
        initial_amount, final_amount = split(dx=dx, dt=dt, align=False,
                                                    source=source,
                                                    hybrid=hybrid)

    else:
        initial_amount, final_amount = line(dx=dx, dt=dt, source=source,
                                                   hybrid=hybrid)
   
    data = pd.DataFrame({
        "hybrid": hybrid,
        "dx": dx,
        "dt": dt,
        "source": source,
        "initial": initial_amount,
        "final_amount": final_amount,
        "diff": initial_amount-final_amount,
        "ratio": 1-final_amount/initial_amount
    },index=[f"{source}_{dx}_{dt}_{hybrid}"])
    data.to_sql(model, conn, if_exists="append", index=False)


