# import numpy as np
from neuron import h

# pi = np.pi


def find_skeleton_maxPoints():
    # if not _3d and nodelist is None:
    xmax = max(max(sec.x3d(i) for i in range(sec.n3d())) for sec in h.allsec())
    ymax = max(max(sec.y3d(i) for i in range(sec.n3d())) for sec in h.allsec())
    zmax = max(max(sec.z3d(i) for i in range(sec.n3d())) for sec in h.allsec())


    print(f"MAX NODES: x={xmax}, y={ymax}, z={zmax}")

    return [xmax, ymax, zmax]


def find_skeleton_minPoints():
    # if nodelist is None and not _3d:
    xmin = min(min(sec.x3d(i) for i in range(sec.n3d())) for sec in h.allsec())
    ymin = min(min(sec.y3d(i) for i in range(sec.n3d())) for sec in h.allsec())
    zmin = min(min(sec.z3d(i) for i in range(sec.n3d())) for sec in h.allsec())
    print(f"MIN NODES: x={xmin}, y={ymin}, z={zmin}")

    return [xmin, ymin, zmin]

def find_3d_extrema(nodelist):
    xmin = min(node.x3d for node in nodelist.nodes)
    ymin = min(node.y3d for node in nodelist.nodes)
    zmin = min(node.z3d for node in nodelist.nodes)
    xmax = max(node.x3d for node in nodelist.nodes)
    ymax = max(node.y3d for node in nodelist.nodes)
    zmax = max(node.z3d for node in nodelist.nodes)

    return [[xmin, ymin, zmin], [xmax, ymax, zmax]]


def find_sec_r(section):
    sec_r = max((seg.diam/2) for seg in section)
    return sec_r


def section_extrema(sec):
    xmin = min(sec.x3d(i) for i in range(sec.n3d()))
    ymin = min(sec.y3d(i) for i in range(sec.n3d()))
    zmin = min(sec.z3d(i) for i in range(sec.n3d()))
    xmax = max(sec.x3d(i) for i in range(sec.n3d()))
    ymax = max(sec.y3d(i) for i in range(sec.n3d()))
    zmax = max(sec.z3d(i) for i in range(sec.n3d()))

    return [[xmin, ymin, zmin], [xmax, ymax, zmax]]


def qualify_sections(all_sections, extrema):
    # extrema = [[xmin, ymin, zmin], [xmax, ymax, zmax]] <-- skeleton extrema
    qualifying_sections = set()
    for sec in all_sections:
        rmax = find_sec_r(sec)
        i = 0
        sec_extrema = section_extrema(sec)
        for extreme in extrema:
            for axis in range(3):
                if i == 0:  # min
                    if (extrema[0][axis] - rmax) <= sec_extrema[0][axis] <= extrema[0][axis]:
                        qualifying_sections.add(sec)
                        break
                else:   # max
                    if extrema[1][axis] <= sec_extrema[1][axis] <= (extrema[1][axis] + rmax):
                        qualifying_sections.add(sec)
                        break
            i+=1

    return qualifying_sections


def find_volume(list_of_sections, _3d=False):
    volume = 0

    if not _3d:
        for sec in list_of_sections:
            for seg in sec:
                volume += seg.volume()
    else:
        sum(list_of_sections.volume)

    print(f"neuron volume: {volume}, 3d? = {_3d}")
    # print(f"Returned list of xs, ys, zs, each {len(xs)} long")

    return volume


def find_bounding_box_volume(maxs, mins, _3d=False):
    vol = 1
    for i in range(3):
        vol *= (maxs[i] - mins[i])
        print(f"dimension: [0=x, 1=y, 2=z] = {i}, length: {maxs[i] - mins[i]}")
    print(f"bounding box volume: {vol}, _3d?={_3d}")
    return vol


def percent(a, b):
    print(f"percentage: {a/b}")
    return a / b
