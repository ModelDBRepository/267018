from neuron import h, rxd
import sys
from volume_functions_truebound import percent, find_volume, find_bounding_box_volume, \
    find_skeleton_maxPoints, find_skeleton_minPoints, qualify_sections, find_3d_extrema

h.load_file('stdrun.hoc')
h.load_file('stdlib.hoc')

h.load_file('import3d.hoc')

class Cell:
    def __init__(self, morphology):
        x = h.Import3d_SWC_read()
        x.input(morphology)
        w = h.Import3d_GUI(x, False)
        w.instantiate(self)


def main():
    """ Takes one command line argument:
    the name of the file containing the morphology definition
    Directs output to filename_outputs.txt; file is overwritten if exists"""
    ogstdout = sys.stdout
    if len(sys.argv) > 1:
        args = sys.argv[1:]
        if args is not None:
            filename = "".join(["truebound_outputs/", args[0][4:], "_outputs_tb.txt"])
            with open(filename, 'w') as f:
                sys.stdout = f
                print(f"\n")

            with open(filename, 'a') as f:
                sys.stdout = f
                morphology = args[0]
                print(f"MORPHOLOGY: {str(morphology)}")
                print(f"MORPHOLOGY: {str(morphology)}", file=ogstdout)

                # vol=0
                x = Cell(morphology)
                skeleton_mins = find_skeleton_minPoints()
                skeleton_maxs = find_skeleton_maxPoints()
                vol = find_volume(h.allsec(), _3d=False)

                qualifying_sections = qualify_sections(x.all, [skeleton_mins, skeleton_maxs])

                print(f"Out of {len(x.all)}, {len(qualifying_sections)} sections qualify for voxelization", file=ogstdout)
                print(f"Those sections are: {qualifying_sections}", file=ogstdout)

                rxd.set_solve_type(qualifying_sections, dimension=3)
                reg = rxd.Region(qualifying_sections, dx=0.10)
                s = rxd.Species(reg)

                extrema = find_3d_extrema(nodelist=s)
                mins = extrema[0]
                maxs = extrema[1]
                box_volume = find_bounding_box_volume(maxs=maxs, mins=mins, _3d=True)   #sed _3d=False if skeleton
                print(f"approx extrema (just qualifying sections included): min: {mins}, max: {maxs}")
                fin_res = percent(vol, box_volume)  # vol here is skeleton volume

                # rxd.set_solve_type(x.all, dimension=3)
                # reg2 = rxd.Region(x.all, dx=0.10)
                # s2 = rxd.Species(reg2)
                # extrema2 = find_3d_extrema(nodelist=s)
                # truemins = extrema2[0]
                # truemaxs = extrema2[1]

                # print(f"true extrema: min: {truemins}, max: {truemaxs}")

                print("FINAL RESULTS:\n")
                print(f"Total volume of {str(morphology)} cell: {vol}")
                print(f"Total {str(morphology)} bounding box volume: {box_volume}\n\n")
                print(f"Morphology {str(morphology)} cell in bounding box: {fin_res * 100:.5f} %\n")

                sys.stdout = ogstdout
    else:
        raise ValueError


if __name__ == "__main__":
    main()
