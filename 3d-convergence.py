from neuron import h, rxd
from neuron.units import ms, mV, µm, mM
import plotnine as p9
import pandas as pd
import multiprocessing as mp
import math
from functools import lru_cache
h.load_file('stdrun.hoc')

# centered at 100
start_x = 95 * µm
stop_x = 105 * µm
initial_concentration = 1 * mM
D = 1

def fundamental_solution(t, D, distance):
    return 1 / math.sqrt(4 * math.pi * D * t) * math.exp(-distance ** 2 / (4 * D * t))

@lru_cache(maxsize=None)
def solution_from_interval(x, t, D, start_x, stop_x, initial_concentration, interval_steps=100):
    """solution of diffusion on an infinite line from a finite interval source"""
    assert(start_x < stop_x)
    dx = (stop_x - start_x) / interval_steps
    return initial_concentration * sum(dx * fundamental_solution(t, D, x - (start_x + (i + 0.5) * dx)) for i in range(interval_steps))

def solution_from_interval_with_boundaries(x, t, D, start_x, stop_x, initial_concentration, interval_steps=100):
    # 200 µm is the total length of the line
    # in principle, the mathematical solution involves an infinite sum
    return (
        solution_from_interval(x, t, D, start_x, stop_x, initial_concentration, interval_steps)
        + solution_from_interval(x, t, D, -stop_x, -start_x, 0, interval_steps)
        + solution_from_interval(x, t, D, 200-stop_x, 200-start_x, 0, interval_steps)
    )

def main(dx=0.25):
    print("starting setup...")
    axon3 = h.Section(name="axon3")

    axon3.pt3dclear()
    axon3.pt3dadd(0, 0, 0, 1 * µm)
    axon3.pt3dadd(200 * µm, 0, 0, 1 * µm)
    
    rxd.set_solve_type(axon3, dimension=3)

    # initialization rule
    def init_concentration(node):
        return initial_concentration if start_x < node.x3d < stop_x else 0

    # set up the model
    cytosol = rxd.Region([axon3], name='cytosol', nrn_region="i", dx=dx)
    ca = rxd.Species(cytosol, name='ca', d=D, charge=2, initial=init_concentration)

    # run the simulation
    print("initializing...")
    h.finitialize(-65 * mV)
    initial_mass = sum(node.concentration * node.volume for node in ca.nodes)
    print("initial mass:", initial_mass)
    print("running...")
    h.continuerun(100 * ms)
    ending_mass = sum(node.concentration * node.volume for node in ca.nodes)
    print("ending mass:", ending_mass)


    # prepare the data
    all_nodes = ca.nodes
    data = pd.DataFrame({
        'x': [node.x3d for node in all_nodes],
        'cai': [node.concentration for node in all_nodes],
        "true_cai": [solution_from_interval_with_boundaries(node.x3d, h.t, D, start_x, stop_x, initial_concentration) for node in all_nodes]
    })
    data["abs_error"] = (data["cai"] - data["true_cai"]).abs()
    data["rel_error"] = data["abs_error"] / data["true_cai"]
    data["dx"] = f"{dx} µm"

    # Plot
    p9.options.figure_size = (4, 1)
    g = (
        p9.ggplot(data, p9.aes(x='x', y='cai'))
        + p9.geom_line()
        + p9.labs(x='Position (µm)', y='[Ca$^{2+}$] (mM)')
    )
    g.save(f'3d-convergence-dx-{dx}.pdf')
    g.save(f'3d-convergence-dx-{dx}.png')

    # Plot
    p9.options.figure_size = (4, 3)
    g = (
        p9.ggplot(data, p9.aes(x='x', y='abs_error'))
        + p9.geom_line()
        + p9.labs(x='Position (µm)', y='Absolute error (mM)', title=f"dx={dx}")
        + p9.scale_y_continuous(trans="log10")
    )
    g.save(f'3d-convergence-abs-error-dx-{dx}.png')

    g = (
        p9.ggplot(data[data["true_cai"] > 1e-8], p9.aes(x='x', y='rel_error'))
        + p9.geom_line()
        + p9.labs(x='Position (µm)', y='relative error', title=f"dx={dx}")
        + p9.scale_y_continuous(trans="log10")
    )
    g.save(f'3d-convergence-rel-error-dx-{dx}.png')

    return data

if __name__ == "__main__":
    with mp.Pool() as pool:
        all_data = pd.concat(pool.map(main, [0.25, 0.5, 0.125]))
    all_data["dx"] = all_data["dx"].astype("category")
    p9.options.figure_size = (4, 3)
    g = (
        p9.ggplot(all_data, p9.aes(x='x', y='abs_error', color='dx'))
        + p9.geom_line()
        + p9.scale_y_continuous(trans="log10")
        + p9.labs(x='Position (µm)', y='Absolute error (mM)')
    ).save("3d-convergence.pdf")

    pd.set_option('precision', 12)

    print("largest absolute errors:")
    print(all_data.groupby(["dx"])["abs_error"].max())
    
