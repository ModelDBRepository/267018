import sqlite3
import pandas as pd
import plotnine as p9
import matplotlib.pyplot as plt

DB_FILENAME = "simple_geometry_convergence.db"


def plot_it():
    with sqlite3.connect(DB_FILENAME) as conn:
        data = pd.read_sql("SELECT * FROM data", conn)
    data["Refinements"] = data["resolution"].astype("category")
    for error_field in ["surface_area_relative_error", "volume_relative_error"]:
        data[error_field] = data[error_field].abs()
    original_data = data
    data = pd.melt(
        data,
        id_vars=["dx", "Refinements"],
        value_vars=["surface_area_relative_error", "volume_relative_error"],
        value_name="Relative Error",
        var_name="Measurement",
    )
    data["Measurement"] = data["Measurement"].apply(
        lambda name: {
            "surface_area_relative_error": "Surface Area",
            "volume_relative_error": "Volume",
        }[name]
    )
    print(
        p9.ggplot(data, p9.aes(color="Refinements", x="dx", y="Relative Error"))#, linetype="Measurement"))
        + p9.geom_line(size=1.2)
        + p9.scale_x_log10()
        + p9.scale_y_log10()
        + p9.facet_wrap("Measurement")
        + p9.theme(subplots_adjust={'right': 0.7})
    )
    print(
        p9.ggplot(original_data, p9.aes(color="Refinements", x="runtime", y="volume_relative_error"))
        + p9.geom_line(size=1.2)
        + p9.geom_point()
        + p9.scale_x_log10()
        + p9.scale_y_log10()
        + p9.theme(subplots_adjust={'right': 0.8})
    )
    import matplotlib.pyplot as plt
    fig = plt.figure()
    axes = []
    all_res = [2, 4, 6, 8, 10]
    colors = ["blue", "red", "darkgreen", "brown", "black"]
    for i, (res, color) in enumerate(zip(all_res, colors)):
        axis = fig.add_subplot(5, 1, i + 1)
        axis.set_xscale("log")
        axis.set_yscale("log")
        for j, res2 in enumerate(all_res):
            my_data = original_data[original_data["resolution"] == res2].sort_values(by="runtime")
            if j != i:
                axis.plot(my_data["runtime"], my_data["volume_relative_error"], color='lightgray')
        my_data = original_data[original_data["resolution"] == res].sort_values(by="runtime")
        axis.plot(my_data["runtime"], my_data["volume_relative_error"], color=color, label=f"{res} VRs")
        axis.legend()
        if i == 2:
            axis.set_ylabel("Relative volume error")
        if i == len(all_res) - 1:
            axis.set_xlabel("Runtime (s)")
        else:
            axis.set_xticks([])
        axis.set_yticks([1e-5, 1e-3, 1e-1])
        #axis.set_ylabel(f"{res} VRs")
        
        axes.append(axis)
    plt.show()


if __name__ == "__main__":
    plot_it()
