import plotnine as p9
import pandas as pd
import sqlite3
import itertools

if __name__ == "__main__":
    conn = sqlite3.connect("discretization.db")
    data = pd.read_sql("SELECT * FROM morphology", conn)
    data = data.drop_duplicates(subset=["morphology", "dx"])
    data["length_ratio"] = data["sum_lengths"] / data["dx"]
    best_data = {}
    for row in data.itertuples():
        if row.morphology in best_data and row.dx >= best_data[row.morphology].dx:
            continue
        best_data[row.morphology] = row
    not_best_rows = [row.dx != best_data[row.morphology].dx for row in data.itertuples()]
    data["best_volume"] = [best_data[row.morphology].volume for row in data.itertuples()]
    data["best_surface_area"] = [best_data[row.morphology].surface_area for row in data.itertuples()]
    data["abs_volume_error"] = abs(data["volume"] - data["best_volume"])
    data["abs_surface_area_error"] = abs(data["surface_area"] - data["best_surface_area"])
    data["morphology"] = data["morphology"].apply(lambda filename: filename[4:-12])
    data["relative_volume_error"] = 100 * data["abs_volume_error"] / data["best_volume"]
    data["relative_surface_area_error"] = (
        100 * data["abs_surface_area_error"] / data["best_surface_area"]
    )
    
    p9.options.figure_size = (3.5, 3.5)

    # plot volume errors
    for y_var in [
        "abs_volume_error",
        "relative_volume_error",
        "abs_surface_area_error",
        "relative_surface_area_error",
    ]:
        for mode in ["", "_smooth"]:
            plot = (
                p9.ggplot(
                    data[not_best_rows], p9.aes(x="dx", y=y_var, color="morphology")
                )
                + p9.geom_point()
            )
            if mode == "":
                plot = plot + p9.geom_line()
            else:
                plot = plot + p9.geom_smooth(method="lm", se=False)


            plot = (
                plot
                + p9.scale_x_log10()
                + p9.scale_y_log10()
                + p9.xlab("dx (µm)")
            )
            if y_var == "relative_volume_error":
                plot = (
                    plot
                    + p9.geom_abline(slope=1, intercept=2, size=1)
                    + p9.geom_abline(slope=2, intercept=0.5, size=1)
                    + p9.geom_abline(slope=3, intercept=-1, size=1)
                )
            if y_var == "relative_surface_area_error":
                plot = (
                    plot
                    + p9.geom_abline(slope=1, intercept=3, size=1)
                    + p9.geom_abline(slope=2, intercept=1, size=1)
                    + p9.geom_abline(slope=3, intercept=-1, size=1)
                )
                        
            if y_var == "relative_volume_error":
                plot = plot + p9.ylab("Estimated Relative Volume Error (%)")
            elif y_var == "relative_surface_area_error":
                plot = plot + p9.ylab("Estimated Relative Surface Area Error (%)")

            plot.save(f"measurements/{y_var}{mode}.pdf")

    (
        p9.ggplot(
            data[not_best_rows], p9.aes(x="discretization_time", y="relative_volume_error", color="morphology")
        )
        #+ p9.geom_line()
        + p9.geom_point()
        + p9.scale_x_log10()
        + p9.scale_y_log10()
        + p9.geom_abline(slope=-1, intercept=3, size=1)
        + p9.geom_abline(slope=-2, intercept=-2, size=1)
        + p9.geom_smooth(method="lm", se=False)
        + p9.ylab("Estimated Relative Volume Error (%)")
        + p9.xlab("Discretization Time (s)")
    ).save(f"measurements/volume_error_vs_time.pdf")

    my_vars = [
        "volume",
        "surface_area",
        "num_voxels",
        "num_surface_voxels",
        "discretization_time",
        "length_ratio",
        "num_sections"
    ]

    # skipping num_sections because it does not depend on dx
    for y_var in set(my_vars) - {"num_sections"}:
        plot = (
            p9.ggplot(data, p9.aes(x="dx", y=y_var, color="morphology"))
            + p9.geom_point()
            + p9.geom_line()
            + p9.scale_x_log10()
            + p9.scale_y_log10()
            + p9.xlab("dx (µm)")
        )
        if y_var == "discretization_time":
            plot = (
                plot + p9.ylab("Discretization Time (s)")
                #+ p9.geom_abline(slope=-1, intercept=0, size=1)
                + p9.geom_abline(slope=-2, intercept=1, size=1)
                #+ p9.geom_abline(slope=-3, intercept=2, size=1)            
            )
        elif y_var == "relative_volume_error":
            plot = plot + p9.ylab("Estimated Relative Volume Error (%)")
        elif y_var == "relative_surface_area_error":
            plot = plot + p9.ylab("Estimated Relative Surface Area Error (%)")
        plot.save(f"measurements/{y_var}_vs_dx.pdf")
    

    data["dx"] = data["dx"].astype("category")

    for x_var, y_var in itertools.combinations(my_vars, 2):
        if x_var == "discretization_time":
            x_var, y_var = y_var, x_var
        plot = (
            p9.ggplot(data, p9.aes(x=x_var, y=y_var, color="morphology", shape="dx"))
            + p9.geom_point()
            + p9.scale_x_log10()
            + p9.scale_y_log10()
        )
        if y_var == "discretization_time":
            plot = plot + p9.ylab("Discretization Time (s)")
        
        plot.save(f"measurements/{y_var}_vs_{x_var}.pdf")
    
    for dx in data["dx"].unique():
        if dx != min(data["dx"]):
            data_at_dx = data[data["dx"] == dx]
            #print(data_at_dx)
            below_1 = len(data_at_dx[data_at_dx["relative_surface_area_error"] < 1])
            below_01 = len(data_at_dx[data_at_dx["relative_surface_area_error"] < 0.1])
            print(f"At dx={dx}, {below_1} out of {len(data_at_dx)} morphologies had a rel surface area error < 1%")
            print(f"At dx={dx}, {below_01} out of {len(data_at_dx)} morphologies had a rel surface area error < 0.1%")
            below_1 = len(data_at_dx[data_at_dx["relative_volume_error"] < 1])
            below_01 = len(data_at_dx[data_at_dx["relative_volume_error"] < 0.1])
            print(f"At dx={dx}, {below_1} out of {len(data_at_dx)} morphologies had a rel volume error < 1%")
            print(f"At dx={dx}, {below_01} out of {len(data_at_dx)} morphologies had a rel volume error < 0.1%")