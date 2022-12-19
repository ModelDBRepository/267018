import sqlite3
import pandas as pd
from math import pi
import plotnine as p9

with sqlite3.connect("cylinder_convergence.db") as conn:
    data = pd.read_sql("SELECT * FROM data", conn)

true_volume = pi * 5
true_surface_area = pi * 2 * 5 + 2 * pi

data["rel_volume_error"] = 100 * (data["volume"] - true_volume) / true_volume
data["rel_surface_error"] = 100 * (data["area"] - true_surface_area) / true_surface_area

print(data)

for dx in data["dx"].unique():
    my_data = data[data["dx"] == dx]
    (
        p9.ggplot(p9.aes(x="rel_volume_error")) + p9.geom_histogram(my_data, position="identity")
        + p9.ggtitle(f"dx = {dx}")
        + p9.xlab("Relative Volume Error (%)")
    ).save(f"cylinder_convergence_volume_error_{dx}.pdf")

    (
        p9.ggplot(p9.aes(x="rel_surface_error")) + p9.geom_histogram(my_data, position="identity")
        + p9.ggtitle(f"dx = {dx}")
        + p9.xlab("Relative Surface Area Error (%)")
    ).save(f"cylinder_convergence_surface_error_{dx}.pdf")

for dx in data["dx"].unique():
    my_data = data[data["dx"] == dx]
    print(f"At dx={dx}, average % volume error = {my_data['rel_volume_error'].abs().mean()} ± {my_data['rel_volume_error'].abs().std()}")
    print(f"At dx={dx}, average % surface error = {my_data['rel_surface_error'].abs().mean()} ± {my_data['rel_surface_error'].abs().std()}")

