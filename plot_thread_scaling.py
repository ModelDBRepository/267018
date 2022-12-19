import sqlite3
import pandas as pd
import plotnine as p9

DB_FILENAME = "thread_scaling.db"

with sqlite3.connect(DB_FILENAME) as conn:
    data = pd.read_sql("""
        SELECT nthread, morphology, kinetics, dx, MIN(runtime)
        FROM data
        GROUP BY nthread, morphology, kinetics, dx
    """, conn)

data["dx"] = data["dx"].astype("category")
data["Model"] = [f"{morph} {kinetics}" for morph, kinetics in zip(data["morphology"], data["kinetics"])]

print(
    data[data["Model"] == "cell cawave"]
)

print(
    p9.ggplot(data, p9.aes(x="nthread", y="MIN(runtime)", color="Model", linetype="dx"))
    + p9.geom_line(size=1.2)
    + p9.geom_point()
    + p9.scale_x_continuous(trans="log2")
    + p9.scale_y_continuous(trans="log10")
    + p9.theme(subplots_adjust={'right': 0.7})
    + p9.labs(x="Number of threads", y="Minimum simulation time (s)")
)