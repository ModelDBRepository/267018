import pandas as pd
import plotnine as p9
import sqlite3
import matplotlib.pyplot as plt

with sqlite3.connect("wave_time_3d.db") as conn:
    data = pd.read_sql("SELECT * FROM data", conn)

# it's possible there could be duplicates in the data; remove them
data = data.drop_duplicates(subset=["alpha", "dx", "phi", "theta"])



print(data["dx"].unique())
data["alpha"] = data["alpha"].astype("category")
data["relative_error"] *= 100
data035 = data[data["alpha"] == 0.25][(data["dx"] < 0.5) & (data["dx"] > 0.25)]
data["dx"] = data["dx"].astype("category")
data05 = data[data["alpha"] == 0.25][data["dx"] == 0.5]
data025 = data[data["alpha"] == 0.25][data["dx"] == 0.25]
my_data = data[data["alpha"] == 0.25]



print(f"We have a total of {len(data)} items.")

print(
    p9.ggplot(p9.aes(x="relative_error", fill="dx", color="dx")) + p9.geom_histogram(data05, binwidth=0.5, position="identity") + p9.geom_density(data05, alpha=0.1, size=2)
    + p9.geom_histogram(data025, binwidth=0.5, position="identity", alpha=0.7) + p9.geom_density(data025, alpha=0.1, size=2)
)
print(
    p9.ggplot(p9.aes(x="relative_error", color="dx")) + p9.geom_density(data05, alpha=0.1, size=2)
    + p9.geom_density(data025, alpha=0.1, size=2)
)

print(
    p9.ggplot(p9.aes(x="relative_error", color="dx", linetype="alpha")) + p9.geom_density(data[data["dx"] == 0.5], alpha=0.1, size=2)
    + p9.geom_density(data[data["dx"] == 0.25], alpha=0.1, size=2)
)


p9.options.figure_size = (3, 1.5)
(
    p9.ggplot(p9.aes(x="relative_error")) + p9.geom_histogram(data05, binwidth=0.5, position="identity")
    + p9.ggtitle("dx = 2 ^ -1")
    + p9.xlim(0, 8)
    + p9.ylim(0, 28)
).save("dx2-1.pdf")



(
    p9.ggplot(p9.aes(x="relative_error")) + p9.geom_histogram(data035, binwidth=0.5, position="identity")
    + p9.ggtitle("dx = 2 ^ -1.5")
    + p9.xlim(0, 8)
    + p9.ylim(0, 28)
).save("dx2-15.pdf")


(
    p9.ggplot(p9.aes(x="relative_error")) + p9.geom_histogram(data025, binwidth=0.5, position="identity")
    + p9.ggtitle("dx = 2 ^ -2")
    + p9.xlim(0, 8)
    + p9.ylim(0, 28)
).save("dx2-2.pdf")


print(f"for dx=0.5, relative error (%): {data05['relative_error'].mean()} ± {data05['relative_error'].std()}")
print(f"for dx=0.3535, relative error (%): {data035['relative_error'].mean()} ± {data035['relative_error'].std()}")
print(f"for dx=0.25, relative error (%): {data025['relative_error'].mean()} ± {data025['relative_error'].std()}")


plt.show()