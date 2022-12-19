import os
import sys
import pandas as pd
import sqlite3
base_dir = "swc"

dx = float(sys.argv[1])
for filename in os.listdir(base_dir):
    # TODO: skip if already in the db
    # skip hidden files
    with sqlite3.connect("discretization.db") as conn:
        old_data = pd.read_sql("SELECT * FROM morphology", conn)

    if any((old_data["dx"] == dx) & (old_data["morphology"] == os.path.join(base_dir, filename))):
        print(f"skipping: dx: {dx}, morph: {filename}")
        continue

    if not filename.startswith("."):
        true_filename = os.path.join(base_dir, filename)
        os.system(f"python time_discretization.py {true_filename} {dx}")
