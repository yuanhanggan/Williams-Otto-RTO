import glob, math, os, sys
import matplotlib.pyplot as plt 
import pandas as pd 

# Load
names = sys.argv[1:]
runs = {}
for name in names: 
    csv = glob.glob(os.path.join(os.path.join(os.getcwd(), "sims"), name, "*.csv"))[0]
    runs[name] = pd.read_csv(csv, index_col=0)

# Plot
cols = runs[names[0]].columns
nrows = math.ceil(len(cols) / 2)
fig, axes = plt.subplots(nrows, 2, figsize = (8.27, 11.69), squeeze=False)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    for name, df in runs.items():
        ax.plot(df.index, df[col], label=name)
    ax.set_title(col)
    ax.grid(True)
for ax in axes[len(cols):]:
    ax.set_visible(False)
axes[0].legend()
fig.tight_layout()

# Save
os.makedirs(os.path.join(os.getcwd(), "plots", "_".join(names) + ".png"), dpi=600)
plt.show()

# Use
# python plot.py name1 name2 ...