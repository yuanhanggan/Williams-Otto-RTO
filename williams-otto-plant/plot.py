import glob, math, os, sys
import matplotlib.pyplot as plt 
import pandas as pd 
plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.formatter.use_mathtext": True,   
    "font.size": 10,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.major.size": 4, "ytick.major.size": 4,
})

# Load
names = sys.argv[1:]
runs = {}
for name in names: 
    csv = glob.glob(os.path.join(os.path.join(os.getcwd(), "sims"), name, "*.csv"))[0]
    runs[name] = pd.read_csv(csv, index_col=0)

# Plot
cols = runs[names[0]].columns
nrows = math.ceil(len(cols) / 2)
fig, axes = plt.subplots(nrows, 2, figsize = (10, 3 * nrows), squeeze=False)
axes = axes.flatten()
for ax, col in zip(axes, cols):
    for name, df in runs.items():
        ax.plot(df.index, df[col], label=name)
    ax.set_title(col)
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
for ax in axes[len(cols):]:
    ax.set_visible(False)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', ncol=len(names), frameon=False)
fig.tight_layout(h_pad=3, rect=[0, 0.5 / fig.get_figheight(), 1, 1])

# Save
os.makedirs(os.path.join(os.getcwd(), "plots"), exist_ok=True)
fig.savefig(os.path.join(os.path.join(os.getcwd(), "plots"), "_".join(names) + ".pdf"), dpi=600)
plt.show()

# Use
# python plot.py name1 name2 ...