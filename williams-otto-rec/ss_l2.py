import os, sys
import matplotlib.pyplot as plt, matplotlib as mpl, pandas as pd, numpy as np

# Formatting
plt.rcParams.update({
    'xtick.direction': 'out', 'ytick.direction': 'out',
    'xtick.top': False, 'ytick.right': False,
    'axes.spines.top': False, 'axes.spines.right': False,
})
figsize = (14, 6)
spine_offset = {'left': 6, 'bottom': 8}
y_fmt = '%.2e'                                        
x_i = ['x_1', 'x_2', 'x_3', 'x_4', 'x_5', 'x_6']
f_i = ['fb']
dx_i = ['dx_dt_1', 'dx_dt_2', 'dx_dt_3', 'dx_dt_4', 'dx_dt_5', 'dx_dt_6']
Tr_i = ['Tr']
groups = [(x_i, 'Blues_r'), (f_i, 'Greens_r'), (dx_i, 'Oranges_r'), (Tr_i, 'Purples_r')]
cmap_range = (0.15, 0.70)                             # _r maps: first run darkest, avoid near-white end                              

# Load
names = sys.argv[1:]
runs = {}
for name in names:
    csv = os.path.join(os.getcwd(), 'sims', name, 'wo.csv')
    runs[name] = pd.read_csv(csv)

# Calculate L2 norms
l2_by_run = {}
for run_name, df in runs.items():
    l2_by_run[run_name]={}
    for col_name, series in df.items():
        x = series.to_numpy(dtype=float)
        l2_by_run[run_name][col_name] = np.linalg.norm(x - x[0])

# Plot
fig, axes = plt.subplots(
    1, len(groups), figsize=figsize, layout='constrained',
    gridspec_kw={'width_ratios': [len(cols) + 0.8 for cols, _ in groups]},
)
shades = np.linspace(*cmap_range, len(list(l2_by_run)))
for ax, (cols, cmap) in zip(axes, groups):
    pos = np.arange(len(cols))
    for k, run in enumerate(list(l2_by_run)):
        vals = [l2_by_run[run][c] for c in cols]
        ax.bar(pos + (k - (len(list(l2_by_run)) - 1) / 2) * 0.75 / len(list(l2_by_run)), vals, 0.75 / len(list(l2_by_run)),
               facecolor=mpl.colormaps[cmap](shades[k]), edgecolor='black', linewidth=0.8)
    ax.set_xticks(pos, cols)
    ax.tick_params(axis='x', labelrotation=45)
    plt.setp(ax.get_xticklabels(), ha='right', rotation_mode='anchor')
    ax.set_xlim(-0.6, len(cols) - 0.4)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(5))
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter(y_fmt))
    for side, pts in spine_offset.items():
        ax.spines[side].set_position(('outward', pts))      

fig.supylabel(r'$\|z - z_0\|_2$')
swatches = [tuple(mpl.patches.Patch(facecolor=mpl.colormaps[cmap](s), edgecolor='black') for _, cmap in groups)
            for s in shades]
fig.legend(swatches, list(l2_by_run), loc='outside upper left', ncol=len(shades),
           handler_map={tuple: mpl.legend_handler.HandlerTuple(ndivide=None, pad=0)},
           handlelength=4)

fig.savefig(os.path.join(os.getcwd(), 'plots', '_'.join(names) + '.pdf'), dpi=900)
plt.show()
# Use
# python ss_l2.py name1 name2 ...
