import sys

if len(sys.argv) < 5:
	print('Usage: python clustermap.py expression_data DEG out_table DEG_cluster.png gene_cluster.txt')
	exit(1)

from collections import OrderedDict
colors = OrderedDict({'C1': 'olive', 'C2': 'darkorange', 'C3': 'darkgreen', 'C4': 'deepskyblue', 'C5': 'purple',
		'C6': 'royalblue', 'C7': 'goldenrod', 'C8': 'chocolate',
		'Y': 'maroon', 'N': 'lightgrey'})
from matplotlib.colors import ListedColormap
tf_cmap = ListedColormap(['lightgrey', 'maroon'])
grp_cmap = list(colors.values())[:-2]
grp_cmap = ListedColormap(grp_cmap)

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

expr = pd.read_csv(sys.argv[1], index_col=0, sep='\t')
deg = pd.read_csv(sys.argv[2], sep='\t')
deg = list(deg.iloc[:,0].unique())
expr = expr.loc[deg]

name_color = {}
tf_color = {}
tf_col = []
cluster_order = list(colors.keys())
grp_col = []
with open(sys.argv[5]) as infile:
	infile.readline()
	for line in infile:
		gene, cluster, isTF = line.rstrip().split('\t')
		name_color[gene] = colors[cluster]
		tf_color[gene] = colors[isTF]
		grp_col.append(cluster_order.index(cluster))
		if isTF == 'Y':
			tf_col.append(1)
		else:
			tf_col.append(0)
row_color1 = [name_color[x] for x in expr.index]
row_color2 = [tf_color[x] for x in expr.index]
ans = sns.clustermap(expr, z_score=0, col_cluster=False, cmap='seismic', yticklabels=False,
		row_colors=[row_color1, row_color2],
		cbar_kws={"orientation": "horizontal"}, vmax=2.5, vmin=-2.5,
		figsize=(3.5, 8), dendrogram_ratio=0.15)
plt.subplots_adjust(bottom=0.1, top=0.95)
ans.ax_heatmap.set_xticklabels(ans.ax_heatmap.get_xmajorticklabels(), rotation=90)
ans.ax_heatmap.set_ylabel('')
reorder = [expr.index[idx] for idx in ans.dendrogram_row.reordered_ind]
ans.ax_heatmap.vlines([0, 5, 9], -1, len(deg), colors='white', linewidth=3)
ans.ax_cbar.set_position((0.22, 0.9, 0.4, 0.02))
for spine in ans.ax_cbar.spines:
    ans.ax_cbar.spines[spine].set_color('black')
    ans.ax_cbar.spines[spine].set_linewidth(1)
ans.ax_cbar.tick_params(labelsize=10)
ans.savefig(sys.argv[4], dpi=800)
plt.close()


expr = expr.loc[reorder]
scaler = StandardScaler()

temp = scaler.fit_transform(expr.values)
expr = pd.DataFrame(temp, index=expr.index, columns=expr.columns)
fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(3,8), sharey='row',
		gridspec_kw={'height_ratios':[0.025, 1], 'width_ratios': [0.07, 0.07, 1]})
ax[0,0].axis('off') # not use
ax[0,1].axis('off') # not use

grp_col = np.array([grp_col, grp_col]).T
ax[1,0].pcolor(grp_col, cmap=grp_cmap)
ax[1,0].set_xticks([])
ax[1,0].axis('off')

tf_col = np.array([tf_col, tf_col]).T
ax[1,1].pcolor(tf_col, cmap=tf_cmap)
ax[1,1].set_xticks([])
ax[1,1].axis('off')

sns.heatmap(expr, vmax=1, vmin=-1, robust=True, cmap="seismic", yticklabels=False, ax=ax[1,2], cbar=False)
ax[1,2].set_ylabel('')
ax[1,2].vlines([0, 5, 9], -1, len(deg), colors='white', linewidth=3)
# show colorbar of clustermap to ax[0,2]
fig.colorbar(ax[1,2].get_children()[0], cax=ax[0,2], orientation="horizontal")
fig.subplots_adjust(wspace=0.0, bottom=0.1, top=0.95)
ax[0,2].set_position( (0.22, 0.935, 0.5, 0.02) )

plt.savefig('heatmap.png', dpi=500)

