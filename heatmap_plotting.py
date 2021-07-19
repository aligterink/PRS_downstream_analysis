# This is a not very elegant script for plotting PRS-phenotype associations. 

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cm

### Settings
toolkit_output = R"C:\Users\anton\Documents\UMCU\Joost_data\OLINK_SlideToolkit\OUTPUT\metagrs_comp_proteins.xlsx" # Output file of the Slide toolkit
r2_root = True # True: plot will be made using the root of the r2 value, False: make plot using regular r2
beta_fontsize = 6 # Fontsize of the beta value in the plot
pvalue_thres1 = 0.0167 # P-value at which to mark an association

# Read the Slide toolkit associations
df = pd.read_excel(toolkit_output, header=0)

# Take the square root of r2
if r2_root:
    df['r^2']=df['r^2']**(1/2)
    ylabel = "√r2"
else: 
    ylabel = "r2"

# Remove _rankNorm from protein names
df['Trait'] = df['Trait'].replace(regex=['_rankNorm'], value='')

df = df.sort_values(by='Trait')
df = df.sort_values(by='Predictor')
rownames = df['Trait'].unique()
colnames = df['Predictor'].unique()

# Create empty dataframes
r2_df = pd.DataFrame(pd.np.empty((0, len(colnames))))
pv_df = pd.DataFrame(pd.np.empty((0, len(colnames))))
beta_df = pd.DataFrame(pd.np.empty((0, len(colnames))))

# Set column names of all three dataframes
r2_df.columns = colnames
pv_df.columns = colnames
beta_df.columns = colnames

# Fill dataframes
for trait in rownames:
    trait_rows = df.loc[df['Trait'] == trait]
    trait_rows = trait_rows.sort_values(by='Predictor')

    r2_df.loc[len(r2_df.index)] = trait_rows['r^2'].tolist()
    pv_df.loc[len(pv_df.index)] = trait_rows['P-value'].tolist()
    beta_df.loc[len(beta_df.index)] = trait_rows['Beta'].tolist()

# Set rownames of all three dataframes
r2_df = r2_df.rename(index = lambda x: rownames[x])
pv_df = pv_df.rename(index = lambda x: rownames[x])
beta_df = beta_df.rename(index = lambda x: rownames[x])    

# Get min and max r2 values
minval = r2_df.min().min()
maxval = r2_df.max().max()

# Sort proteins by combined r2
r2_df['sum'] = r2_df.sum(axis=1).to_list()
r2_df = r2_df.sort_values(by='sum', ascending=False)
r2_df = r2_df.drop(['sum'], axis=1)

# Set the proteins that occur in each of the 3 sub-heatmaps
r2_df1 = r2_df.iloc[0:88, ]
r2_df2 = r2_df.iloc[88:176, ]
r2_df3 = r2_df.iloc[176:, ]

rownames_df1 = r2_df1.index.values
rownames_df2 = r2_df2.index.values
rownames_df3 = r2_df3.index.values

fig = plt.figure(figsize=(17, 35))
# fig = plt.figure(figsize=(10, 30))

ax1 = fig.add_axes([0.1, 0.1, 0.3, 0.8])
ax2 = fig.add_axes([0.4, 0.1, 0.3, 0.8])
ax3 = fig.add_axes([0.7, 0.1, 0.3, 0.8])

im = ax1.imshow(r2_df1, cmap=cm.plasma, vmin=minval, vmax=maxval)
im = ax2.imshow(r2_df2, cmap=cm.plasma, vmin=minval, vmax=maxval)
im = ax3.imshow(r2_df3, cmap=cm.plasma, vmin=minval, vmax=maxval)

ax1.set_aspect(0.5)
ax2.set_aspect(0.5)
ax3.set_aspect(0.5)

# Setting the labels
ax1.set_xticks(np.arange(len(colnames)))
ax1.set_yticks(np.arange(len(rownames_df1)))
ax1.set_xticklabels(colnames)
ax1.set_yticklabels(rownames_df1)

ax2.set_xticks(np.arange(len(colnames)))
ax2.set_yticks(np.arange(len(rownames_df2)))
ax2.set_xticklabels(colnames)
ax2.set_yticklabels(rownames_df2)

ax3.set_xticks(np.arange(len(colnames)))
ax3.set_yticks(np.arange(len(rownames_df3)))
ax3.set_xticklabels(colnames)
ax3.set_yticklabels(rownames_df3)

# Rotate the tick labels and set their alignment.
plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.setp(ax3.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")

### Add beta value to each association
for i in range(len(rownames_df1)):
    for j in range(len(colnames)):
        pvalue = pv_df.loc[rownames_df1[i], colnames[j]]
        if pvalue <= pvalue_thres1:
            ax1.text(j, i, round(beta_df.loc[rownames_df1[i], colnames[j]], 2), ha="center", va="center", color="green", fontsize=beta_fontsize, weight='bold')
        else:
            ax1.text(j, i, round(beta_df.loc[rownames_df1[i], colnames[j]], 2), ha="center", va="center", color="w", fontsize=beta_fontsize)

for i in range(len(rownames_df2)):
    for j in range(len(colnames)):
        pvalue = pv_df.loc[rownames_df2[i], colnames[j]]
        if pvalue <= pvalue_thres1:
            ax2.text(j, i, round(beta_df.loc[rownames_df2[i], colnames[j]], 2), ha="center", va="center", color="green", fontsize=beta_fontsize, weight='bold')
        else:
            ax2.text(j, i, round(beta_df.loc[rownames_df2[i], colnames[j]], 2), ha="center", va="center", color="w", fontsize=beta_fontsize)

for i in range(len(rownames_df3)):
    for j in range(len(colnames)):
        pvalue = pv_df.loc[rownames_df3[i], colnames[j]]
        if pvalue <= pvalue_thres1:
            ax3.text(j, i, round(beta_df.loc[rownames_df3[i], colnames[j]], 2), ha="center", va="center", color="green", fontsize=beta_fontsize, weight='bold')
        else:
            ax3.text(j, i, round(beta_df.loc[rownames_df3[i], colnames[j]], 2), ha="center", va="center", color="w", fontsize=beta_fontsize)

plt.colorbar(im, ax=(ax1, ax2, ax3), label=ylabel)
plt.tight_layout()
plt.show()
