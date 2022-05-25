# importing pandas package
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
import os
from matplotlib.ticker import FormatStrFormatter
from scipy.stats import pearsonr

#load previously generated measurement results
df1= pd.read_csv('../3-SpC-geometry/SpGeometry.csv')
print(df1)
df2= pd.read_csv('../4-SpC-Surr-FreeNCP-density/condensateEdge_FreeNCP_Density.csv')
print(df2)

# format df2 as df1
df2clean=pd.DataFrame(columns=df2.columns)
print(df2clean)
unique_vals = df2['uniq'].drop_duplicates().tolist()
for val in unique_vals:
    temserie=df2[df2['uniq']==val].mean(axis = 0)
    temserie['sub_group_type']=df2[df2['uniq']==val].iloc[0]['sub_group_type']
    df2clean=df2clean.append(temserie, ignore_index=True)
print(df2clean[['selgroup','sub_group_type','idd','shell_density_numsel_noisecluster']])

# merge df2 to df1
df1['Surrounding Noise'] = df2clean['shell_density_numsel_noisecluster']

print(df1)
df3=df1[['sub_group_type','idd','eccentricity','Avg Size (nm)','Surrounding Noise','numPar']]
df3['idd2']=df3['sub_group_type']+df3['idd'].astype('str')

print(df3)
df3clean=pd.DataFrame(columns=['sub_group_type', 'idd', 'eccentricity', 'Avg Size (nm)','numPar',
        'Surrounding Noise', 'idd2','eccentricity_sd','Avg Size (nm)_sd','Surrounding Noise_sd'])
unique_vals = df3['idd2'].drop_duplicates().tolist()
for val in unique_vals:
    temserieM=df3[df3['idd2']==val].mean(axis = 0)
    temserieSD=df3[df3['idd2']==val].std(axis = 0)
    #first copy the mean
    temserie=temserieM
    temserie['sub_group_type']=df3[df3['idd2']==val].iloc[0]['sub_group_type']
    temserie['idd2']=df3[df3['idd2']==val].iloc[0]['idd2']
    ## next copy the sd
    temserie['eccentricity_sd']=temserieSD['eccentricity']
    temserie['Avg Size (nm)_sd']=temserieSD['Avg Size (nm)']
    temserie['Surrounding Noise_sd']=temserieSD['Surrounding Noise']
    temserie['numPar_sd']=temserieSD['numPar']

    #print(temserie)

    df3clean=df3clean.append(temserie, ignore_index=True)

df4=df3clean.drop(columns=['idd'])
print(df4)


################

dfNeighbour_dist=pd.read_csv('../5-SpC-neighbour-Dist/SpC-neighborDist.csv')
dfNeighbour_dist['uniq']=dfNeighbour_dist['sub_group_type'] + dfNeighbour_dist['idd'].astype(str)
#print(dfNeighbour_dist)
df5=dfNeighbour_dist.groupby('uniq').mean().reset_index()
print(df5)


### merge with df4 
df4['Neighbour Dist']=df5['Distance']
print(df4)

df4.columns=['sub_group_type', 'Eccentricity', 'Avg Size (nm)', 'NCP Number',
       'Free Par Conc. (\u03BCM)', 'idd2', 'eccentricity_sd', 'Avg Size (nm)_sd',
       'Surrounding Noise_sd', 'numPar_sd', 'Neighbour Dist. (nm)']

###############################################################################
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, axes = plt.subplots(nrows=5, ncols=5)
fig.set_size_inches(12,12)
plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.12,wspace=0.1, hspace=0.1)
#fig.tight_layout(pad=1)

#set color pallete
sns.set_palette(sns.color_palette(['#efd9c1','#c7d8c6','#a9b7c0']))
sns.set_palette(sns.color_palette("Set2"))

alist=pd.DataFrame(columns=['q1','q2','q3','q4'])
#n=0
#for nameA in ['Eccentricity','Free Par Conc. (\u03BCM)','Avg Size (nm)', 'NCP Number','Neighbour Dist. (nm)']:
#    m=0
#    for nameB in  ['Eccentricity','Free Par Conc. (\u03BCM)','Avg Size (nm)', 'NCP Number','Neighbour Dist. (nm)']:

namelist=['Eccentricity','Free Par Conc. (\u03BCM)','Avg Size (nm)', 'NCP Number','Neighbour Dist. (nm)']
xrangelist=[[0.5,0.7],[1,5],[50,80],[50,150],[40,60,80,100]]
for n in range(0,5):
    for m in range(1+n,5):
        nameA=namelist[n]
        nameB=namelist[m]
        xaxisrange=xrangelist[n]

        sns.scatterplot(x=nameA, y=nameB, data=df4, hue="sub_group_type",ax=axes[m,n],s=50)
        slope, intercept, r_value, p_value, std_err = stats.linregress(df4[nameA].values,df4[nameB].values)

        if r_value <0:
            axes[m,n].text(0.9, 0.92, "r="+"{:.1f}".format(r_value), ha="right", va="top", transform=axes[m,n].transAxes,fontsize=14)
        else:
            axes[m,n].text(0.1, 0.92, "r="+"{:.1f}".format(r_value), ha="left", va="top", transform=axes[m,n].transAxes,fontsize=14)

        print(nameA,nameB)
        #axes[m,n].set_xlabel(nameA,fontsize=14)
        #axes[m,n].set_ylabel(nameB,fontsize=14)
    #    axes[m,n].set(xlim=(0.8,1))
    #    axes[m,n].set(ylim=(20,90))
        axes[m,n].tick_params(labelsize=14)
        axes[m,n].get_legend().remove()
        # plot the regression line on the extended canvas
        xlims = axes[m,n].get_xlim()
        new_x = np.arange(xlims[0], xlims[1],(xlims[1]-xlims[0])/250.)
        axes[m,n].plot(new_x, intercept + slope *  new_x, color='gray', linestyle='--', lw = 2.5)
        axes[m,n].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axes[m,n].tick_params(labelsize=14)
        axes[m,n].set_xticks(xaxisrange) 

        axes[m,n].set_xlabel(namelist[n], fontsize=12)
        axes[m,n].set_ylabel(namelist[m], fontsize=12)


        if m < 4:
           axes[m,n].xaxis.set_visible(False)
        if n > 0:
           axes[m,n].yaxis.set_visible(False)



for ia in range(0,5):
    for ib in range(0+ia,5):
            axes[ia,ib].axis('off')






fig.savefig("Correlation-Plot.png",dpi=400)

plt.show()

