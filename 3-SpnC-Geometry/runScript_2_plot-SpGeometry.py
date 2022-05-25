
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
import os
from matplotlib.ticker import FormatStrFormatter

df= pd.read_csv('sp-statistics.csv')

df['Avg Size (nm)']=(df['PC1-axis'] + df['PC3-axis'])/2/10
dfsp= df.loc[df['sub_group_type'].isin(['SP1','SP2','SP3'])]

print(dfsp)


# remove the group label with max id (noise group)
unique_vals = dfsp['tomo_set'].drop_duplicates().tolist()
for val in unique_vals:
    noiseClusterid=max(dfsp.loc[dfsp['tomo_set']==val,'clustersubgroup'])
    dfsp = dfsp.drop(dfsp[(dfsp.tomo_set==val) & (df.clustersubgroup==noiseClusterid)].index)


temptranpose1=dfsp[['PC1-axis','idd','sub_group_type']]
temptranpose1.columns=['pcval','type','sub_group_type']
temptranpose1['type']='PC1'


temptranpose3=dfsp[['PC3-axis','idd','sub_group_type']]
temptranpose3.columns=['pcval','type','sub_group_type']
temptranpose3['type']='PC3'
print(temptranpose3)

tranposedf=pd.concat([temptranpose1,temptranpose3])
print(tranposedf)
tranposedf['pcval']=tranposedf['pcval']/10
print(tranposedf)


dfsp.loc[dfsp['sub_group_type']=='SP1','sub_group_type']='LS'
dfsp.loc[dfsp['sub_group_type']=='SP2','sub_group_type']='SP1'
dfsp.loc[dfsp['sub_group_type']=='SP3','sub_group_type']='SP2'

tranposedf.loc[tranposedf['sub_group_type']=='SP1','sub_group_type']='LS'
tranposedf.loc[tranposedf['sub_group_type']=='SP2','sub_group_type']='SP1'
tranposedf.loc[tranposedf['sub_group_type']=='SP3','sub_group_type']='SP2'

dfsp.to_csv('SpGeometry.csv', index=False)


############################################### plot figure
fig, axes = plt.subplots(nrows=2, ncols=3)
fig.tight_layout(pad=3.5)
fig.set_size_inches(15,10)
#plt.rcParams["font.weight"] = "bold"
#plt.rcParams["axes.labelweight"] = "bold"
#matplotlib.rcParams.update({'font.size': 22})


### leave the board 
#plt.gcf().subplots_adjust(left=0.25)
#plt.gcf().subplots_adjust(bottom=0.35)
sns.set_palette(sns.color_palette("Set2"))

sns.barplot(x="sub_group_type", y="Avg Size (nm)", data=dfsp, capsize=.2,ax=axes[0,0])
sns.barplot(x="sub_group_type", y="eccentricity", data=dfsp, capsize=.2,ax=axes[0,1])
sns.barplot(x="type", y="pcval", data=tranposedf,hue="sub_group_type", capsize=.2,ax=axes[0,2])
sns.barplot(x="sub_group_type", y="numPar", data=dfsp, capsize=.2,ax=axes[1,2])




def report_plot_stat(axx):
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    for p in axx.lines:
        print(p)
        xy = p.get_xydata()
        ## get the value of vertical error bar
        if not xy[0][1]-xy[1][1] ==0:
            print('-----' + str((abs(xy[0][1]-xy[1][1])/2)))
    print('--------') 
    for q in axx.patches:
        height = q.get_height()
        print(height)

report_plot_stat(axes[0,0])
report_plot_stat(axes[0,1])
report_plot_stat(axes[0,2])
report_plot_stat(axes[1,2])



axes[0,0].set_ylabel('Length (nm)',fontsize=14)
axes[0,0].set_xlabel('Condensate Type',fontsize=14)
axes[0,0].set(ylim=(0,170))
axes[0,0].tick_params(labelsize=14)
axes[0,0].set_title( label='Cluster Average Size',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top

axes[0,1].set_ylabel('Eccentricity',fontsize=14)
axes[0,1].set_xlabel('Condensate Type',fontsize=14)
axes[0,1].set(ylim=(0.4,1.1))
axes[0,1].tick_params(labelsize=14)
axes[0,1].set_title( label='Cluster Eccentricity',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top
axes[0,1].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))


axes[0,2].set_ylabel('Length (nm)',fontsize=14)
axes[0,2].set_xlabel('PC Axis',fontsize=14)
axes[0,2].set(ylim=(0,300))
axes[0,2].tick_params(labelsize=14)
axes[0,2].set_title( label='Cluster Dimensions',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top
axes[0,2].get_legend().remove()

axes[1,2].set_ylabel('Count',fontsize=14)
axes[1,2].set_xlabel('Condensate Type',fontsize=14)
axes[1,2].set(ylim=(0,200))
axes[1,2].tick_params(labelsize=14)
axes[1,2].set_title( label='NCP Number',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top



plt.savefig("distPlot.png",dpi=400)
plt.show()


