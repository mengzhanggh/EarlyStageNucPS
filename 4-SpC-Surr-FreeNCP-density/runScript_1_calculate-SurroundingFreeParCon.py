import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats
import glob
import matplotlib.gridspec as gridspec


df1= pd.read_csv('Chimera-generate-FreeNCP-Density.csv')
plt.gcf().subplots_adjust(left=0.25)
plt.gcf().subplots_adjust(bottom=0.35)
df=df1.loc[df1['sub_group_type'].isin(['SP1','SP2','SP3'])]
print(df)

df['uniq']=df['sub_group_type']+'-'+df['idd'].astype(int).astype(str)+'-'+df['selgroup'].astype(int).astype(str)
print(df)
clusterinTomo = df['uniq'].drop_duplicates().tolist()
print(clusterinTomo)
#'SP1-1-0', 'SP1-1-1', 'SP1-1-2', 'SP1-1-3', 'SP1-1-4', 


df_align=pd.DataFrame(columns=df.columns) 
for uni in clusterinTomo:
    tempdf=df.loc[df['uniq']==uni]
    ## Manually find the pixelRato2initalMask that close to 1 (here, zero shell here defined as the outer most shell, where FreeNCP density start to increase)
    manualdenfined_zero_shell =tempdf.iloc[(tempdf['pixelRato2initalMask']-1).abs().argsort()[:1],]['sradius'].values
    ## Align the data to the zero shell
    if tempdf.iloc[0]['sub_group_type']=='SP1':
        tempdf['sradius']=tempdf['sradius']- manualdenfined_zero_shell
    if tempdf.iloc[0]['sub_group_type']=='SP2':
        tempdf['sradius']=tempdf['sradius']- manualdenfined_zero_shell+40
    if tempdf.iloc[0]['sub_group_type']=='SP3':
        tempdf['sradius']=tempdf['sradius']- manualdenfined_zero_shell

    df_align=df_align.append(tempdf, ignore_index=True)


df_align.loc[df_align['sub_group_type']=='SP1','sub_group_type']='LS'
df_align.loc[df_align['sub_group_type']=='SP2','sub_group_type']='SP1'
df_align.loc[df_align['sub_group_type']=='SP3','sub_group_type']='SP2'

df_align['sradius']=df_align['sradius']/10 ##A to nm conversion


fig, axes = plt.subplots(nrows=2, ncols=3)
fig.tight_layout(pad=2.0)
fig.set_size_inches(15,10)

#########################################################################

fig, axes = plt.subplots(nrows=2, ncols=3)
fig.tight_layout(pad=3.5)
fig.set_size_inches(15,10)
#set  color pallete
sns.set_palette(sns.color_palette(['#efd9c1','#c7d8c6','#a9b7c0']))
sns.set_palette(sns.color_palette("Set2"))

for m,na in zip([0,1,2],['LS','SP1','SP2']):
    # twin objects for two different y-axis on the sample plot

    # make a plot with different y-axis using second axis object
    sns.lineplot(data=df_align.loc[df_align['sub_group_type']==na], x="sradius", y="shell_density_numsel_noisecluster",err_style="bars",err_kws={'capsize':3},style="sub_group_type", markers=True, dashes=False,ax=axes[0,m],color='black')





    #axes[0,m].axvline(-60, color='gray', linestyle='--')
    axes[0,m].set_ylabel('Condensate Conc. Gradient (\u03BCM)',color="#396d9a",fontsize=14)
    axes[0,m].set_xlabel('Delta Mask Extension (nm)',fontsize=14)
    axes[0,m].set(xlim=(-1,22))
    axes[0,m].set(ylim=(0,20))
    axes[0,m].tick_params(labelsize=14)
    axes[0,m].set_title( label=na,fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top

    axes[0,m].get_legend().remove()



condensateEdge_FreeNCP= df_align.loc[(df_align['sradius']>0) & (df_align['sradius']<20)]
condensateEdge_FreeNCP.to_csv('condensateEdge_FreeNCP_Density.csv', index=False)




sns.barplot(x="sub_group_type", y="shell_density_numsel_noisecluster", data=condensateEdge_FreeNCP, capsize=.2,ax=axes[1,0])
axes[1,0].set_ylabel('Concentration (\u03BCM)',fontsize=14)
axes[1,0].set_xlabel('Condensate Type',fontsize=14)
axes[1,0].set(ylim=(0,10))
axes[1,0].tick_params(labelsize=14)
axes[1,0].set_title( label='Surrounding Free Par Conc.',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top



sns.barplot(x="sub_group_type", y="shell_density_numsel_noisecluster", data=condensateEdge_FreeNCP, capsize=.2,ax=axes[1,1])
axes[1,0].set_ylabel('Concentration (\u03BCM)',fontsize=14)
axes[1,0].set_xlabel('Condensate Type',fontsize=14)
#axes[1,0].set(ylim=(0,7))
axes[1,0].tick_params(labelsize=14)
axes[1,0].set_title( label='Surrounding Free Par Conc.',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top



from statannot import add_stat_annotation
add_stat_annotation(axes[1,1], data=condensateEdge_FreeNCP, x="sub_group_type", y="shell_density_numsel_noisecluster",
                    box_pairs=[("LS", "SP1"), ("SP1", "SP2"), ("LS", "SP2")],
                    test='t-test_ind', text_format='star', loc='inside', line_offset_to_box=-1, line_offset=0.05, line_height=0.1, verbose=2)




fig.savefig('SpC-Surr-FreeNCP-density.png', dpi=400)

def report_plot_stat(axx):
    print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
    for p in axx.lines:
        print(p)
        xy = p.get_xydata()
        ## if the y value is not equal, then it is the vertical error bar
        if not xy[0][1]-xy[1][1] ==0:
            print('-----' + str((abs(xy[0][1]-xy[1][1])/2)))
    print('--------') 
    for q in axx.patches:
        height = q.get_height()
        print(height)

report_plot_stat(axes[1,0])




plt.show()
