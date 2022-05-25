import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats
import glob
import matplotlib.gridspec as gridspec


df1= pd.read_csv('Chimera-generate-shellDensity.csv')
print(df1)
#selgroup        volseg  pixelRato2initalMask  sradius  numsel  ...  shelldistance  vol_density_numsel  vol_density_numsel_noise  filename sub_group_type idd
####################################### plot all
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
    df_align=df_align.append(tempdf, ignore_index=True)


df_align.loc[df_align['sub_group_type']=='SP1','sub_group_type']='LS'
df_align.loc[df_align['sub_group_type']=='SP2','sub_group_type']='SP1'
df_align.loc[df_align['sub_group_type']=='SP3','sub_group_type']='SP2'

## remove the outlier (there is a density spike when the shell volume approach to the size of a NCP)
df_align = df_align[df_align.vol_density_numsel_Currentcluster < 1300]
df_align = df_align[df_align.shell_density_numsel_noisecluster < 1300]

df_align['sradius']=df_align['sradius']/10

#df_align.to_csv('zz-SP-peak-density-mask-info.csv', index=False) 


######################################################################### plot 
#fig, axes = plt.subplots(nrows=3, ncols=3)
#fig.tight_layout(pad=2.0)
#fig.set_size_inches(15,15)
matplotlib.rcParams['legend.handlelength'] = 0
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
fig, axes = plt.subplots(nrows=3, ncols=3)
fig.set_size_inches(10,10)
fig.subplots_adjust(wspace=0.3, hspace=0.3) # internal space

sns.set_palette("Set2")

for m,na,ct2 in zip([0,1,2],['LS','SP1','SP2'],[7,7,7]):

    sns.lineplot(data=df_align.loc[df_align['sub_group_type']==na], x="shell_idd", y="shell_density_numsel_Currentcluster", style="sub_group_type", markers=True,ax=axes[m,0],color=sns.color_palette("Set2").as_hex()[0])
    sns.lineplot(data=df_align.loc[df_align['sub_group_type']==na], x="shell_idd", y="vol_density_numsel_Currentcluster", style="sub_group_type", markers=True,ax=axes[m,0],color=sns.color_palette("Set2").as_hex()[1])



    mean_vol_conc=df_align.loc[(df_align['sub_group_type']==na) & (df_align['shell_idd']==ct2),'vol_density_numsel_Currentcluster' ].mean()


    axes[m,0].axhline(mean_vol_conc, color='gray', linestyle='--',linewidth=2.5)
    #axes[m,0].axhline(400, color='gray', linestyle='-.')
    axes[m,0].axvline(6, color='gray', linestyle='--',linewidth=2.5)
    axes[m,0].set_ylabel('Concentration (\u03BCM)',fontsize=14)
    axes[m,0].set_xlabel('Contour Number',fontsize=14)
    axes[m,0].set(xlim=(0,60))
    axes[m,0].set(ylim=(0,1200))
    axes[m,0].tick_params(labelsize=14)
    #axes[m,0].set_title( label=na,fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top

    x = np.arange(-100,100,0.1)
    y1=x*0 -10
    y2=x*0 +200
    y3=x*0 +400
    y4=x*0 +1300
    #axes[m,0].fill_between(x, y1,y2, facecolor='#f69431', alpha=0.2)
    #axes[m,0].fill_between(x, y2,y3, facecolor='#f1dd40', alpha=0.2)
    #axes[m,0].fill_between(x, y3,y4, facecolor='#c6c752', alpha=0.2)
    axes[m,0].get_legend().remove()
    #ax2.get_legend().remove()

fig.savefig('SpC-density-avg.png', dpi=400)

plt.show()




