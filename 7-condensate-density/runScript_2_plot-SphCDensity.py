import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats
import glob
import matplotlib.gridspec as gridspec
import ruptures as rpt

###plot ice slab density
df1= pd.read_csv('Chimera-generate-shellDensity.csv')

print(df1)
#selgroup        volseg  pixelRato2initalMask  sradius  numsel  ...  shelldistance  vol_density_numsel  vol_density_numsel_noise  filename sub_group_type idd

####################################### plot all
plt.gcf().subplots_adjust(left=0.25)
plt.gcf().subplots_adjust(bottom=0.35)
df=df1.loc[df1['sub_group_type'].isin(['nucleiA','dropS','dropH1'])]
print(df)

df['uniq']=df['sub_group_type']+'-'+df['idd'].astype(int).astype(str)+'-'+df['selgroup'].astype(int).astype(str) # get unique condensate cluster within each tomoset
print(df)
clusterinTomo = df['uniq'].drop_duplicates().tolist()
print(clusterinTomo)


df_align=pd.DataFrame(columns=df.columns) 
for uni in clusterinTomo:
    tempdf=pd.DataFrame()
    if uni in ['nucleiA-1-0', 'nucleiA-1-1', 'nucleiA-1-2', 'nucleiA-1-3', 'nucleiA-1-4', 'nucleiA-1-5', 'nucleiA-1-6','nucleiA-2-1', 'nucleiA-2-2', 'nucleiA-2-3', 'nucleiA-2-4','nucleiA-3-0', 'nucleiA-3-1', 'nucleiA-3-2', 'nucleiA-3-3', 'nucleiA-3-4', 'nucleiA-3-5','nucleiA-4-0', 'nucleiA-4-1', 'nucleiA-4-2']:
        tempdf=df.loc[df['uniq']==uni] # since there are multiple small nucleus within the same tomo set, split the condensates
    else:
        tempdf=df.loc[(df['uniq']==uni) & (df['selgroup']==0)] # for tomogram contain large condensate (only one droplet), select cluster 0


    if len(tempdf): # if select group contain cluster
        print(len(tempdf))


        ## select the smaller droplets which were not flatten by ice for analysis:
        if uni == 'dropH1-1-0':
            tempdf['sub_group_type']='NC-H1-L'
        if uni == 'dropH1-2-0':
            tempdf['sub_group_type']='NC-H1-S'
        if uni == 'dropH1-4-0':
            tempdf['sub_group_type']='NC-H1-S'
        if uni == 'dropH1-5-0':
            tempdf['sub_group_type']='NC-H1-L'
        if uni == 'dropH1-6-0':
            tempdf['sub_group_type']='NC-H1-S'



        if uni == 'dropS-1-0':
            tempdf['sub_group_type']='NC2-S'
        if uni == 'dropS-2-0':
            tempdf['sub_group_type']='NC2-S'
        if uni == 'dropS-3-0':
            tempdf['sub_group_type']='NC2-S'
        if uni == 'dropS-4-0':
            tempdf['sub_group_type']='NC2-S'
        if uni == 'dropS-5-0':
            tempdf['sub_group_type']='NC2-S'



        df_align=df_align.append(tempdf, ignore_index=True)

## remove the outlier (there is a density spike when the shell volume approach to the size of a NCP)
df_align = df_align[df_align.vol_density_numsel_Currentcluster < 1300]
df_align = df_align[df_align.shell_density_numsel_noisecluster < 1300]
df_align['sradius'] = df_align['sradius']/10 ## convert to nm unit
###############################################################################################################################################################################
#fig, axes = plt.subplots(nrows=3, ncols=3)
#fig.tight_layout(pad=2.0)
#fig.set_size_inches(15,15)
matplotlib.rcParams['legend.handlelength'] = 0
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

fig, axes = plt.subplots(nrows=3, ncols=3)
fig.set_size_inches(10,10)
fig.subplots_adjust(wspace=0.3, hspace=0.3) 
m=0
n=0
## loop through the sub_group_type
for subtype,ct2 in zip(['nucleiA','NC2-S','NC-H1-S'],[10,43,32]):
    subdata= df_align.loc[df_align['sub_group_type']==subtype]

    #print(subdata)
    ################################################## fill the plot position

    sns.lineplot(data=subdata, x="shell_idd", y="shell_density_numsel_Currentcluster", style="sub_group_type", markers=True,ax=axes[m,n],color=sns.color_palette("Set2").as_hex()[0])
    sns.lineplot(data=subdata, x="shell_idd", y="vol_density_numsel_Currentcluster", style="sub_group_type", markers=True,ax=axes[m,n],color=sns.color_palette("Set2").as_hex()[1])


    mean_vol2_conc=subdata.loc[(subdata['shell_idd']==ct2),'vol_density_numsel_Currentcluster' ].mean()
    axes[m,n].axhline(mean_vol2_conc, color='gray', linestyle='--',linewidth=2.5)


    axes[m,n].axvline(6, color='gray', linestyle='--',linewidth=2.5)
    axes[m,n].set_ylabel('Concentration (\u03BCM)',fontsize=14)
    axes[m,n].set_xlabel('Contour Number',fontsize=14)
    axes[m,n].set(xlim=(0,60))
    axes[m,n].set(ylim=(0,1200))
    axes[m,n].tick_params(labelsize=14)
    #axes[m,n].set_title( label=str(subtype),fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top


    x = np.arange(-100,100,0.1)


    axes[m,n].get_legend().remove()
    axes[0,1].yaxis.label.set_visible(False)
    axes[1,1].yaxis.label.set_visible(False)
    #ax2.get_legend().remove()


    m+=1
    if m==3:
        n+=1
        m=0

fig.savefig('SphC-density-avg.png', dpi=400)
plt.show()

