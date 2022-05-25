
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats


import glob
readfilelist = glob.glob('Tomo_Gr.csv')

readfilelist.sort()
print(readfilelist)

for filename in readfilelist:
    df= pd.read_csv(filename)

    #make a 3x2 plot
    matplotlib.rcParams['legend.handlelength'] = 0
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    fig, axes = plt.subplots(nrows=3, ncols=2)
    fig.set_size_inches(10,10)
    fig.subplots_adjust(wspace=0.3, hspace=0.3) 

    m=0
    n=0



    types=df['sub_group_type'].drop_duplicates().tolist()
    for eachtype in types:
        print(eachtype)

        ########################################################  do the normalization of gr

        fig3 = plt.figure(2)  # temp figure to collect information of bar height
        ax3 = fig3.gca()
        barss=[h.get_height() for h in sns.barplot(data=df[(df['sub_group_type']==eachtype) & (df['r']<=600)], x="r", y="gr",ax=ax3,label=filename+eachtype).patches] # get bars height

        plt.close(2) # close the figure 2, for next iteration
        print(len(barss))
        print(len(df.loc[(df['sub_group_type']==eachtype) & (df['r']<=600) & (df['idd']==1),'r'].values))

        newdf=pd.DataFrame(columns=['r','gr'])
        newdf['r']=df.loc[(df['sub_group_type']==eachtype) & (df['r']<=600) & (df['idd']==1),'r'].values/10
        newdf['gr']=barss

        newdf['gr']=newdf['gr']/np.max(newdf['gr'].values)
        newdf.loc[newdf['r']==0.8,'gr']=0

        ######################################################### normalization done 



        print(newdf)


        sns.lineplot(data=newdf, x="r", y="gr",ax=axes[m,n],label=eachtype, linewidth = 3, color='black')
        axes[m,n].set(xlim=(0,60),ylim=(0,1.21)) 
        axes[m,n].axvline(8,ymax=1.2,color=sns.color_palette("Set2").as_hex()[2],ls='--',linewidth=2.5)
        axes[m,n].axvline(21,ymax=1.2,color=sns.color_palette("Set2").as_hex()[3],ls='--',linewidth=2.5)
        axes[m,n].set_ylabel('g(r)',fontsize=14)
        axes[m,n].set_xlabel('r (nm)',fontsize=14)
        #axes[m,n].xaxis.set_ticks([9,21,40,60])
        axes[m,n].yaxis.set_ticks(np.arange(0,1.2,0.3))
        axes[m,n].tick_params(labelsize=14)
        axes[m,n].get_legend().remove()
        m+=1
        if m==3:
            n+=1
            m=0


        #break

    fig.savefig("gr-plot.png",dpi=400)
plt.show()


