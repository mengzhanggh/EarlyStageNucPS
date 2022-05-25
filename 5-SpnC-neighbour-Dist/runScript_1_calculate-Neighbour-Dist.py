import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats
import glob
import matplotlib.gridspec as gridspec
from sklearn.neighbors import NearestNeighbors
from matplotlib import collections  as mc




df_dist=pd.DataFrame(columns=['Distance','sub_group_type','idd'])
df1= pd.read_table('data-setList.txt', sep='\s*',header=None)
df1.columns=['sub_group_type','idd','Center_slice_id','lpn']
#df1.loc[df1['sub_group_type']=='SP1','sub_group_type']='LS'
#df1.loc[df1['sub_group_type']=='SP2','sub_group_type']='SP1'
#df1.loc[df1['sub_group_type']=='SP3','sub_group_type']='SP2'




print(df1)
df1['tomo_set']=df1['sub_group_type']+'-'+df1['idd'].astype(str)
for prefix in df1['tomo_set'].drop_duplicates().tolist():
    print(prefix)

    
    tempdf1= pd.read_table(prefix + '-area.txt', sep='\s*',header=None)
    tempdf1.columns=['contourID','pt_num','area']
    tempdf2= pd.read_table(prefix + '-coord.txt', sep='\s*',header=None)
    tempdf2.columns=['contourID','x','y','n']


    print(tempdf1['area'].mean())
    print(tempdf1['area'].std())
    tempdf1['x_c']=0
    tempdf1['y_c']=0

    ### for each contour, calc weight center
    for index, row in tempdf1.iterrows():
        tempdf1['x_c'].iloc[index] = tempdf2.loc[tempdf2['contourID']==row['contourID'], 'x'].mean()
        tempdf1['y_c'].iloc[index] = tempdf2.loc[tempdf2['contourID']==row['contourID'], 'y'].mean()

    ## convert tempdf1 into nparray style

    X = tempdf1[['x_c','y_c']].values
    nbrs = NearestNeighbors(n_neighbors=2, algorithm='ball_tree').fit(X)
    distances, indices = nbrs.kneighbors(X)
    print(indices)
    print(distances)

    tempdf3=pd.DataFrame(distances, columns=['x','Distance'])
    tempdf3['sub_group_type']=df1.loc[df1['tomo_set']==prefix,'sub_group_type'].values[0]
    tempdf3['idd']           =df1.loc[df1['tomo_set']==prefix,'idd'].values[0]

    
    df_dist=df_dist.append(tempdf3, ignore_index=True)

    ## loop though the np array 
    lines=[]
    for ind in range(0,len(indices)):
        ind1=indices[ind][0]
        ind2=indices[ind][1]
        ## get x1,y1      x2,y2
        val=[(tempdf1['x_c'].iloc[ind1],tempdf1['y_c'].iloc[ind1]),(tempdf1['x_c'].iloc[ind2],tempdf1['y_c'].iloc[ind2])]
        print(val)
        lines.append(val)


    ################################ plot the line segment and circles ######################
    fig, axes = plt.subplots()
    fig.tight_layout(pad=1)
    fig.set_size_inches(2.5,2.5)
    sns.set_palette(sns.color_palette("Set2"))
    #plt.rcParams["font.weight"] = "bold"
    #plt.rcParams["axes.labelweight"] = "bold"
    #matplotlib.rcParams.update({'font.size': 22})
   
    ### leave the board 
    #plt.gcf().subplots_adjust(left=0.25)
    #plt.gcf().subplots_adjust(bottom=0.35)
    sns.scatterplot(x="x", y="y", data=tempdf2,ax=axes,s=5,linewidth=0,color='#637eb6')
    sns.scatterplot(x="x_c", y="y_c", data=tempdf1,ax=axes,s=10)
    # draw line segment
    lc = mc.LineCollection(lines, linewidths=2)
    axes.add_collection(lc)
    axes.xaxis.set_visible(False)
    axes.yaxis.set_visible(False)
    ## save the fig
    fig.savefig('Contour-DistLine-'+df1.loc[df1['tomo_set']==prefix,'sub_group_type'].values[0]+'-'+str(df1.loc[df1['tomo_set']==prefix,'idd'].values[0])+'.png', dpi=400)

    
    #print(df_dist)



#break


    ###plot all the  distance ###################################################################################################################
fig, axes = plt.subplots(nrows=2, ncols=3)
fig.tight_layout(pad=3.5)
fig.set_size_inches(15,10)
sns.set_palette(sns.color_palette("Set2"))
#plt.rcParams["font.weight"] = "bold"
#plt.rcParams["axes.labelweight"] = "bold"
#matplotlib.rcParams.update({'font.size': 22})


### leave the board 
#plt.gcf().subplots_adjust(left=0.25)
#plt.gcf().subplots_adjust(bottom=0.35)
print(df_dist)
df_dist['Distance']=df_dist['Distance']*1.168  ## unit nm/ a pixel size
print(df_dist)
sns.barplot(x="sub_group_type", y="Distance", data=df_dist, capsize=.2,ax=axes[0,0])

axes[0,0].set_ylabel('Length (nm)',fontsize=14)
axes[0,0].set_xlabel('Condensate Type',fontsize=14)
axes[0,0].set(ylim=(0,200))
axes[0,0].tick_params(labelsize=14)
axes[0,0].set_title( label='Nearest Neighbor Dist.',fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top
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

report_plot_stat(axes[0,0])


fig.savefig('SpC-neighborDist.png', dpi=400)
df_dist.to_csv('SpC-neighborDist.csv', index=False)

plt.show()

