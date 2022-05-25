import numpy as np
import matplotlib.pyplot as plt
from utilities import *
from paircorrelation import pairCorrelationFunction_3D_meng
from paircorrelation import pairCorrelationFunction_3D
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
from kneed import KneeLocator
from sklearn.cluster import AgglomerativeClustering

def deleOut(templist, boundar):
    ## delete ouliars outside boundary
    templist = templist[templist['x'] >= boundar.iloc[0]['x1']]  
    templist = templist[templist['y'] >= boundar.iloc[0]['y1']]  
    templist = templist[templist['z'] >= boundar.iloc[0]['z1']]  
    templist = templist[templist['x'] <= boundar.iloc[0]['x2']]  
    templist = templist[templist['y'] <= boundar.iloc[0]['y2']]  
    templist = templist[templist['z'] <= boundar.iloc[0]['z2']] 
    return templist

def expand(templist, boundar):
    ##copy the data, arrage it into up, down left, right, up-left,down-lef..... 8 dirction
    xshift = boundar.iloc[0]['x2']-boundar.iloc[0]['x1']
    yshift = boundar.iloc[0]['y2']-boundar.iloc[0]['y1']
    templistappend= templist

    temp = templist.copy()
    temp['x']= temp['x']+xshift
    templistappend =  templistappend.append(temp)

    temp = templist.copy()
    temp['x']= temp['x']-xshift
    templistappend =  templistappend.append(temp)

    temp = templist.copy()
    temp['y']= temp['y']+yshift
    templistappend =  templistappend.append(temp)

    temp = templist.copy()
    temp['y']= temp['y']-yshift
    templistappend =  templistappend.append(temp)

 #####
    temp = templist.copy()
    temp['x']= temp['x']+xshift
    temp['y']= temp['y']+yshift
    templistappend =  templistappend.append(temp)

    temp = templist.copy()
    temp['x']= temp['x']-xshift
    temp['y']= temp['y']+yshift
    templistappend =  templistappend.append(temp)

    temp = templist.copy()
    temp['x']= temp['x']+xshift
    temp['y']= temp['y']-yshift
    templistappend =  templistappend.append(temp)


    temp = templist.copy()
    temp['x']= temp['x']-xshift
    temp['y']= temp['y']-yshift
    templistappend =  templistappend.append(temp)

    return templistappend

import glob
#readfilelist = glob.glob('./*parse')

df1= pd.read_csv('NCP_3D_coordinates.csv')
df2= pd.read_csv('Tomo-boundary.csv')

df1.loc[df1['sub_group_type']=='SP1','sub_group_type']='LS'
df1.loc[df1['sub_group_type']=='SP2','sub_group_type']='SP1'
df1.loc[df1['sub_group_type']=='SP3','sub_group_type']='SP2'


print(df1)
print(df2)

dr=16
counter=1
grdfall=pd.DataFrame(columns=['gr','r','sub_group_type','idd'])
for group in ['LS','SP1','SP2','nucleiA','dropS','dropH1']:

    #get the group tomo
    df11=df1.loc[df1['sub_group_type']==group]
    unique_tomos = df11['idd'].drop_duplicates().tolist()
    #get each tomo
    for tomo in unique_tomos:

        templist=df11.loc[df11['idd']==tomo]
        ## append selected data to the df3

        boundary=df2.loc[(df2['sub_group_type']==group)&(df2['idd']==tomo)]


        templist=deleOut(templist, boundary)
        templist=expand(templist, boundary)
        df = templist

        domain_size= min([df['x'].max()-df['x'].min(),df['y'].max()-df['y'].min()]) # calc the min length of x and y
        rMax = domain_size / 4  # define rMax to be measured relative to the tomosize
        #   print(domain_size)
        print(rMax)
        zmin=mean(df['z'])-2*np.std(df['z'])
        zmax=mean(df['z'])+2*np.std(df['z']) 
        #deltaz=df['z'].max()-df['z'].min()
        deltaz=10*np.std(df['z']) # thickness of the reference density
        #print(deltaz)
        #print(np.std(df['z']))
        bound=boundary.iloc[0][['x1','x2','y1','y2','z1','z2']]

        g_r, r, reference_indices = pairCorrelationFunction_3D_meng(df['x'].values, 
                                                                    df['y'].values, 
                                                                    df['z'].values,
                                                                    df['x'].min(),
                                                                    df['x'].max(),
                                                                    df['y'].min(),
                                                                    df['y'].max(),
                                                                    zmin,
                                                                    zmax,
                                                                    deltaz, 
                                                                    rMax, 
                                                                    dr, 
                                                                    bound)



        grdf=pd.DataFrame(columns=['gr','r'])
        grdf['gr']=g_r
        grdf['r']=r
        grdf['sub_group_type']=group
        grdf['idd']=tomo
        grdf['gr']=grdf['gr']/np.max(grdf['gr'].values)
        grdfall=grdfall.append(grdf, ignore_index=True)



grdfall.to_csv('Tomo_Gr.csv', index=False)












