import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import math
from scipy import stats
import glob
import matplotlib.gridspec as gridspec

from scipy import optimize
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,mark_inset)



def figure_adj(axes, icoa, elem):
    axes[icoa].xaxis.set_ticks( [0,30,60,90])

    axes[icoa].set_title( label=elem,fontsize=8, loc='center',  y=0.98, fontweight='bold') ## title on top
    if icob != 0 :
        axes[icoa].yaxis.set_visible(False) # hide both ticks and label
    else:
        y_axis = axes[icoa].axes.get_yaxis()
        #y_axis.set_label_text('foo')
        y_axis.get_label().set_visible(False)
    if icoa != 2 :
        axes[icoa].xaxis.set_visible(False)
    else:

        axes[icoa].set_xlabel('Angle (\N{DEGREE SIGN})',fontsize=11)

def nomalized_distribution(bars, linex, liney, numofBins,axe_in1 ,par_num):
    #axe_in 1-> bar, fitted displot
    #axe_in 2-> normalized, inset

    #### adjustment of the fitted line
    df2=pd.DataFrame(columns=['x','y','adj_factor_basedOn_angle'])   # df store the adjustmnet by cos() function
    df2['x']=linex
    df2['y']=liney
    df2['adj_factor_basedOn_angle']= [math.cos(math.radians(x)) for x in df2['x'].values]
    print(df2)


    def test_func(x, a):
        return a * np.cos(x)
    radians_val= [math.radians(x) for x in df2['x'].values]
    params, params_covariance = optimize.curve_fit(test_func, radians_val, df2['y'],  p0=[0.01]) #p0 initail guess of paramneter a
    # normalize the area under curve to 1
    params[0]=np.pi/180/2  ### 0.00873




    axe_in1.plot(df2['x'], test_func(radians_val, params[0]),label='Fitted function', linewidth=2, linestyle='dashed',color='black')
    #axe_in1.plot(df2['x'], df2['adj_factor_basedOn_angle'], color="green")
    #axe_in1.bar(x="x", height="y",width=180/numofBins/3, data=df2, color="green")
    axe_in1.yaxis.set_ticks( [0,0.005,0.01,0.015])
    #axe_in1.set(xlim=(0,90),ylim=(0,0.01))

    # calculated the differencd, by using df2 - cos() function fitting
    df2['y_adj_factor_basedOn_angle2']= [params[0] * math.cos(math.radians(x)) for x in df2['x'].values] 
    df2['y_adj2'] =  df2['y']-  df2['y_adj_factor_basedOn_angle2']


    poly = np.polyfit(df2['x'],df2['y_adj2'],10)   ###### 8th power polynimial fitting
    df2['y_adj3'] = np.poly1d(poly)(df2['x'])

    axe_in2 = axe_in1.inset_axes([0.52, 0.52, 0.44, 0.44]) #create inset

    #sns.lineplot(x="x", y="y_adj2", data=df2, ax=axe_in2)
    sns.lineplot(x="x", y="y_adj3", color='red', data=df2, ax=axe_in2)
    #axe_in2.yaxis.set_ticks( [0,0.003,0.006,0.009])
    axe_in2.set(xlim=(0,90),ylim=(-0.002,0.002))
    axe_in2.axhline(y=0,color='black', linestyle='dashed')

    y_axis = axe_in2.axes.get_yaxis()
    #y_axis.set_label_text('foo')
    y_axis.get_label().set_visible(False)
    x_axis = axe_in2.axes.get_xaxis()
    #y_axis.set_label_text('foo')
    x_axis.get_label().set_visible(False)
    axe_in2.xaxis.set_ticks([0,45, 90])
    axe_in1.xaxis.set_ticks([0,30,60, 90])

    axe_in2.set_xlabel('', fontsize=2)
    axe_in2.set_ylabel('', fontsize=2)


################################################################################

df1 = pd.read_table('TomoNameList.txt',sep='\s+', header=None)  
print(df1)
df1.columns=['type','idd']

alldf= pd.DataFrame(columns=['model_ID','cluster_ID','disk_angle','dist2sk','disk2surf','vol','type','idd'])

for index, row in df1.iterrows():
    prefix=row['type']+'-'+ str(row['idd'])
    tempdf = pd.read_csv(prefix +'-Suface-angles.csv')  
    tempdf.columns=['model_ID','cluster_ID','disk_angle','dist2sk','disk2surf','vol']
    tempdf['type']=row['type']
    tempdf['idd']=row['idd']
    alldf=alldf.append(tempdf, ignore_index=True)


print(alldf)


plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"

fig1, axes1 = plt.subplots(nrows=1, ncols=5)
fig1.subplots_adjust(hspace=0.5)
fig1.suptitle('Distributions of Angles')
fig1.set_size_inches(19,3)
fig1.subplots_adjust(top=0.92,wspace=0.21, bottom=0.18) 

typelist= alldf.type.drop_duplicates()
print(typelist)
numBins=20
m=0
n=0
for typ in ['SP2','SP3','nucleiA','dropS','dropH1']: ## SP1 too small to distinguish the interior/surfaceNCP
    # select the condensates' inner NCPs which > 11 nm (a NCP diameter) far from the outer most shell
    ddf= alldf.loc[(alldf.type== typ) & (alldf.disk2surf > 110)] 
    #print(ddf)
    print(m)

    #sns.distplot(ddf['disk_angle'], rug=False, hist=True, color ='blue', bins = numBins, ax=axes0[m])
    axes1[m].set_title( label=typ,fontsize=14, loc='center',  y=1.1, fontweight='bold') ## title on top


    ## mirrow the dist to generate a full cos func for curve fitting
    dftemp=ddf['disk_angle'].copy()
    dftemp=dftemp.append(0-dftemp)

    #sns.distplot(dftemp, rug=False, hist=True, bins = numBins, ax=axes0[m])
    fig3 = plt.figure(3)  # templ figure to collect information of bar height
    ax3 = fig3.gca()
    barss=[h.get_height() for h in sns.distplot(dftemp, rug=False, hist=True, bins = numBins, ax=ax3).patches] # get bars height
    linexx,lineyy = sns.distplot(dftemp, rug=False, hist=True, bins = numBins, ax=ax3).get_lines()[0].get_data() # get the fitted line x list and y list
    print(barss)   
    plt.close(3) # close the figure 2, ready for next iteration

    nomalized_distribution(barss, linexx, lineyy, numBins, axes1[m], len(dftemp))

    #sns.histplot(dftemp, ax=axes1[m],bins= numBins,stat="density")


    co1=sns.color_palette("Set2").as_hex()[0]
    co2=sns.color_palette("Set2").as_hex()[1]
    sns.distplot(dftemp, hist=True, rug=False, color=co1, ax=axes1[m],bins= numBins,kde_kws=dict(linewidth=2, color=co2),hist_kws=dict(edgecolor="white", linewidth=2))   ### here check the
    axes1[m].set(xlim=(-0,90),ylim=(0,0.017))

    #figure_adj(axes1, m, n, typ)
    axes1[m].set_ylabel("Probability Density",fontsize=14)
    if m != 0 :

        axes1[m].axes.yaxis.set_ticklabels([])
        axes1[m].set_ylabel(" ")
    axes1[m].set_xlabel("Angle (°)", fontsize=14)


    plt.setp(axes1[m].get_xticklabels(), fontsize=14)
    plt.setp(axes1[m].get_yticklabels(), fontsize=14)

        #axes1[m].yaxis.set_visible(False) # hide both ticks and label
    #break
    m+=1


plt.show()
fig1.savefig('InteriorNCP_Orientation.png', dpi=400)



