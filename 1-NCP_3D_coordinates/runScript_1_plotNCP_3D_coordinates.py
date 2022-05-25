import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats


df= pd.read_csv('NCP_3D_coordinates.csv')
df['tomo_set']=df['sub_group_type']+'_'+df['idd'].astype(str)
unique_group = df['tomo_set'].drop_duplicates().tolist()
print(unique_group)

targets = [df.loc[df['tomo_set'] == val] for val in unique_group] 
for target in targets:
    templist=target
    fig = plt.figure(figsize=(3,3))


    ax = fig.add_subplot(projection='3d')
    ax.set_title( label=target.iloc[1]['sub_group_type']+'-'+str(target.iloc[1]['idd']),fontsize=10, loc='left', x=0.6, fontweight='bold') ## title on top

    ## flip data around axis for display only
    if target.iloc[1]['sub_group_type']=='SP3' and target.iloc[1]['idd']==3 :
        templist['x']= 512*11.68 - templist['x']   


    ax.scatter(templist['x'], templist['y'], templist['z'], c= templist['color'], marker='o',s=4)

    ###loop though each label, and add text
    for uniqulabel in sorted(templist['label'].drop_duplicates().tolist())[:-1]:
        xmean= templist.loc[templist['label']== uniqulabel, 'x' ].mean()   
        ymean= templist.loc[templist['label']== uniqulabel, 'y' ].mean()   
        zmean= templist.loc[templist['label']== uniqulabel, 'z' ].mean()   
        ax.text(xmean, ymean, zmean, str(uniqulabel), color='black', size=5)




    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.set_xlim3d(0,6000)
    ax.set_ylim3d(0,6000)
    ax.set_zlim3d(0,3000)
    ax.set_xticks([]) 
    ax.set_yticks([]) 
    ax.set_zticks([])

    ax.xaxis.pane.set_edgecolor('#f0f0f0')  ### light gray
    ax.yaxis.pane.set_edgecolor('#f0f0f0')
    ax.zaxis.pane.set_edgecolor('#f0f0f0')
    ax.xaxis.pane.set_alpha(0.4)
    ax.yaxis.pane.set_alpha(0.4)
    ax.zaxis.pane.set_alpha(0.4)  # plane transparency


    plt.legend(numpoints=1 , loc='upper left')

    ax.view_init(70, -60)   ### (  rotate xy plane up/down,    rotatez)

    #plt.show()
    #plt.savefig(target.iloc[1]['sub_group_type']+str(target.iloc[1]['idd'])+target.iloc[1]['tomo_set'].replace("./","-")+'-DBSCAN-cluster'+str(iterx).zfill(2)+'.png',dpi=200)
    plt.savefig(target.iloc[1]['sub_group_type']+str(target.iloc[1]['idd'])+target.iloc[1]['tomo_set'].replace("./","-")+'-DBSCAN-cluster.png',dpi=200)
    #np.savetxt(target.iloc[1]['tomo_set']+'-DBSCAN-cluster.txt',templist[['x','y','z','label']] , fmt='%d') 


    ## rotated around y for displayback
    if target.iloc[1]['sub_group_type']=='SP3' and target.iloc[1]['idd']==3 :
        templist['x']= 512*11.68 - templist['x']   


    templist['tomo_set']=target.iloc[0]['tomo_set']
    templist['sub_group_type']=target.iloc[0]['sub_group_type']
    templist['idd']=target.iloc[0]['idd']
    print('=========================================================')


  ###### save data
    
    ##### now convert to cmm files
    ## get unique group id
    unique_g = sorted(templist['label'].drop_duplicates().tolist())
    print(unique_g)
    subgroups = [templist.loc[templist['label'] == xx] for xx in unique_g] ## get each sub dataframe with unique cluster id

    ### open a file to over wirte
    f = open(target.iloc[1]['tomo_set']+'-DBSCAN-cluster.cmm', "w")
    f.write("<marker_sets>\n")
    

    for subgroup in subgroups:
        print(subgroup)
        f.write('<marker_set name="marker set '+str(subgroup.iloc[1]['label'])+'">\n') ### the marker set cluster id

        ## loop through each point in subgroup
        coun=1
        for index, row in subgroup.iterrows():
            ### marker set cluster id and color hexcode

            h = row['color'].lstrip('#')
            colorRGB = tuple(int(h[i:i+2], 16) for i in (0, 2, 4))
            
            f.write('<marker id="'+str(coun)+'" x="'+str(row['x'])+'" y="'+str(row['y'])+'" z="'+str(row['z'])+'" r="'+str(colorRGB[0]/255)+'" g="'+str(colorRGB[1]/255)+'" b="'+str(colorRGB[2]/255)+'" radius="40"/>\n') 
            coun+=1
        

        
        f.write("</marker_set>\n")
    # finished writteing 
    f.write("</marker_sets>")
    f.close()


    #break


