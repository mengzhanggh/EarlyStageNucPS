
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
import os
from sklearn.decomposition import PCA


dataframeall= pd.read_csv('NCP_3D_coordinates.csv')
dataframeall['tomo_set']=dataframeall['sub_group_type']+'_'+dataframeall['idd'].astype(str)
unique_vals = dataframeall['tomo_set'].drop_duplicates().tolist()

result_df = pd.DataFrame(columns=['tomo_set','sub_group_type','idd','clustersubgroup','PC1-axis','PC2-axis','PC3-axis','eccentricity', 'numPar'])
print(dataframeall)


targets = [dataframeall.loc[dataframeall['tomo_set'] == val] for val in unique_vals] #list of dataframe
for target in targets:  ## loop through each tomo set

    numNuc= len(target.index)
    print(target.iloc[1]['tomo_set'])

    unique_clusters = sorted(target['label'].drop_duplicates().tolist())

    ## loop through each subgroup condensate cluster
    for unique_cluster in unique_clusters:
        print(unique_cluster)
        ### do PCA in 3D ##########################################
        pca = PCA(n_components=3)
        temp_xyz= target.loc[target['label']==unique_cluster, ['x','y','z']]
        temp_xyz =temp_xyz.reset_index() #  reset index
        temp_xyz =temp_xyz.drop('index', 1) #remove the cluster, which grouped as noise
        numNCP_in_Cluster = len(temp_xyz)
        #print(temp_xyz)
        #        x     y     z
        #0    2840   356   837
        #1    2766  3510   822
        #2    4740  5430   896
        #3    3372  2958  1032
        #4    3414  5923   773


        pca.fit(temp_xyz)
        result=pd.DataFrame(pca.transform(temp_xyz), columns=['PCA%i' % i for i in range(3)], index=temp_xyz.index) 
        #        print(result)
        #            PCA0         PCA1        PCA2
        #0   -1542.705926  1827.311376   -8.622210
        #1     227.485138  -784.025618  -36.184770
        #2    2947.741896 -1201.356085 -156.975715
        #3     397.164369    18.290508 -244.841887

        #print(np.var(result['PCA0']))
        #print(np.var(result['PCA1']))
        #print(np.var(result['PCA2']))

        pc0_dist= abs(result['PCA0'].max() - result['PCA0'].min())
        pc1_dist= abs(result['PCA1'].max() - result['PCA1'].min())
        pc2_dist= abs(result['PCA2'].max() - result['PCA2'].min())
        #determine the min-max distance 
        if pc0_dist> pc2_dist:
            l=pc0_dist
            s=pc2_dist
        else:
            l=pc2_dist
            s=pc0_dist
        #ecc = np.sqrt(np.square(l)-np.square(s))/l
        ecc = 1 - s/l
        print(ecc)
        ##now append the data to the large df
        result_df.loc[len(result_df)] = [target.iloc[1]['tomo_set'],target.iloc[1]['sub_group_type'],target.iloc[1]['idd'],unique_cluster,pc0_dist,pc1_dist,pc2_dist,ecc,numNCP_in_Cluster]

print(result_df)
result_df.to_csv('sp-statistics.csv',index=False)

