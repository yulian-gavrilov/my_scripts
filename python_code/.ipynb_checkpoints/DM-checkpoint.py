import numpy as np
import pandas as pd
import math as m
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm

def make_matrices(coord,total_CA):
    nframes = int(len(coord)/total_CA)
    #print(nframes)
    dm = np.empty(shape=(total_CA,total_CA))
    dm1 = np.empty(shape=(total_CA,total_CA))
    dm_sum = np.empty(shape=(total_CA,total_CA))
    dm_avg = np.empty(shape=(total_CA,total_CA))
    for_std = np.empty(shape=(total_CA,total_CA))
    sum_for_std = np.empty(shape=(total_CA,total_CA))
    std = np.empty(shape=(total_CA,total_CA))
    
    
    for frame in range(0,nframes,1):
        for i in range(0,total_CA,1):
            CAi = total_CA*frame + i
            for j in range(0,total_CA,1):
                CAj = total_CA*frame + j
                dm[i,j] = np.sqrt((coord[0].loc[CAi]-coord[0].loc[CAj])**2 + (coord[1].loc[CAi]-coord[1].loc[CAj])**2 + (coord[2].loc[CAi]-coord[2].loc[CAj])**2)   
    
        dm_sum = dm_sum + dm
    
    dm_avg = dm_sum / nframes
    
    
    for frame in range(0,nframes,1):
        for i in range(0,total_CA,1):
            CAi = total_CA*frame + i
            for j in range(0,total_CA,1):
                CAj = total_CA*frame + j
                dm1[i,j] = np.sqrt((coord[0].loc[CAi]-coord[0].loc[CAj])**2 + (coord[1].loc[CAi]-coord[1].loc[CAj])**2 + (coord[2].loc[CAi]-coord[2].loc[CAj])**2)   
                for_std = (dm1-dm_avg)**2
        
        sum_for_std = sum_for_std + for_std
        
    std = np.sqrt(sum_for_std/nframes)
    
    return dm_avg, std



total_CA = 58
block=total_CA*10

coord = pd.read_csv("./ca_trj/md_nopbc_fit_dt10000.txt",sep='\s+',header=None)
std_all = np.zeros(shape=(total_CA,total_CA,1))
#co = np.zeros(shape=(3,total_CA))


for i in range(0,31,10):
    #print(total_CA*block, total_CA*block+total_CA)
    df_new = coord.iloc[total_CA*i:total_CA*i+block].reset_index(drop=True)
    #print(df_new)
    
    dm_avg, std = make_matrices(df_new,total_CA)
    
    print (std,"\n#####")
#     if block == 0:
#         std_all = std
#         std_all = np.concatenate((std_all[None], [std]))
#     else:
#     std_all = np.dstack((std_all, std))
    
# std_all.shape
