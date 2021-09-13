import numpy as np
import os
from natsort import natsorted
import sys
#import pandas as pd
#import matplotlib.pyplot as plt
#import seaborn as sns
#import re
#import subprocess

if len(sys.argv) != 4:
        print(f"Three argument are needed, run {sys.argv[0]} system_name output_folder pdb_folder\n");
        exit()

system_name = sys.argv[1]
output_folder = sys.argv[2]
pdb_folder = sys.argv[3]


os.system(f"mkdir {output_folder}")
os.chdir(f"./{output_folder}")
#os.system("touch test")
#os.system("ls")


files = np.array([])
filename_vol_all= ""
a_directory = f"/media/yulian/MyBook/yulian/projects/gr/md_out/{system_name}/{pdb_folder}/"
for filename in os.listdir(a_directory):

    filepath = os.path.join(a_directory, filename)

    if filename[-3:]=="pdb":

        os.system(f"voronota-volumes -i {filepath} --per-residue  --sum-at-end > {filename[:-3]}vol")
        #os.system(f"voronota-volumes -i {filepath} --per-residue  --sum-at-end > {filepath[:-3]}vol")

        filename_vol = f"{filename[:-3]}vol"
        #filename_vol = f"{filepath[:-3]}vol"

        filename_vol_all = filename_vol_all+' '+filename_vol

filename_vol_all_list = filename_vol_all.split(" ")
filename_vol_all_list = filename_vol_all_list[1:]

filename_vol_all_list_sorted = natsorted(filename_vol_all_list)

#print(sorted_filename_vol_all_list_sorted)

os.system("rm all_part1.vol all_part2.vol all.vol ")

j=0
for i in filename_vol_all_list_sorted:

    #!tail -1 {i} >> all.vol
    os.system(f"tail -1 {i}| cut -c 18-40 >> all_part2.vol")
    os.system(f"echo {j} >> all_part1.vol")
    j=j+1
    os.system(f"paste -d\" \" all_part1.vol all_part2.vol > {system_name}_all.vol")



