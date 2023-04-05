import numpy as np
from tqdm import tqdm
import os.path

file_name = "data/230403/230403_d3_"

job_arr_start = 0
job_arr_end = 9
job_arr_step = 1
job_arr_l = np.arange(job_arr_start,job_arr_end+1,job_arr_step)

cnt = 0

for job_id in tqdm(job_arr_l):
    job_name = file_name + str(job_id)
    if os.path.isfile(job_name+".csv") == False:
        print(job_id)
        continue

    if cnt == 0:
        data = np.loadtxt(job_name+".csv",skiprows=1,delimiter=",")
        x_l = data[:,0]
        y_l = data[:,1]

    else:
        data = np.loadtxt(job_name+".csv",skiprows=1,delimiter=",")
        x_l = x_l + data[:,0]
        y_l = y_l + data[:,1]

    cnt = cnt + 1

x_l = x_l / cnt
y_l = y_l / cnt

np.savez(file_name+"pp.npz",x_l=x_l,y_l=y_l)