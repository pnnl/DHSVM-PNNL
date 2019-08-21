# --------------------------------------------------------------------------------------------------------------
# This is the Python3 script used to post-process RBM output (<ProjectName>.temp file)
# before running this script, just modify the hard-coded value for your output
# Run the script in Linux:
#       >>> python3 postprocess_rbm_output.py
# Note: the output segment ID is the DHSVM ID, not RBM ID!
#       the output format is consistent with DHSVM Outflow.Only file
# --------------------------------------------------------------------------------------------------------------

import numpy as np
import os
import copy
import math
import re
import csv
from shutil import copyfile
import calendar
import datetime
import shutil


# --------------------------------------------------------------------------------------------------------------
def npnan(x,y):
	# this function creates the np.nan 2d-array (np.nan should be float)
 	array_2d = np.zeros((x,y), float) 
 	array_2d[:] = np.nan
 	return array_2d




# --------------------------------------------------------------------------------------------------------------
# 1. hard-coded values here
# only modify the script here
met_resolution = 3  # the temporal resolution for RBM and DHSVM simulation, in hours

# RBM starting time: dhsvm + <met_resolution> hours (date the 1st line in .Only files + 3 hours)
# RBM ending time:   dhsvm + <met_resolution> hours
s_year  = 2001;    e_year  = 2014
s_month = 10;      e_month = 1
s_day   = 2;       e_day   = 1
s_hour  = 3;       e_hour  = 0

s_date = datetime.datetime(s_year, s_month, s_day, s_hour, 0, 0, 0)
e_date = datetime.datetime(e_year, e_month, e_day, e_hour, 0, 0, 0)

# total steps
num_step = int((e_date - s_date).total_seconds()/60/60/met_resolution + 1)

# number of DHSVM segment
nseg = 797

# RBM output directory
temp_dir = './noaa.temp'

# RBM<->DHSVM segment ID conversion file
seg_dir = './noaa.segmap'

# outflow.only dir
outflowonly_dir = '../dhsvm/output_for_RBM/Outflow.Only'

# output dir
outputfile_dir = './Tw_output'









# --------------------------------------------------------------------------------------------------------------
# 2. load the RBM output file .temp
# 4th column: RBM segment ID
# 6th column: simulated temperature

lines = [line.rstrip('\n') for line in open(temp_dir)]	

# prepare a np.arrary to catch the RBM output
data_rbm = npnan(len(lines), 2)   # RMB seg, water T

# load the data into data_rbm
count = 0
for line in lines:
	item = line.split()  
	data_rbm[count, 0] = float(item[3])
	data_rbm[count, 1] = float(item[5])
	count += 1





# --------------------------------------------------------------------------------------------------------------
# 3. load the RBM<->DHSVM segment conversion
lines = [line.rstrip('\n') for line in open(seg_dir)]	
lines = lines[1:]  # remove the 1st line

# segment link file
seg_link = npnan(nseg, 2)   # 1st column: RBM segment ID; 2nd column DHSVM segment ID

# load the data into seg_link
count = 0
for line in lines:
	item = line.split()  
	seg_link[count, 0] = float(item[1])
	seg_link[count, 1] = float(item[3])
	count += 1




# --------------------------------------------------------------------------------------------------------------
# 4. output data for each DHSVM segment

# prepare a np.array to catch the output
output_dhsvm = npnan(num_step, nseg)   # record temperature for each dhsvm segment

# ordered dhsvm segment id
ascend_dhsvm_id = np.sort(seg_link[:,1])

# loop through each segment
for i in range(nseg):

	# rbm segment id
	rbm_id = seg_link[i,0]

	# find the assoicated dhsvm segment ID
	dhsvm_seg_id = seg_link[seg_link[:,0]==rbm_id, 1]

	# find the dhsvm segment id index in the array
	index_in_array = np.where(ascend_dhsvm_id == dhsvm_seg_id)
	index_in_array = index_in_array[0][0]

	# record data only for the ith RBM segment
	temp_data = npnan(num_step,1)
	temp_data[:,0] = data_rbm[data_rbm[:,0]==rbm_id, 1]  
	output_dhsvm[:,index_in_array] = temp_data[:,0]

# release memory
del data_rbm




# ----------------------------------------------------------------------------------------------------------------------------------------
# 5. prepare the output format, consistent with the Outflow.Only file
output = []
output.append([])
output[0].append('%02d/%02d/%04d-%02d:00:00 %02d/%02d/%04d-%02d:00:00 %d' %(s_month, s_day, s_year, s_hour, e_month, e_day, e_year, e_hour, met_resolution))
output.append([])
lines = [line.rstrip('\n') for line in open(outflowonly_dir)]	
output[1].append(lines[1])

for i in range(num_step):
	output.append([])

	temp_year = float((s_date + datetime.timedelta(seconds=i*3600*met_resolution)).year)
	temp_month = float((s_date + datetime.timedelta(seconds=i*3600*met_resolution)).month)
	temp_day = float((s_date + datetime.timedelta(seconds=i*3600*met_resolution)).day)
	temp_hour = float((s_date + datetime.timedelta(seconds=i*3600*met_resolution)).hour)

	output[i+2].append('%02d.%02d.%04d-%02d:00:00' %(temp_month, temp_day, temp_year, temp_hour))

	for j in range(nseg):
		output[i+2].append('%4.2f' %(output_dhsvm[i,j]))


# write out the final data
myfile = open(outputfile_dir,'w')
wr = csv.writer(myfile, delimiter=' ')
wr.writerows([' '.join(c[:])] for c in output)
myfile.close()

# remove the "" in each line
lines = [line.rstrip('\n') for line in open(outputfile_dir)]	
new_lines = copy.deepcopy(lines)
count = 0
for line in lines:
	word = []                    #a list to contain the string item
	for item in line.split():    #split string into item
		word.append(item)
	word[0] = word[0][1:]
	word[-1] = word[-1][:-1]
	new_lines[count] = ' '.join(word)
	count += 1

with open(outputfile_dir, mode='wt', encoding='utf-8') as myfile:
	myfile.write('\n'.join(new_lines))
myfile.close()	


