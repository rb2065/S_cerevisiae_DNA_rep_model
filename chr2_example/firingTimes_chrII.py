# This finds the times at which each origin fires over the different simulations
# The results are saved in dataframes (which can later be plotted)
# It does not include running the Beacon Calculus (bcs) script

##################################################################################
import numpy as np
import pandas as pd
import datetime
import os
import time
#------------------ 01 SETUP ------------------
start_time = time.time() # for tracking how long the script takes to run
todaysDate = datetime.date.today().strftime("%d%m%Y")
bcsOutputName = '15FF0p05d_fitted_s500_chrII'# 500 simulations of a fitted model for chromosome II with F=15 and recycling rates of 0.05
label = '15FF0p05d_chrII'  # to identify which version of the model is being used

# making a new folder to store the output
new_folder = f"{todaysDate}_fireTime_{label}"
folderPath = f'output/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# getting the positions of origins
inputFilePath = "data/origin_positions.csv" # contains the position ('i') and chromosome number ('ch') of each origin
ori_df = pd.read_csv(inputFilePath,header=0)

chrII_ori_df = ori_df[ori_df['ch']==2].reset_index(drop=True) # getting the origins just on chromosome II

# making a dictionary which maps origin positions to a list of times at which that origin fired
ori2times = {}
for ori in chrII_ori_df['i']:
    ori2times[ori]=[]

#------------------ 02.1 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations

for index, iteration in enumerate(iterations):  # iterating over each simulation
    for line in iteration.splitlines(): # iterating through each line of the simulation output
        splitLine = line.split('\t')
        
        if len(splitLine) == 7 and splitLine[1] in ["factor"] and splitLine[5] in ["fire"]: # lines where origins fire
            origin = int(splitLine[4])
            fireTime = int(round(float(splitLine[0])))
            ori2times[origin].append(fireTime) # recording the time at which the origin fired

#------------------ 02.2 FURTHER_ANALYSIS ------------------
            
# finding the mean and std for the times at which each origin fires

summery_data = []
for ori, times in ori2times.items():
    mean_time = np.mean(times)
    sd_time = np.std(times)
    data = (ori, mean_time, sd_time)
    summery_data.append(data)

# saving the dataframe containing summary information on origin firing times  
df = pd.DataFrame(summery_data, columns=['i', 'mean_fireTime', 'sd_fireTime'])
df.to_csv(f'{folderPath}/{todaysDate}firingTimesSummery_{label}.csv', index=False)

# saving a dataframe of the times at which origins on chromosome II fired

data = [(origin, time) for origin, times in ori2times.items() for time in times]# Create a list of tuples from the dictionary
df = pd.DataFrame(data, columns=['i', 'time'])
df.to_csv(f'{folderPath}/{todaysDate}chrII_firingTimes_{label}.csv', index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
