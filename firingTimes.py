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
bcsOutputName = "200FF0p05d_fitted_s500" # 500 simulations of a fitted model with F=500 and recycling rates of 0.05
label = '200FF0p05d'  # to identify which version of the model is being used

# making a new folder to store the output
new_folder = f"{todaysDate}_fireTime_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# getting the positions of origins
inputFilePath = "file/path/to/origin_positions.csv" # contains the position ('i') and chromosome number ('ch') of each origin
ori_df = pd.read_csv(inputFilePath,header=0)

chr2origins = {}
grouped_dataframe = ori_df.groupby('ch')
for group, group_data in grouped_dataframe:
    group_data = np.array(group_data['i'])
    chr2origins[group]=group_data

# making a dictionary which maps each chromosome to a nested dictionary which maps origin positions to a list of times at which that origin fired
chr2ori2times = {}
for chromosome in range(1,17):
    chr2ori2times[chromosome]={}
    for ori in chr2origins[chromosome]:
        chr2ori2times[chromosome][ori]=[]

#------------------ 02.1 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations

for index, iteration in enumerate(iterations):  # iterating over each simulation
    for line in iteration.splitlines(): # iterating through each line of the simulation output
        splitLine = line.split('\t')
        
        if len(splitLine) == 11 and splitLine[1] in ["factor"] and splitLine[9] in ["fire"]: # lines where origins fire
            origin = int(splitLine[4])
            fireTime = int(round(float(splitLine[0])))
            chromosome = int(splitLine[6])
            chr2ori2times[chromosome][origin].append(fireTime) # recording the time at which the origin fired

#------------------ 02.2 FURTHER_ANALYSIS ------------------
            
# finding the mean and std for the times at which each origin fires
summery_data = []
for chromosome, ori2times in chr2ori2times.items():  
    for ori, times in ori2times.items():
        mean_time = np.mean(times)
        sd_time = np.std(times)
        data = (chromosome, ori, mean_time, sd_time)
        summery_data.append(data)
# saving the dataframe containing summary information on origin firing times  
df = pd.DataFrame(summery_data, columns=['ch','i', 'mean_fireTime', 'sd_fireTime'])
df.to_csv(f'{folderPath}/{todaysDate}firingTimesSummery_{label}.csv', index=False)

# for each chromosome, saving a dataframe of the times at which origins on that chromosome fired
for chromosome, ori2times in  chr2ori2times.items():
    data = [(origin, time) for origin, times in ori2times.items() for time in times]# Create a list of tuples from the dictionary
    df = pd.DataFrame(data, columns=['i', 'time'])
    df.to_csv(f'{folderPath}/{todaysDate}chr{chromosome}firingTimes_{label}.csv', index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
