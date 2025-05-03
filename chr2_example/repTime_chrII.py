 # This finds the mean and std replication timing at each position
# The results are saved in dataframes (which can later be plotted)
# It does not include running the Beacon Calculus (bcs) script

##################################################################################

import pandas as pd
import numpy as np
import os
import datetime
import time
#------------------ 01 SETUP ------------------
start_time = time.time() # for tracking how long the script takes to run
todaysDate = datetime.date.today().strftime("%d%m%Y")
bcsScriptName = "15FF0p05d_fitted_s500_chrII" # 500 simulations of a fitted model for chromosome II with F=15 and recycling rates of 0.05
label = '15FF0p05d_chrII'  # to identify which version of the model is being used

# setting up a dictionary mapping each position to a list of the times at which that position was replicated (over the different simulations)

pos2repTimes = {}

#------------------ 02 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"bcs_output/{bcsScriptName}.simulation.bcs"

with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting the bcs output into seperate simulations (iterations)

for index, iteration in enumerate(iterations):  # iterating over each simulation

    alreadyDone = []
    for line in iteration.splitlines(): # iterating over each line in the simulation
        splitLine = line.split('\t')
        if len(splitLine) == 5 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated
            repPos= int(splitLine[4]) # position being replicated
            repTime = int(round(float(splitLine[0]))) 
            if repPos not in alreadyDone:
                pos2repTimes[repPos].append(repTime)
            else:
                pos2repTimes[repPos] = repTime
                alreadyDone.append(repPos)

# finding the mean, std, median, and upper and lower qurtiles in replication timing at aech kb

data_list = []
for i, repTimes in pos2repTimes.items():
    mean_repTime = np.mean(repTimes)
    sd_repTime = np.std(repTimes)
    median_repTime = np.percentile(repTimes, 50)
    lq_repTime = np.percentile(repTimes, 25)
    uq_repTime = np.percentile(repTimes, 75)
    row_data = (i, mean_repTime, sd_repTime, lq_repTime, median_repTime, uq_repTime)
    data_list.append(row_data)

repTime_df = pd.DataFrame(data_list, columns=['i', 'repTime', 'repTime_sd', 'repTime_lq', 'repTime_mid', 'repTime_uq'])    
repTime_df.sort_values(by=['i'], inplace=True)

#------------------ 03 SAVING_RESULTS ------------------

# making a new folder to store the output
new_folder = f"{todaysDate}_repTime_std_{label}"
folderPath = f'output/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving the data frame
repTime_df.to_csv(f"{folderPath}/{todaysDate}RepTime_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
