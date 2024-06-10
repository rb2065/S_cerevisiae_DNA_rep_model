# This finds the mean number of firing factors available at each time point.
# The results are saved in dataframes (which can later be plotted)
# It does not include running the bcs script

##################################################################################
import pandas as pd
import numpy as np
import math
import datetime
import os
import time
#------------------ 01 SETUP ------------------
start_time = time.time() # for tracking how long the script takes to run
todaysDate = datetime.date.today().strftime("%d%m%Y")
bcsOutputName = "200FF0p05d_fitted_s500" # 500 simulations of a fitted model with F=500 and recycling rates of 0.05
label = '200FF0p05d'  # to identify which version of the model is being used
maxFF = 200 # maximum number of firing factors included in the bcs model

#------------------ 02 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"{bcsOutputName}.bcs"

with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

### bcs output analysis 
iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations (iterations)

time2numberOfFF = {} # dictionary mapping time to the number of available firing factors at that time
maxTimes = [] # for capping the time plotted
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyDone = []
    numberOfFF = maxFF # tracking number of free firing factors (they are all available at that start of each simulation)
    repTimesTrack = [] # tracking replication timings to find the max replication times for each simulation
    if 0.0 not in time2numberOfFF:
        time2numberOfFF[0.0] = [numberOfFF]
    else:
        time2numberOfFF[0.0].append(numberOfFF)

    for line in iteration.splitlines(): # iterating through each line of the simulation output
        splitLine = line.split('\t')

        if len(splitLine) == 9 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated
            repTime = math.ceil(float(splitLine[0])*2)/2 # rounded to the 0.5 above
            repPos= int(splitLine[4])
            chromosome = int(splitLine[6])
            chrPos = (chromosome,repPos)
            if chrPos not in alreadyDone:
                repTimesTrack.append(repTime)
                alreadyDone.append(chrPos)

        if len(splitLine) == 3 and splitLine[1] in ["factor"]:  # lines where firing factors are used (when origins fire)
            repTime = math.ceil(float(splitLine[0])*2)/2 # rounded to the 0.5 above
            numberOfFF-=1 # decreasing the number of available firing factors when an origin fires
            if repTime not in time2numberOfFF:
                time2numberOfFF[repTime] = [numberOfFF]
            else:
                time2numberOfFF[repTime].append(numberOfFF)
            
        elif len(splitLine) == 3 and splitLine[1] in ["dwell"]:  # lines where firing factors are recycled 
            repTime = math.ceil(float(splitLine[0])*2)/2
            numberOfFF +=1 # increasing the number of available firing factors when they recycle
            if repTime not in time2numberOfFF:
                time2numberOfFF[repTime] = [numberOfFF]                
            else:
                time2numberOfFF[repTime].append(numberOfFF)


    maxTimes.append(max(repTimesTrack)) # tracking the duration if S-phase for each simulation (the latest time a position is replicated)

#------------------ 03 FURTHER_ANALYSIS ------------------

# cutting of the dictionaries to not include data for values above the mean maximom time (to miss out the noise at the end)
meanMaxTime = np.round(np.mean(maxTimes), 2)
time2numberOfFF = {key: value for key, value in time2numberOfFF.items() if key <= meanMaxTime}

# finding the mean number of avalable firing factors at each time point over all of the simulations
for repTime, numberOfFF in time2numberOfFF.items():
    time2numberOfFF[repTime]=np.mean(numberOfFF)

# storing the number of available firing factors at each time in a dataframe
numberOfFF_df = pd.DataFrame(list(time2numberOfFF.items()), columns=['repTime', 'numberOfFF'])
numberOfFF_df.sort_values(by='repTime', inplace=True)

#------------------ 04 SAVING_RESULTS ------------------

# making a new folder to store the output
new_folder = f"{todaysDate}_freeFF_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving the dataframe
numberOfFF_df.to_csv(f"{folderPath}/numberOfFF_df_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
