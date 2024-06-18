# this analysis the output from an alterntive version of the model in which each fork has an additional parameter "ori" which is the position of the origin from which it initiated
# (this allows individual replication forks to be identified)
# this finds the mean number of replication forks active at each 30 second time points over the simulations
# The results are saved in dataframes (which can later be plotted)
# It does not include running the Beacon Calculus (bcs) script

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
#------------------ 02 OUTPUT_ANALYSIS ------------------

sim2data = {} # dictionary maping the simulation number to data relating to that simulation

bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

### bcs output analysis 
iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations (itterations)

for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyDone = []
    fork2times = {} # dictionary mapping each individual replication fork to a list of times at which it replicated DNA

    for line in iteration.splitlines(): # iterating over each line in the simulation
        splitLine = line.split('\t')
        if len(splitLine) == 11 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated 
            repPos= int(splitLine[4]) # position being replicated
            repTime = int(round(float(splitLine[0]))) # time of replication
            chromosome = int(splitLine[6])
            ori = int(splitLine[10]) # the origin the replication fork originated from
            chrPos = (chromosome,repPos)
            fork = (splitLine[2],ori,chromosome) # identifying the individual replication fork
            if chrPos not in alreadyDone: # prevents replication of the same position being accounted for multiple times within the same simulation
                alreadyDone.append(chrPos)
                if fork in fork2times:
                    fork2times[fork].append(repTime) # tracking the times at which each individual replication fork moves
                else:
                    fork2times[fork]=[repTime]

    # finding the initiation and termination times of each fork (as the earliest and latest times of fork movement)
    data_list = []
    for fork, times_list in fork2times.items():
        initiate_time = min(times_list)
        termination_time = max(times_list)
        row_data1 = ("start",initiate_time)
        row_data2 = ("end",termination_time)
        # making a list of labelled initiation and termination times of forks
        data_list.append(row_data1)
        data_list.append(row_data2)
    # making a dataframe containing fork initiation and termination times for each simulation
    df = pd.DataFrame(data_list, columns=['status', 'time'])
    df = df.sort_values(by=['time']).reset_index(drop=True)
    sim2data[index]=df

#------------------ 03 FURTHER_OUTPUT ------------------

time2activeForks = {} # dictionary maping each time point to a list of the numbers of active replication forks at that time from each simulation
for sim, df in  sim2data.items():
    active_forks = 0
    # increasing and decreasing the count of active forks for times when they initiate and terminate respectively
    for row_tuple in df.itertuples():
        repTime = row_tuple[2]
        if row_tuple[1] == "start":
            active_forks +=1
        elif row_tuple[1] == "end":
            active_forks -=1
        if repTime in time2activeForks:
            time2activeForks[repTime].append(active_forks)
        else:
            time2activeForks[repTime]=[active_forks]

data_list = []    
for repTime, forks_list in time2activeForks.items():
    active_forks = np.mean(forks_list) # finding the mean number of active replication forks at each time
    row_data = (repTime, active_forks)
    data_list.append(row_data)

activeForks_df = pd.DataFrame(data_list, columns=['time', 'activeForks'])
activeForks_df = activeForks_df.sort_values(by=['time']).reset_index(drop=True)

#------------------ 03 SAVING_RESULTS ------------------

# making a new folder to store the output
new_folder = f"{todaysDate}_activeForks_{label}"
folderPath = f'ath/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)
# saving the data frame
activeForks_df.to_csv(f"{folderPath}/{todaysDate}activeForks_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")


