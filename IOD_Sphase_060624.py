
# This calculates the duration of each simulation (simulated length of S-phase) and the inter-origin distances (IODs)
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
bcsOutputName = "200FF0p05d_fitted_s500" # 500 simulations of a fitted model with F=500 and recycling rates of 0.05
label = '200FF0p05d'  # to identify which version of the bcs model is being used

# for storing IODs for each chromosome
# dictionary mapping each chromosome to a list of IODs for that chromosome
chr2iod_list = {}
for chromosome in range(1, 17):
    chr2iod_list[chromosome]=[]

# for finding IODs
# dictionary mapping each chromosome to a nested dictionary mapping each simulation number to a list of positions of origins which fired on that chromosome during that simulation
chr2iteration2firePos_list = {}
for chromosome in range(1, 17):
    chr2iteration2firePos_list[chromosome] = {}

#------------------ 02 ANALYSIS_OUTPUT ------------------
bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

### bcs output analysis 
iterations = bcs_output.split(">")[1:] # splitting it into seperate itterations (simulations)

numberOfIterations = len(iterations) 

for chromosome, iteration2firePos_list in chr2iteration2firePos_list.items():
    for i in range(0,numberOfIterations):
        chr2iteration2firePos_list[chromosome][i]=[] # for each chromosome for each simulation, making a list to store the positions of origins which fire

maxRepTimes = [] # for recording the simulated lengths of S-phase
totalFireCount = [] # for tracking the number of origins which fire i each simualtion 
for index, iteration in enumerate(iterations):  # iterating over each itteration
    repTimesTrack = []
    oriFireCount = 0
    alreadyDone = []
    chr2firePos_list = {}
    for key in range(1, 17):
        chr2firePos_list[key] = []

    for line in iteration.splitlines():
        splitLine = line.split('\t')

        if len(splitLine) == 9 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated
            repPos= int(splitLine[4]) # the position being replicated
            repTime = int(round(float(splitLine[0]))) # time of replication 
            chromosome = int(splitLine[6]) 
            chrPos = (chromosome,repPos) 
            if chrPos not in alreadyDone: # prevents accounting for replication of the same position multiple times within the same simulation
                repTimesTrack.append(repTime) # recording the replication timings (to later calculate the maximum duration of each simulation)
                alreadyDone.append(chrPos)
        
        if len(splitLine) == 11 and splitLine[1] in ["factor"] and splitLine[9] in ["fire"]: # lines where origins fire
            firePos = int(splitLine[4]) # the position of the origin which is firing          
            chromosome = int(splitLine[6])
            chr2iteration2firePos_list[chromosome][index].append(firePos)
            oriFireCount+=1
    totalFireCount.append(oriFireCount)
    maxRepTimes.append(max(repTimesTrack)) # finding the duration of each simulation (which are the simulated lengths of S-phase)

    # finding the IODs
    for chromosome, iteration2firePos_list in chr2iteration2firePos_list.items():
        firePos_list = sorted(iteration2firePos_list[index])
        for i in range(len(firePos_list)-1):
            dist =abs(firePos_list[i+1]-firePos_list[i]) # finding the distance between consecutive origins which fire
            chr2iod_list[chromosome].append(dist)


#------------------ 03 ANALYSIS ------------------

# Finding the mean and std in IODs over the whole genome
iod_list = [item for sublist in chr2iod_list.values() for item in sublist]
meanIOD = np.mean(iod_list).round(decimals=3)
sdIOD = np.std(iod_list).round(decimals=3)

# Finding the mean and std in the number of origins which fire in each simulation
meanFireCount = np.mean(totalFireCount).round(decimals=3)
sdFireCount = np.std(totalFireCount).round(decimals=3)

# Finding the mean and std in the simulation durations (durations of S-phase)
meanRepTime = np.mean(maxRepTimes).round(decimals=3)
sdRepTime = np.std(maxRepTimes).round(decimals=3)

# (these values can be printed if required)

#------------------ 04 SAVING_DATAFRAMES ------------------
new_folder = f"{todaysDate}_IOD_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving IOD dataframe 
iod_data = [(chromosome, iod) for chromosome, iod_l in chr2iod_list.items() for iod in iod_l] # Creating a list of tuples (key, value)
iod_df = pd.DataFrame(iod_data, columns=['ch', 'iod'])
iod_df.to_csv(f"{folderPath}/{todaysDate}iod_df_{label}.csv", index=False)

# saving dtaaframe containing maximum simulation durations (S-phase durations)
sPhase_df = pd.DataFrame(maxRepTimes,columns=['time'])
sPhase_df.to_csv(f"{folderPath}/{todaysDate}sPhase_df_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")

