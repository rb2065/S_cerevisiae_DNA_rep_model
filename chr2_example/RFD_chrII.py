# This finds the replication fork directionality (RFD) at each kb position
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
bcsOutputName = '15FF0p05d_fitted_s500_chrII'# 500 simulations of a fitted model for chromosome II with F=15 and recycling rates of 0.05
label = '15FF0p05d_chrII'  # to identify which version of the model is being used
todaysDate = datetime.date.today().strftime("%d%m%Y")

# making a new folder to store the output
new_folder = f"{todaysDate}_RFD_{label}"
folderPath = f'output/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

chrII_length = 813 # length of chromosome II (kb)

# initiating a list where indices represent the positions and values count the number of times that position was replicated by a right fork over all the simulations
rightFork_list = [0.0]* (chrII_length +1)

#------------------ 02.1 ANALYSIS_OUTPUT ------------------

bcsOutputPath = f"bcs_output/{bcsOutputName}.simulation.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()
 
iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations

numberOfIterations = len(iterations)
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyDone = []
    for line in iteration.splitlines(): # iterating over each line in the simulation
        splitLine = line.split('\t')
        if len(splitLine) == 5 and splitLine[2] in ["FR"]:  # lines where DNA is replicated by a rightward moving fork
            repPos= int(splitLine[4])
            if repPos not in alreadyDone: # prevents replication of the same position being accounted for multiple times within the same simulation
                rightFork_list[repPos]+=1 # counting how many times each position is replicated by a rightward moving fork
                alreadyDone.append(repPos)

#------------------ 02.2 FURTHER_ANALYSIS ------------------

# Creating a dataframe storing the RFDs of each position

rightFork_list = [x / numberOfIterations for x in rightFork_list] # finding the fraction of simulations each position is replicated by a rightward moving fork
simRFD_df=pd.DataFrame({'fork_right':rightFork_list})
simRFD_df['fork_right']=(simRFD_df['fork_right']*2)-1 # rescaling to fall in a range between -1 and 1
simRFD_df['i']=simRFD_df.index
simRFD_df = simRFD_df.sort_values(by= 'i').reset_index(drop=True)


# saving a dataframe containing simulated RFD data
simRFD_df.to_csv(f'{folderPath}/{todaysDate}RFD_{label}.csv', index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
