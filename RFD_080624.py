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
bcsOutputName = "loading200FF0p05d_fitted_s500"
label = "loading200FF0p05d"
todaysDate = datetime.date.today().strftime("%d%m%Y")

# making a new folder to store the output
new_folder = f"{todaysDate}_RFD_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

chr_lengths = [230, 813,317,1532,577,270,1100,563,440,746,667,1078,924,784,1091,948] # lengths of chromosomes (kb)
# dictionary mapping each chromosome to a list the length of the chromosome where indices represent the positions and values count the number of times that position was replicated by a right fork over all the simulations
chr2FR = {chromosome_num: [0.0] * (length+1) for chromosome_num, length in zip(range(1, 17), chr_lengths)}

#------------------ 02.1 ANALYSIS_OUTPUT ------------------

bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()
 
iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations

numberOfIterations = len(iterations)
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyDone = []
    for line in iteration.splitlines(): # iterating over each line in the simulation
        splitLine = line.split('\t')
        if len(splitLine) == 9 and splitLine[2] in ["FR"]:  # lines where DNA is replicated by a rightward moving fork
            repPos= int(splitLine[4])
            chromosome = int(splitLine[6])
            chrPos = (chromosome,repPos)
            if chrPos not in alreadyDone: # prevents replication of the same position being accounted for multiple times within the same simulation
                chr2FR[chromosome][repPos]+=1 # counting how many times each position is replicated by a rightward moving fork
                alreadyDone.append(chrPos)

#------------------ 02.2 FURTHER_ANALYSIS ------------------

# dictioonary mapping each chromosome to a dataframe storing the RFDs of each position on that chromosome
chr2simRFD_df = {}
for chromosome in chr2FR.keys():
    chr2FR[chromosome] = [x / numberOfIterations for x in chr2FR[chromosome]] # finding the fraction of simulations each position is replicated by a rightward moving fork
    df=pd.DataFrame({'fork_right':chr2FR[chromosome]})
    df['fork_right']=(df['fork_right']*2)-1 # rescaling to fall in a range between -1 and 1
    df['i']=df.index
    df['ch']=chromosome
    df = df.sort_values(by= 'i').reset_index(drop=True)
    chr2simRFD_df[chromosome]=df

# saving a dataframe containing simulated RFD data
rfd_df = pd.concat(chr2simRFD_df.values(), ignore_index=True)
rfd_df = rfd_df.sort_values(by=['ch', 'i']).reset_index(drop=True)
rfd_df.to_csv(f'{folderPath}/{todaysDate}RFD_{label}.csv', index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
