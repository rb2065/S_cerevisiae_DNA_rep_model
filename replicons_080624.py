# This analysis the output from an alterntive version of the model in which each fork has an additional parameter "ori" which is the position of the origin from which it initiated
# This finds the mean replicon length for each origin
# The results are saved in dataframes (which can later be plotted)
# It does not include running the Beacon Calculus (bcs) script

################################################################################
import pandas as pd
import numpy as np
import os
import math
import datetime
import time
#------------------ 01 SETUP ------------------
start_time = time.time() # for tracking how long the script takes to run
bcsOutputName = "200FF0p05d_fitted_s500" # 500 simulations of a fitted model with F=500 and recycling rates of 0.05
label = '200FF0p05d'  # to identify which version of the bcs model is being used
todaysDate = datetime.date.today().strftime("%d%m%Y")

# getting the positions of origins
inputFilePath = "file/path/to/origin_positions.csv" # contains the position ('i') and chromosome number ('ch') of each origin
ori_df = pd.read_csv(inputFilePath,header=0)
ori_df = ori_df.assign(emptyColumn=0) # adds a column where every value is zero

# making dictionary mapping each chromosome to a dataframe containgin origin positions
grouped_oris = ori_df.groupby('ch')
chr2ori_df = {}
for group, group_data in grouped_oris:
    group_data = group_data.reset_index(drop=True)
    chr2ori_df [group]=group_data

# seting up a dictionary mapping each chromosome to a nested dictionary mapping each origin on that chromosome to a count of the number of times DNA is nreplicated by forks initiating from that origin (over all simulations)
chr2ori2replicons = {} 
# setting the initial vale to zero for each origin
for key in range(1, 17):
    chr2ori2replicons[key]=dict(zip(chr2ori_df[key]['i'], chr2ori_df[key]['emptyColumn']))

#------------------ 04 OUTPUT_ANALYSIS ------------------
bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations (iterations)
numberOfIterations = len(iterations)
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyDone = []
    
    for line in iteration.splitlines(): # iterating through each line of the simulation output
        splitLine = line.split('\t')
        if len(splitLine) == 11 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated
            repPos= int(splitLine[4]) # position being replicated
            chromosome = int(splitLine[6])
            ori = int(splitLine[10]) # identifying the origin which the replication fork initiated from
            chrPos = (chromosome,repPos)
            if chrPos not in alreadyDone: # prevents replication of the same position being accounted for multiple times within the same simulation
                alreadyDone.append(chrPos)
                chr2ori2replicons[chromosome][ori]+=1 # countig the number of positions replicated by forks from that origin             

# setting up a dictionary mapping each chromosome to a datareame containing the mean replicon length of each origin
chr2replicons_df = {}
for chromosome, ori2replicons in chr2ori2replicons.items():
    for ori, replicon in ori2replicons.items():
        chr2ori2replicons[chromosome][ori]= replicon/numberOfIterations # finding the mean replicon length  
    chr2replicons_df[chromosome] = pd.DataFrame(list(chr2ori2replicons[chromosome].items()), columns = ['i','replicons']) # making the dataframe

# combining the replicon dataframes for each chromosome into one dataframe
fullReplicon_df=pd.DataFrame()
for chromosome, df in chr2replicons_df.items():
    df['ch']=chromosome
    fullReplicons_df = pd.concat([fullReplicon_df,df],ignore_index=True)

#------------------ 03.0 SAVING_RESULTS ------------------

# making a new folder to store the output
new_folder = f"{todaysDate}_replicons_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving dataframes
fullReplicon_df.to_csv(f"{folderPath}/fullReplicon_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
