# this finds the efficiencies of each origin
# The results are saved in dataframes (which can later be plotted)
# It does not include running the Beacon Calculus (bcs) script

##################################################################################
import pandas as pd
import os
import datetime
import time
#------------------ 01 SETUP ------------------
start_time = time.time() # for tracking how long the script takes to run
todaysDate = datetime.date.today().strftime("%d%m%Y")

bcsOutputName = '200FF0p05d_fitted_s500' # 500 simulations of a fitted model with F=500 and recycling rates of 0.05
label = '200FF0p05d'  # to identify which version of the bcs model is being used

# getting the positions of origins
inputFilePath = "file/path/to/origin_positions.csv" # contains the position ('i') and chromosome number ('ch') of each origin
ori_df = pd.read_csv(inputFilePath,header=0)

# making dictionary mapping each chromosome to a dataframe containgin origin positions
grouped_oris = ori_df.groupby('ch')
chr2ori_df = {}
for group, group_data in grouped_oris:
    group_data = group_data.reset_index(drop=True)
    chr2ori_df [group]=group_data

# setting up a dictionary to map chromosome to a nested dictionary which will map each origin to the number of times it fires
chr2ori2efficiency = {}
for chromosome in range(1, 17):
    chr2ori2efficiency[chromosome] = {}

for chromosome, df in chr2ori_df.items():
    for ori in df['i']:
        chr2ori2efficiency[chromosome][ori]=0 # the count for the number of times each origin fires starts off at 0

#------------------ 02 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"{bcsOutputName}.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations
numberOfIterations = len(iterations)
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyFired = []
    for line in iteration.splitlines():  # iterating through each line of the simulation output
        splitLine = line.split('\t')
        if len(splitLine) == 11 and splitLine[1] in ["factor"] and splitLine[9] in ["fire"]: # lines where origins fire
            firePos = int(splitLine[4])
            chromosome = int(splitLine[6])
            chrPos = (chromosome,firePos)
            if chrPos not in alreadyFired: # prevents firing of the same origin being counted twice during the same simulation
                alreadyFired.append(chrPos)
                if firePos in chr2ori2efficiency[chromosome].keys():
                    chr2ori2efficiency[chromosome][firePos]+=1 # counting the number of times that origin fires (over all of the simulations)

#------------------ 03 FURTHER_ANALYSIS ------------------

chr2efficiency_df = {} # dictionary mapping each chromosome to a dataframe storing the efficiencies of origins on that chromosome
for chromosome, ori2simEfficiency in chr2ori2efficiency.items():
    for ori, simEfficiency in ori2simEfficiency.items():
        chr2ori2efficiency[chromosome][ori]=simEfficiency/numberOfIterations # averaging over all the simulations
    efficiency_df = pd.DataFrame({'i': chr2ori2efficiency[chromosome].keys(), 'efficiency': chr2ori2efficiency[chromosome].values()})
    chr2efficiency_df[chromosome] = efficiency_df

#------------------ 04 SAVING_RESULTS ------------------
    
new_folder = f"{todaysDate}_efficiencies_{label}"
folderPath = f'path/to/new/folder/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving the results dataframe
totalOriEfficiency_df = pd.concat(chr2efficiency_df.values(), ignore_index=True)
totalOriEfficiency_df.to_csv(f"{folderPath}/{todaysDate}_efficiency_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")