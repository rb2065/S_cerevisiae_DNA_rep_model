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

bcsOutputName = '15FF0p05d_fitted_s500_chrII' # 500 simulations of a fitted model for chromosome II with F=15 and recycling rates of 0.05
label = '15FF0p05d'  # to identify which version of the bcs model is being used

# getting the positions of origins
inputFilePath = "data/origin_positions.csv" # contains the position ('i') and chromosome number ('ch') of each origin
ori_df = pd.read_csv(inputFilePath,header=0)
chrII_ori_df = ori_df[ori_df['ch']==2].reset_index(drop=True)

# setting up a dictionary to map each origin to the number of times it fires

ori2efficiency = {}
for ori in chrII_ori_df['i']:
    ori2efficiency[ori]=0 # the count for the number of times each origin fires starts off at 0

#------------------ 02 OUTPUT_ANALYSIS ------------------

bcsOutputPath = f"bcs_output/{bcsOutputName}.simulation.bcs"
with open(bcsOutputPath, 'r') as file:
    bcs_output = file.read()

iterations = bcs_output.split(">")[1:]   # splitting it into seperate simulations
numberOfIterations = len(iterations)
for index, iteration in enumerate(iterations):  # iterating over each simulation
    alreadyFired = []
    for line in iteration.splitlines():  # iterating through each line of the simulation output
        splitLine = line.split('\t')
        if len(splitLine) == 7 and splitLine[1] in ["factor"] and splitLine[5] in ["fire"]: # lines where origins fire
            firePos = int(splitLine[4])
            if firePos not in alreadyFired: # prevents firing of the same origin being counted twice during the same simulation
                alreadyFired.append(firePos)
                if firePos in ori2efficiency.keys():
                    ori2efficiency[firePos]+=1 # counting the number of times that origin fires (over all of the simulations)

#------------------ 03 FURTHER_ANALYSIS ------------------

# Creating dataframe storing the efficiencies of origins

for ori, simEfficiency in ori2efficiency.items():
    ori2efficiency[ori]=simEfficiency/numberOfIterations # averaging over all the simulations
efficiency_df = pd.DataFrame({'i': ori2efficiency.keys(), 'efficiency': ori2efficiency.values()})

#------------------ 04 SAVING_RESULTS ------------------
    
new_folder = f"{todaysDate}_efficiencies_{label}"
folderPath = f'output/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# saving the results dataframe
efficiency_df.to_csv(f"{folderPath}/{todaysDate}_efficiency_{label}.csv", index=False)

# calculating how long the script took to run
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")

print("Finished")
