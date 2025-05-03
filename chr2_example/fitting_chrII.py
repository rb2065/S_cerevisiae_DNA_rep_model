# This fits the firing rates of the origins in the Beacon Calculus (bcs) script to experimental replication timing data (Muller et al., 2014)
# The bcs script from the final fitting iteration and the bcs script with the lowest error are saved
# The simulated replication timings and errors over the fitting are saved in  dataframes

##################################################################################

import pandas as pd
import numpy as np
import ast
import os
import datetime

#------------------ 01 SETUP_GENERAL ------------------

todaysDate = datetime.date.today().strftime("%d%m%Y")
number_of_runs = 15 # number of fitting iterations
simulations_per_run = 500 # number of model simulations per iteration of fitting
bcsScriptName = '15FF0p05d_unfitted_chrII.bc'  # to identify which version of the model is being used (eg: chromosome 2 model with 15 firing factors and recycling rates of 0.05)

runs = 0 # tracking the fitting iteration number

# tracking the errors over the fitting iterations
meanAbsoluteError_track = []
meanSquaredError_track = []

chrII_length = 813 # length of chromosome II (kb)

# making output folder
new_folder = f"{todaysDate}_fitting_{bcsScriptName}"
folderPath = f'output/{new_folder}'
os.makedirs(folderPath, exist_ok=True)

# making a reset script (to save a copy of the original bcs script)
bcsScriptPath = f'chr2_example/{bcsScriptName}.bc'
with open(bcsScriptPath, 'r') as file:
    original_bcs_code = file.read()

bcsScriptPath_reset = f'{folderPath}/{bcsScriptName}_reset.bc'      
with open(bcsScriptPath_reset, 'w') as file:
    file.write(original_bcs_code)

#------------------ 02 SETUP_DATAFRAMES ------------------

# making a dictionary which maps iteration run number to a timeStore list

run2timeStore = {}
for key in range(number_of_runs):
    run2timeStore[key] = [0.0 for x in range(0,chrII_length+1)]


# importing data frames
expRepTime_df = pd.read_csv("data/completeExpRepTime.csv", header=0) # contains experimental replication times for all positions
ori_positions_df = pd.read_csv("data/origin_positions.csv", header=0) # contains the positions of origins

# getting experimental replication times just at origins
expOriRepTime_df = ori_positions_df.merge(expRepTime_df[['ch', 'i', 'expRepTime']], on=['ch', 'i'], how='left')

# making dictioary mapping chromosome to experimental replication time data frame

# Getting the experimental replication timings for chromosome 2
chrIIrepTime_df = expRepTime_df[expRepTime_df['ch']==2].reset_index(drop=True)

# making dictioary mapping chromosome to experimental origin replication time data frame

# Getting the experimental replication timings for just chromosome 2 origins
chrIIoriRepTime_df = expOriRepTime_df[expOriRepTime_df['ch']==2].reset_index(drop=True)

#------------------ 03 RUNNING_BCS------------------

while runs < number_of_runs:
    print(f"itteration {runs}")
    
    bcs_command = f"bcs/bin/bcs -s {str(simulations_per_run)} -o bcs_output/{bcsScriptName} " + bcsScriptPath
    os.system(bcs_command)

#------------------ 04 FINDING_REPLICATION_TIMING ------------------

    # Opening the bcs output
    bcsOutputPath = f"bcs_output/{bcsScriptName}.simulation.bcs"
    f = open(bcsOutputPath,'r')

    numberOfSimulations = 0 # for counting the number of simulations in the bcs output

    # iterating through each line of the bcs output
    for line in f:
        if line[0] == '>': # marks the start of a new simulation
            alreadyDone = []
            numberOfSimulations += 1
            continue

        splitLine = line.split('\t')

        if len(splitLine) == 5 and splitLine[2] in ["FL","FR"]:  # lines where DNA is replicated
            position = int( splitLine[4] )
            time = float(splitLine[0] )
            if position not in alreadyDone: 
                run2timeStore[runs][position]+=time
                alreadyDone.append(position) # prevents replication of the same position being counted twice in the same simulation
            
    f.close()   

    # finding the average replication time of each position over all of the simulations
    for i, j in enumerate(run2timeStore[runs]):
        run2timeStore[runs][i] = float(j) / float(numberOfSimulations)

#------------------ 05 EXTRACTING_ORI_PARAMETERS ------------------


    with open(bcsScriptPath, 'r') as file:
        bcs_code = file.read()

    # defining the parts of the Beacon Calculus script which the origin processes parameters start and finish
    start_marker = "initial processes in the system\nORI["
    end_marker = "] || FF[]"

    # finding the idicies in the bcs script where the origin processes parameters start and finish
    start_index = bcs_code.find(start_marker)
    end_index = bcs_code.find(end_marker, start_index)

    # extracting the origin parameters from the script
    if start_index != -1 and end_index != -1:
        ori_parameters = bcs_code[start_index+len(start_marker)-4: end_index +1 ]

    # converting the origin parameters to a list
    parameters_list = ori_parameters.split(" || ")
    # removing the string "ORI" from the parameters
    parameter_data = [ast.literal_eval(lst_str.replace("ORI", "")) for lst_str in parameters_list]

    # making a data frame containing the origin paramiters from the Beacon Calculus script
    parameter_df = pd.DataFrame(parameter_data, columns=["i","fire"]).astype(float)
    parameter_df['i'] = parameter_df['i'].round().astype(int)


    # adding simulated origin replication timings to chr2oriRepTime_df 
    timeStore = run2timeStore[runs]
    chrIIoriRepTime_df['repTime_sim'] = [timeStore[position] for position in parameter_df['i'] ]


#------------------ 06 CHANGING_FIRING_RATES ------------------

    for index, fire in enumerate(parameter_df['fire']):
        new_fire = fire*(chrIIoriRepTime_df.loc[index,'repTime_sim']/chrIIoriRepTime_df.loc[index,'repTime_exp'])**1.2
        new_fire=round(new_fire,5)
        # safeguarding against firing rates <= 0
        if new_fire <= 0:
            new_fire = 0.0001
        parameter_df.loc[index, 'fire'] =new_fire

#------------------ 07 UPDATING_BCS_SCRIPT ------------------

    # putting the parameters with the new firing rates back into the Beacon Calculus script
    modified_parameters = " || ".join(f"ORI[{','.join(map(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() else f'{x:.5f}', row))}]" for row in parameter_df.values.tolist())

    modified_bcs_code = bcs_code[:start_index+len(start_marker)-4] + modified_parameters + bcs_code[end_index+1:]

    # Save the modified Beacon Calculus script
    with open(bcsScriptPath, 'w') as file:
        file.write(modified_bcs_code)

#------------------ 08 ERROR ------------------
    
    # combining the simulated replication times from all chromosomes into one list

    absolute_errors = abs(expRepTime_df['repTime_exp']-timeStore)
    squared_errors=(expRepTime_df['repTime_exp']-timeStore)**2

    new_MAE = np.mean(absolute_errors)
    new_MSE = np.mean(squared_errors)

    meanAbsoluteError = new_MAE
    meanAbsoluteError_track.append(meanAbsoluteError)

    # saving the bcs script if it has the lowest error so far
    if meanAbsoluteError == min(meanAbsoluteError_track):
        bcsMinErrorScript = modified_bcs_code

    meanSquaredError = new_MSE
    meanSquaredError_track.append(meanSquaredError)
   
    runs+=1

#------------------ 09 SAVING_SCRIPTS ------------------

# saving the bcs script from the final iteration of fitting
with open(f'{folderPath}/{bcsScriptName}_final.bc', 'w') as file:
    file.write(modified_bcs_code)

# saving the bcs script with the lowest error over the fitting iterations
with open(f'{folderPath}/{bcsScriptName}_minError.bc', 'w') as file:
    file.write(bcsMinErrorScript)    

#------------------ 10 SAVING_DATAFRAMES ------------------

# saving the simulated replication timing dataframe
timeStore_dfs_list = []

for run, timeStore in run2timeStore.items():
    df=pd.DataFrame({"repTime":timeStore})
    df["i"]=df.index
    df["run"]=run
    timeStore_dfs_list.append(df)

timeStore_df = pd.concat(timeStore_dfs_list,ignore_index=True)

timeStore_df.to_csv(f"{folderPath}/{todaysDate}timeStore_{bcsScriptName}.csv", index=False)

# saving the errors dataframe
error_df = pd.DataFrame(meanAbsoluteError_track, columns=['MAE'])
error_df.reset_index(inplace=True)
error_df.rename(columns={'index': 'run'}, inplace=True)
error_df["MSE"]=meanSquaredError_track

error_df.to_csv(f"{folderPath}/{todaysDate}error_{bcsScriptName}.csv", index=False)

print("Finished")
