# This script writes Beacon Calculus (bcs) scripts with custom parameters
# To set origin firing rates based on experimental replication timings, replication timiing data from Muller et al. (2014) must be downloaded and processed.
# This creates a version of the script which records the origin each replication fork is from (allowing replicon lengths and the of active replication forks to be tracked)

##################################################################################
import pandas as pd
import datetime

todaysDate = datetime.date.today().strftime("%d%m%Y") # Formatted like DDMMYYYY

#------------------ 01 SETTING GLOBAL PARAMETERS ------------------

number_of_FF = 200 # number of firing factors

recycling_rate = 0.05 # recycling rate of firing factors

fork_speed = 1.4 # constant replication fork speed (kb/min)

#------------------ 02 SETTING ORIGIN PARAMETERS ------------------

# Importing the DataFrame containing the positions of origins (processed data from Siow et al. 2012)
filePath = "data/origin_positions.csv"
ori_df = pd.read_csv(filePath, header=0) 

# Importing the DataFrame containing the experimentally determined replication timings (processed data from Muller et al. 2014)
filePath = "data/completeExpRepTime.csv"
extRepTime_df = pd.read_csv(filePath, header=0) 

# Getting the experimentally determined replication timing of each origin
ori_df = ori_df.merge(extRepTime_df[['ch', 'i', 'expRepTime']], on=['ch', 'i'], how='left')

# Calculating the firing rate of each origin based on their experimental replication timing and the total number of firing factors in the model 
ori_df['fire'] = (1/ori_df['expRepTime'])/number_of_FF

# Adding a column for the length of each chromosome
chr_lengths = [230, 813, 317, 1532, 577, 270, 1100, 563, 440, 746, 667, 1078, 924, 784, 1091, 948] # Define the length of each chromosome (in kb)
chr2length = {chromosome: chr_lengths[chromosome - 1] for chromosome in range(1, 17)} # Create a dictionary mapping each chromosome to its length (adding 1 to include the endpoint)
ori_df['length'] = ori_df['ch'].map(chr2length)
ori_df['ori'] = ori_df['i'].copy()

# Making a dataframe containing the origin parameters to be used in the model
parameter_df = ori_df[['i','ch','length','fire','ori']].copy()

parameters_string = " || ".join(f"ORI[{','.join(map(lambda x: str(int(x)) if isinstance(x, float) and x.is_integer() else f'{x:.5f}', row))}]" for row in parameter_df.values.tolist())

#------------------ 03 PUTTING TOGETHER THE BCS SCRIPT ------------------

FF_chunk = " || FF[]"
firing_factors = FF_chunk * number_of_FF

model_string = f"""//S. cerevisiae whole genome replication model

fast = 100000; //fast rate
v = {fork_speed}; //fork velocity in kilobases per minute
           
//process definitions

FF[] = {{@factor![0],1}}.{{dwell,{recycling_rate}}}.FF[];

ORI[i,ch,length,fire,ori] = {{@factor?[0],fire}}.(FL[i,ch,length,ori]||FR[i,ch,length,ori])
              + {{ch?[i],fast}};
FR[i,ch,length,ori] = {{ch![i],fast}}.[i < length] -> {{~ch?[i+1],v}}.FR[i+1,ch,length,ori];
FL[i,ch,length,ori] = {{ch![i],fast}}.[i > 0] -> {{~ch?[i-1],v}}.FL[i-1,ch,length,ori];

//initial processes in the system
{parameters_string}{firing_factors};

"""

#------------------ 04 SAVING THE THE BCS SCRIPT ------------------

file_path = f"bcs_scripts/200FF0p05d_mapOris.bc" # example file (as shown in other scripts)
# custom name:
# file_path = f"bcs_scripts/{todaysDate}_custom_bcs_script.bc"

with open(file_path, "w") as file:
    file.write(model_string)
