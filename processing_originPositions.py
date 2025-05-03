
import pandas as pd

# Making a dataframe of positions of origins to include in the model

# Path to input file of origins from the OriDB (Siow et al. 2012)
oriDB_filePath = "data/oriDB_S_cerevisiae.txt"
oriDB_df = pd.read_csv(oriDB_filePath, header=0)

# Filtering to only keep 'Confirmed' or 'Likely' origins
filtered_oriDB_df = oriDB_df[oriDB_df['status'].isin(['Confirmed','Likely'])].reset_index(drop=True)

# Identifying the positions of the origins to use in the model by finding their midpoints and converting these from base pairs (bp) to kilobases (kb)
filtered_oriDB_df['i']=(filtered_oriDB_df['start']+filtered_oriDB_df['end'])//2000
filtered_oriDB_df.rename(columns={'chr': 'ch'}, inplace=True)

# Creating a new DataFrame with only the necessary information
ori_df = filtered_oriDB_df[['ch','i','status']]

# Save the final DataFrame to a CSV file
ori_df.to_csv("data/origin_positions.csv", index=False)


print("Finished")
