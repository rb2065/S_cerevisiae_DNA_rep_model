# This script processes the experimentally determined replication timings from Muller et al. 2014.
# It uses linear interpolation to estimate the replication timings at every kb and saves the results in a DataFrame.

##################################################################################
import pandas as pd

# Converting experimental replication timings into a DataFrame

# Path to the input file containing experimental replication timing data from Muller et al. (2014)
expRepTime_filePath = "data/trep.wig"

# Read the file content and split it into chunks, one for each chromosome
with open(expRepTime_filePath, 'r') as file:
    content = file.read()
    chromosomes_list = content.split('variableStep')[1:]  # Split on 'variableStep' to isolate chromosome data

# Dictionary to store a DataFrame for each chromosome with replication timings
chr2expRepTime_df = {}

# Process each chromosome's data
for i, data in enumerate(chromosomes_list):
    chromosome_number = i + 1  # Chromosome number (1-based index)
    chr_data_list = []  # Temporary list to store replication timing data for the chromosome
    
    # Split the data into individual lines and ignore the first line (header)
    lines_list = data.splitlines()[1:]
    
    # Process each line in the chromosome data
    for line in lines_list:
        split_line = line.split("\t")  # Split the line into position and replication timing
        position = int(split_line[0])  # Position in base pairs (bp)
        position_kb = position // 1000  # Convert position from bp to kilobases (kb)
        repTime = float(split_line[1])  # Replication timing value (min)
        chr_data_list.append((position_kb, repTime))  # Append the position and timing as a tuple
    
    # Create a temporary DataFrame from the list of replication timings for the chromosome
    df1 = pd.DataFrame(chr_data_list, columns=["i", "expRepTime"])
    
    # Group by position to calculate a single replication time for each kb (since multiple entries may round to the same kb after converting from bp to kb)
    kb_data_list = [] # Temporary list to store replication timing data for the kb position
    grouped_data = df1.groupby("i")
    for kb, kb_data in grouped_data:
        mean_repTime = kb_data["expRepTime"].mean()  # Calculate the mean replication time for the kb
        kb_data_list.append((chromosome_number, kb, mean_repTime))  # Store chromosome number, position, and mean time
    
    # Create a DataFrame with aggregated replication timings for the chromosome
    df2 = pd.DataFrame(kb_data_list, columns=["ch", "i", "expRepTime"])
    chr2expRepTime_df[chromosome_number] = df2  # Add the DataFrame to the dictionary

# Filling in missing replication timings for kilobases that don't have data

# Define the length of each chromosome (in kb)
chr_lengths = [230, 813, 317, 1532, 577, 270, 1100, 563, 440, 746, 667, 1078, 924, 784, 1091, 948]

# Create a dictionary mapping each chromosome to its length (adding 1 to include the endpoint)
chr2length = {chromosome: chr_lengths[chromosome - 1] + 1 for chromosome in range(1, 17)}

# List to store DataFrames for all chromosomes with missing values filled
df_list = []

# Process each chromosome to fill in missing replication timings
for chr, repTime_df in chr2expRepTime_df.items():
    # Create a DataFrame with all positions for the chromosome
    all_positions = pd.DataFrame({
        "ch": [chr] * chr2length[chr],  # Repeat the chromosome number
        "i": range(0, chr2length[chr])  # Generate all positions from 0 to chromosome length
    })
    
    # Merge the full position DataFrame with the existing data to align positions
    merged_df = all_positions.merge(repTime_df, on=['ch', 'i'], how='left')
    
    # Interpolate missing replication timings linearly (for NaN values within the chromosome with data points either side)
    merged_df['expRepTime'].interpolate(method='linear', inplace=True)
    
    # Fill remaining NaN values at the begining and end of the chromosome with the replication time of the closest position with data
    merged_df['expRepTime'].fillna(method='ffill', inplace=True)  # Forward fill (for NaN at the chromosome end)
    merged_df['expRepTime'].fillna(method='bfill', inplace=True)  # Backward fill (for NaN values at the chromosome start)
    
    # Append the completed DataFrame to the list
    df_list.append(merged_df)

# Concatenate all chromosome DataFrames into a single DataFrame, resetting the index
completeExpRepTime_df = pd.concat(df_list, ignore_index=True)

# Save the final DataFrame to a CSV file
completeExpRepTime_df.to_csv("data/completeExpRepTime.csv", index=False)

print("Finished")