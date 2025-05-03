
# This code just runs the Beacon Calculus (bcs) script

##################################################################################
import os
import time
#------------------ 01 RUNNING ------------------
start_time = time.time() # for tracking how long the simulations take to run

script_name = "200FF0p05d_fitted" # the bcs model being run
bcsScriptPath = f'bcs_scripts/{script_name}.bc'
numberOfsimulations = 10 # number of model simulations being run
outputName = f"{script_name}_s{numberOfsimulations}"
bcs_command = f"bcs/bin/bcs -s {numberOfsimulations} -o bcs_output/{outputName} {bcsScriptPath}"
os.system(bcs_command) # running the bcs script

print(f"output saved as {outputName}.simulation.bcs")
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")
print("Finished")
