
# This code just runs the Beacon Calculus (bcs) script

##################################################################################
import os
import time
#------------------ 01 RUNNING ------------------
start_time = time.time() # for tracking how long the simulations take to run

script_name = "200FF0p05d_fitted" # the bcs model being run
bcsScriptPath = f'file/path/to/{script_name}.bc'
numberOfsimulations = 500 # number of model simulations being run
outputName = f"{script_name}_s{numberOfsimulations}"
bcs_command = f"bin/bcs -s {numberOfsimulations} -o {outputName} {bcsScriptPath}"
os.system(bcs_command) # running the bcs script

print(f"output saved as {outputName}.simulation.bcs")
end_time = time.time()
duration = end_time - start_time
duration/=60 # converting from seconds to minutes
print(f"The script took {duration} minutes to run.")
print("Finished")
