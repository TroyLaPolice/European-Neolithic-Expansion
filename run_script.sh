#!/bin/bash


# File name variables, no need to include file extentions
# ---------------------------------------------------------
# SLiM model that will be run
model="model"
# File that collects errors, memory and time statistics and end of job message 
final_log="final_statistics"
# File that monitors the process as it is being run, outputs the current year etc
process_log="process_monitor"

# Run the slim model and collect statistics on time run and memory usage, log to files
/usr/bin/time -v slim $model.slim 2> $final_log.log 1> $process_log.log

# NOTE: If the job ends with err code 9 in the final log file, 
#         it means that the program used too much RAM and was force killed as a result of the memory not being available
