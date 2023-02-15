#!/bin/bash


# File name variables, no need to include file extentions
# ---------------------------------------------------------
# SLiM model that will be run
model="model"

# Job name attached to output files
output_name="run"

# Run the slim model and collect statistics on time run and memory usage, log to files
    # final_statistics.log is file generated which collects errors, memory and time statistics and end of job message
    # process_monitor.log is file generated which collects and monitors the process as it is being run, outputs the current year etc

# Job1
/usr/bin/time -v slim -d "output_name='$output_name'" $model.slim 2> final_statistics_$output_name.log 1> process_monitor_$output_name.log

# NOTE: If the job ends with err code 9 it means that the program used too much RAM and was force killed as a result of the memory not being available

# NOTE: See as quoted from the SLiM Manual:
    # Defined constants can be of type logical, integer, float, or string; defining string constants
    # probably requires playing quoting games with your Un*x shell, such as:
    # slim -d "foo='bar'" test.txt
