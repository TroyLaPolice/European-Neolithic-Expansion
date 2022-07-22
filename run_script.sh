#!/bin/bash


# File name variables, no need to include file extentions
# ---------------------------------------------------------
# SLiM model that will be run
model="model"
# File that collects errors, memory and time statistics and end of job message
final_log="final_statistics_test"
# File that monitors the process as it is being run, outputs the current year etc
process_log="process_monitor_test"

# Run the slim model and collect statistics on time run and memory usage, log to files
/usr/bin/time -v slim -d "output_name='test'" $model.slim 2> $final_log.log 1> $process_log.log

# NOTE: If the job ends with err code 9 it means that the program used too much RAM and was force killed as a result of the memory not being available

# NOTE: See as written in SLiM Manual:
# Defined constants can be of type logical, integer, float, or string; defining string constants
    # probably requires playing quoting games with your Un*x shell, such as:
    # slim -d "foo='bar'" test.txt
