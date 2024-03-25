# Open and read in file with formatted parameters
with open("param_inputs.txt", 'r') as param_file:

    # Remove new line chars and put into a list with each index corresponding to a job
    params = [line.rstrip() for line in param_file]

# Initialize a list that will contain the file names created for running
job_script_files = []

# Make a job script for each param set
for job_num in range(len(params)):

    # Build the filenames
    filename = "job_sub_script_" + str(job_num+1) + ".pbs"
    job_script_files.append(filename)

    # Generate a new file for each job
    with open(filename, 'w') as job:
        job.write("#!/bin/bash\n\n# Provide job name\n")
        job.write("#PBS -N run")
        job.write(str(job_num+1))
        job.write("\n")
        job.write("#PBS -l walltime=48:00:00\n# Specify queue\n#PBS -A open\n# Request memory allocation\n#PBS -l mem=128gb\n# Email when the job starts and when it terminates or aborts\n#PBS -m bea\n#PBS -M troy.lapolice@psu.edu\n\n# Move to directory where job was submitted from\ncd $PBS_O_WORKDIR\n\n# File name variables, no need to include file extentions\n# ---------------------------------------------------------\n")
        job.write("# SLiM model that will be run\n")
        job.write('model="model"')
        job.write('\n\n  # Job name attached to output files\noutput_name="run')
        job.write(str(job_num + 1))
        job.write('"\n\n')
        job.write("# Run the slim model and collect statistics on time run and memory usage, log to files\n\t# final_statistics.log is file generated which collects errors, memory and time statistics and end of job message\n\t# process_monitor.log is file generated which collects and monitors the process as it is being run, outputs the current year etc\n\n#Job\n")
        job.write("/usr/bin/time -v /storage/home/tml5905/conda/miniconda3/envs/env/bin/slim -d \"output_name='$output_name'\" ")
        job.write(params[job_num])
        job.write(" -d \"wd='/storage/home/tml5905/work/HuberLab/HunterGatherFarmerInteractions'\" $model.slim 2 > final_statistics_$output_name.log 1 > process_monitor_$output_name.log")
        job.write("\n\n# NOTE: If the job ends with err code 9 it means that the program used too much RAM and was force killed as a result of the memory not being available\n\n# NOTE: See as quoted from the SLiM Manual:\n\t# Defined constants can be of type logical, integer, float, or string; defining string constants\n\t# probably requires playing quoting games with your Un*x shell, such as:\n\t# slim -d \"foo='bar'\" test.txt\n")

# Generate a file to submit the jobs
    with open("submit_multi_jobs.sh", 'w') as sub:
        sub.write("#!/bin/bash")
        sub.write("\n\n")

        for job in job_script_files:
            sub.write("qsub ")
            sub.write(job)
            sub.write("\n")
