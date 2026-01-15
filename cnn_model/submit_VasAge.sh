#!/bin/bash

# Define base path variable for frequent reuse
base_path1="path for waveform data here"
base_path2="path for scripts and codes here"

# Current and working directories
curr_dir=$(pwd)
work_dir="${base_path2}"
tmp_cmd_dir="${work_dir}/tmp_cmds"

# Create working and temporary command directories if they do not exist
mkdir -p $work_dir
mkdir -p $tmp_cmd_dir

cd $work_dir

epochs_array=(300 500 200)
lr_array=(0.001 0.0001)
batch_size_array=(256 128)
patience_array=(200 100)
loss_functions=("L1" "MSE")
optimizers=("SGD" "Adam")

# Define the number of parallel jobs per submission
jobs_per_submission=4
counter=0
job_prefix="VA_Job"

# Prepare arrays for storing commands
cmds=()

# Loop over each configuration
for epochs in "${epochs_array[@]}"; do
    for patience in "${patience_array[@]}"; do
        for lr in "${lr_array[@]}"; do
            for batch_size in "${batch_size_array[@]}"; do
                for loss_function in "${loss_functions[@]}"; do
                    for optimizer_name in "${optimizers[@]}"; do

                        # Define the command with the current configuration
                        cmd="python ${base_path2}/train_run4.py \
                            --data_dir '${base_path1}/folds' \
                            --result_path '${base_path2}/loss/' \
                            --model_path '${base_path2}/trained_mdl' \
                            --lr $lr \
                            --batch-size $batch_size \
                            --epochs $epochs \
                            --optimizer_name $optimizer_name \
                            --loss_function $loss_function \
                            --patience $patience &"

                        # Add the command to the array
                        cmds+=("$cmd")
                        ((counter++))

                        # If we've reached the desired number of commands, submit them as a batch
                        if ((counter == jobs_per_submission)); then
                            job_name="${job_prefix}_${epochs}_${lr}_${batch_size}"
                            cmd_file="${tmp_cmd_dir}/${job_name}.sh"
                            {
                                echo "source ~/.bashrc"
                                echo "conda activate Vas_Age3" # Vas_Age3 is the virtual environment on HPC with all necessary packages 
                                for command in "${cmds[@]}"; do
                                    echo "$command"
                                done
                                echo "wait"
                            } >"$cmd_file"

                            # Submit the job using the CBIG submission script
                            csh ${base_path2}/CBIG_pbsubmit -cmd "bash $cmd_file" -walltime 90:00:00 -mem 32G -ncpus 4 -ngpus 1 \
                                -name "$job_name" \
                                -joberr "${base_path2}/logs/myjoberr_${job_name}.log" \
                                -jobout "${base_path2}/logs/myjobout_${job_name}.log"

                            # Reset the counter and command list
                            counter=0
                            cmds=()
                        fi
                    done
                done
            done
        done
    done
done

# Submit any remaining commands if there are any left
if [[ ${#cmds[@]} -gt 0 ]]; then
    job_name="${job_prefix}_remaining"
    cmd_file="${tmp_cmd_dir}/${job_name}.sh"
    {
        echo "source ~/.bashrc"
        echo "conda activate Vas_Age3"
        for command in "${cmds[@]}"; do
            echo "$command"
        done
        echo "wait"
    } >"$cmd_file"

    # Submit the remaining batch
    csh ${base_path2}/CBIG_pbsubmit -cmd "bash $cmd_file" -walltime 90:00:00 -mem 16G -ncpus 4 -ngpus 1 \
        -name "$job_name" \
        -joberr "${base_path2}/logs/myjoberr_${job_name}.log" \
        -jobout "${base_path2}/logs/myjobout_${job_name}.log"
fi
