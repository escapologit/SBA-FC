#!/bin/bash

input_dir="/media/ericsson/verbatim_hdd/ERIKA_TOMAKIC/SEMINAR/CopBET/LSDdata"
output_dir="/media/ericsson/verbatim_hdd/ERIKA_TOMAKIC/mri_vol2surf_gii"

hemispheres=("lh" "rh")

for subject_path in "${input_dir}"/sub-*; do
    # Extract subject ID from the path
    subject=$(basename "${subject_path}")

    # Create an output directory for the subject
    subject_output_dir="${output_dir}/${subject}"
    mkdir -p "${subject_output_dir}"

    # Loop through each session (e.g., ses-LSD, ses-PLCB)
    for session_path in "${subject_path}"/ses-*; do
        session=$(basename "${session_path}")

	session_output_dir="${output_dir}/${subject}/${session}"
	mkdir -p "${session_output_dir}"

        # Path to the functional data
        func_dir="${session_path}/func"

        # Check if func directory exists
        if [[ -d "${func_dir}" ]]; then
            # Loop through each run in the functional data
            for func_file in "${func_dir}"/*_bold.nii.gz; do
                # Extract run information from the file name
                func_filename=$(basename "${func_file}")
                run=$(echo "${func_filename}" | grep -o 'run-[0-9][0-9]')

                # Skip files without "run" in the name
                if [[ -z "${run}" ]]; then
                    continue
                fi

                for hemi in "${hemispheres[@]}"; do
                    # Output file path
                    output_file="${output_dir}/${subject}/${session}/${subject}_${session}_${run}_vol2surf_${hemi}.gii"
                    echo "$func_file" 
		    echo "$output_file"
	            echo "$subject, $hemi, $run, $session"
                    #exit 1
                    # Run the mri_vol2surf command
                    echo "Processing ${subject} ${session} ${run} ${hemi}..."
                    mri_vol2surf --mov "${func_file}" --regheader "${subject}" --hemi "${hemi}" --projfrac 0.5 --o "${output_file}"
                    
                    # Check if the command was successful
                    if [ $? -ne 0 ]; then
                        echo "Error processing ${subject} ${session} ${run} ${hemi}. Skipping..."
                        continue
                    fi
                done
            done
        else
            echo "Functional directory not found for ${subject} ${session}. Skipping..."
        fi
    done
done

echo "Processing completed!"
