#!/bin/bash

# Arguments passed from Snakefile
working_directory=$1
assemblies_list=$2  # All assemblies as a space-separated string
taxid=$3
tools=$4
db=$5

assembly_directory="data/assemblies/"

# Change to working directory
cd "$working_directory" || exit

# Iterate over each assembly
for assembly in $assemblies_list; do
    # Extract sample name from assembly path (remove directory and extension)
    assembly_name=$(basename "$assembly" .fa)

    # Create output directory structure
    mkdir -p "${working_directory}/contaminant_screening_report/fcs_vec/${assembly_name}/"
    mkdir -p "${working_directory}/contaminant_screening_report/fcs_contam/${assembly_name}/"

    # Prepare files for vector screen
    if cp "$assembly" "${tools}/fcsadaptor/input_dir"; then
        # Run vector screen
        "${tools}/fcsadaptor/run_fcsadaptor.sh" \
        --fasta-input "${assembly}" \
        --output-dir "${working_directory}/contaminant_screening_report/fcs_vec/${assembly_name}" \
        --euk

        # Prepare run for contaminant screening
        python3 "${tools}/fcsgx/fcs.py" screen genome \
        --fasta "${assembly}" \
        --out-dir "${working_directory}/contaminant_screening_report/fcs_contam/${assembly_name}" \
        --gx-db "$db" \
        --tax-id "$taxid"
    else
        echo "Failed to copy assembly file: $assembly"
    fi
done

# Create the output flag file to signal completion
touch "${working_directory}/contaminant_screening_done.flag"
