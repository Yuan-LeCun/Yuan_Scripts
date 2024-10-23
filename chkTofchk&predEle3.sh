#!/bin/bash

# Output file
output_file="results.txt"

# Clear the output file
> "$output_file"

# Loop through all .chk files in the current directory
for chk_file in *.chk; do
    # Extract the base name without the extension
    base_name="${chk_file%.chk}"
    # Convert .chk file to .fchk file using formchk
    if ! formchk "$chk_file" "${base_name}.fchk"; then
        echo "Error processing $chk_file with formchk" >> "$output_file"
        continue
    fi
    
    # Run Multiwfn_noGUI analysis as shown in predElectro.sh
    echo "Running:" "${base_name}.fchk"
    if Multiwfn_noGUI "${base_name}.fchk" << EOF > out.txt
12
0
-1
-1
q
EOF
    then
        echo "Successfully processed ${base_name}.fchk with Multiwfn_noGUI"
    else
        echo "Error processing ${base_name}.fchk with Multiwfn_noGUI" >> "$output_file"
        continue
    fi

    minimal=$(grep "Minimal value" out.txt | awk -F: '{print $2}' | awk '{print $1}')
    maximal=$(grep "Maximal value" out.txt | awk -F: '{print $3}' | awk '{print $1}')
    MPI=$(grep "Molecular polarity index (MPI)" out.txt | awk '{print $8}')

    rm -f out.txt

    # Run Multiwfn_noGUI analysis for HOMO and LUMO
    if Multiwfn_noGUI "${base_name}.fchk" << EOF > out.txt
0
q
EOF
    then
        HOMO=$(grep "Orbital" out.txt | grep "HOMO" | awk '{print $7}')
        LUMO=$(grep "Orbital" out.txt | grep "LUMO" | awk '{print $6}')
    else
        echo "Error processing ${base_name}.fchk for HOMO/LUMO" >> "$output_file"
        continue
    fi

    # Write minimal, maximal, MPI, HOMO, and LUMO to the output file
    echo "$minimal $maximal $MPI $HOMO $LUMO" >> "$output_file"

    rm -f out.txt
done
