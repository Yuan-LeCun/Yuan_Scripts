#!/bin/bash

# Loop through all .chk files in the current directory
for chk_file in *.chk; do
    # Extract the base name without the extension
    base_name="${chk_file%.chk}"
    # Convert .chk file to .fchk file using formchk
    if ! formchk "$chk_file" "${base_name}.fchk"; then
        echo "Error processing $chk_file with formchk" >> results.txt
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
        echo "Error processing ${base_name}.fchk with Multiwfn_noGUI" >> results.txt
        continue
    fi

    minimal=$(grep "Minimal value" out.txt | awk -F: '{print $2}' | awk '{print $1}')
    maximal=$(grep "Maximal value" out.txt | awk -F: '{print $3}' | awk '{print $1}')
    echo $minimal $maximal >> results.txt
    rm -f out.txt

done