#!/bin/bash

# Loop through all .chk files in the current directory
for chk_file in *.chk; do
    # Extract the base name without the extension
    base_name="${chk_file%.chk}"
    # Convert .chk file to .fchk file using formchk
    formchk "$chk_file" "${base_name}.fchk"
    
    # Run Multiwfn_noGUI analysis as shown in predElectro.sh
    echo "Running:" "${base_name}.fchk"
    Multiwfn_noGUI "${base_name}.fchk" << EOF > out.txt
12
0
-1
-1
q
EOF

    minimal=$(grep "Minimal value" out.txt | awk -F: '{print $2}' | awk '{print $1}')
    maximal=$(grep "Maximal value" out.txt | awk -F: '{print $3}' | awk '{print $1}')
    echo $minimal $maximal >> results.txt
    echo " " >> results.txt
    rm -f out.txt

done