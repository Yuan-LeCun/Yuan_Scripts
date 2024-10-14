#!/bin/bash
echo "Running:" $1
Multiwfn_noGUI $1 << EOF > out.txt
12
0
-1
-1
q
EOF

t1=$(grep "Minimal value" out.txt |awk -F: '{print$2}' |awk '{print $1}')
t2=$(grep "Maximal value" out.txt |awk -F: '{print$3}' |awk '{print $1}')
echo $t1 $t2 >> results.txt
echo " " >> results.txt
# rm -f out.txt