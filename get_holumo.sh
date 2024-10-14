#!/bin/bash
echo "Running:" $1
Multiwfn_noGUI $1 << EOF > out.txt
0
q
EOF

#  Note: Orbital    57 is HOMO, energy:   -0.331902 a.u.   -9.031509 eV
#        Orbital    58 is LUMO, energy:   -0.005486 a.u.   -0.149290 eV
HOMO=$(grep "Orbital" out.txt | grep "HOMO" | awk '{print $7}')
LUMO=$(grep "Orbital" out.txt | grep "LUMO" | awk '{print $6}')
echo "HOMO energy: $HOMO"
echo "LUMO energy: $LUMO"
echo $HOMO $LUMO >> holumo.txt

rm -f out.txt
