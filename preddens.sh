#!/bin/bash
echo "Running:" $1
Multiwfn $1 << EOF > out.txt
12
0
-1
-1
q
EOF

t1=$(grep "Estimated density" out.txt | awk -F : '{print $2}' | awk '{print $1}')
t2=$(grep "Product of sigma^2_tot and nu" out.txt | awk -F \( '{print $2}')
echo "M/Vm:" $t1 "g/cm^3"
echo "nu*sig2tot:" $t2 "(kcal/mol)^2"
echo "$t1 $t2" | awk '{ print "Predicted density: " a*$1+b*$2+g " g/cm^3"}' a=0.9183 b=0.0028 g=0.0443
rm -f out.txt
