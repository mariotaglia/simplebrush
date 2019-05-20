
salt="1 0.5 0.1 0.05 0.01 0.005 0.001 0.0001"
long="25"
pkaA="7"
pkaB="7"
pKEo="-1.1"

for j in $pKEo ; do
echo $j

mkdir long_13_c_$j
cd long_13_c_$j

for pA in $pkaA ; do
echo $pA

mkdir pkaA_$pA
cd pkaA_$pA

for pB in $pkaB ; do
echo $pB

mkdir pkaB_$pB
cd pkaB_$pB

for i in $salt ; do
echo $i

mkdir salt_$i
cd salt_$i

for z in {1..40} ; do
echo $(( 48-$z))

mkdir dimz_$(( 48-$z))
cd dimz_$(( 48-$z))



echo "
# dimz #
$(( 48-$z))
# delta #
0.5
# cuantas #
400000
# long #
25
# lseg @
0.5
# sigmaA #
0.8
# sigmaB #
0.8
# csalt #
$i
# pKaA #
$pA
# pKaB #
$pB
# pHbulk #
7.0
# st #
0.0
#pKEo#
$j
# Xulimit #
1
# infile #
0" > fort.8

../../../../../../brush

cd ..

done
cd ..

done
cd ..

done

cd ..

done

cd ..

done

