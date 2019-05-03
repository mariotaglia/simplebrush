
salt="1"
long="25"
pkaA="7"
pkaB="7"

for j in $long ; do
echo $j

mkdir long_4_b_$j
cd long_4_b_$j

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

for z in {5..40} ; do
echo $z

mkdir dimz_$z
cd dimz_$z

echo "
# dimz #
$z
# delta #
0.5
# cuantas #
500000
# long #
$j
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
0
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

