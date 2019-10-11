#aprox b# correr
salt=" 0.1"
long="50"
pkaA="5"
pkaB="9"
phbul="7"
pkaANa="-3"
pkaBCl="-3"
npKEo="5"
phb="7"
ll="50"
infile="0"
#sigma="0.01 0.1 0.12 0.15 0.17 0.2 0.22 0.25 0.27 0.3 0.32 0.35 0.37 0.4 0.42 0.45 0.5 0.6 0.7"
sigma="0.6 0.7 "
for j in $long ; do
echo $j

mkdir aprox_1_10o_a_$j
cd aprox_1_10o_a_$j


for pph in $phb ; do
echo $pph

mkdir ph_$pph
cd ph_$pph


for i in $salt ; do
echo $i

mkdir csalt_$i
cd csalt_$i

for pkE in $npKEo ; do
echo $pkE

mkdir pkEo_$pkE
cd pkEo_$pkE

for pA in $pkaA ; do
echo $pA

mkdir pkaA_$pA
cd pkaA_$pA

for pB in $pkaB ; do
echo $pB

mkdir pkaB_$pB
cd pkaB_$pB

for sig in $sigma ; do
echo $sig

mkdir sigma_$sig
cd sigma_$sig

for z in {1..1} ; do
echo $(( $ll-$z))

mkdir dimz_$(( $ll-$z))
cd dimz_$(( $ll-$z))
	
#cp ../../../../../pkEo_-4.5/pkaA_5/pkaB_9/sigma_0.18/dimz_$(( $ll-$z))/outin.001.dat in.txt
#cp ../../../../../../../../aprox_1_18_b_$j/ph_$pph/csalt_$i/pkEo_7/pkaA_5/pkaB_9/sigma_$sig/dimz_$(( $ll-$z))/outin.005.dat in.txt

#cp outin.001.dat in.txt
#sed -i '1d' in.txt
#sed -i $(( $ll-$z+1))'d' in.txt
#sed -i $(( 2*($ll-$z)+2))'d' in.txt
#for iii in $salti; do 
#echo $iii/ph_5/csalt_$i/pkEo_$pkE/pkaA_
#cp ../../../../../pkEo_7/pkaA_5/pkaB_9/sigma_$sig/dimz_$(( $ll-$z))/outin.006.dat in.txt

#mv outin.006.dat in.txt


echo "
# dimz #
$(( $ll-$z))
# delta #
0.5
# cuantas #
400000
# long #
$j
# lseg @
1.0
# sigmaA #
$sig
# sigmaB #
$sig
# csalt #
$i
# pKaA #
$pA
#pkaANA#
$pkaANa
# pKaB #
$pB
#pkaBCl#
$pkaBCl
# pHbulk #
$pph
# st #
0.0
#npKEo#
$npKEo
#pkeos#
5
-1
-3
-4
-5
# Xulimit #
1
# infile #
$infile" > fort.8

../../../../../../../../../brush


#done


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

cd ..

done

cd ..

done

cd ..

done

