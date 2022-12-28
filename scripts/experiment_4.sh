rm -rf complex.text molecule.text new.pdb complex.normal
grep -w ATOM < $1 | awk '$3 == "C" '> $1.inpu
echo "REMARK PCI" > $1.pdb
perl arrangeR.pl $1.inpu $1.inp
awk '{printf ("%-10s %-11s %-7d %-8s %-8s %-8s\n",$1,$3,$5,$6,$7,$8)}' $1.inp > new.pdb
cat new.pdb >> $1.pdb
tail -n1 $1.inp | awk '{ print $2 "\n" 1 "\n" $5 "\n" "'"$1.pdb"'" "\n" "'"$1.out"'" "\n" 1.0}' > $1.txt
./pdbgepol <$1.txt
./gepol <$1.out > $1.output
rm -rf fort.15 $1.out $1.output $1.txt new.pdb $1.inpu  $1.inp 
sed '1,7d' fort.7 > $1.center 
rm -rf fort.8 fort.7  
sed '1d' $1.center > tmpfile
mv tmpfile $1.text
rm -rf $1.center tmpfile

grep -w ATOM < $2 | awk '$3 == "C" '> $2.inpu
echo "REMARK PCI" > $2.pdb
perl arrangeR.pl $2.inpu $2.inp
awk '{printf ("%-10s %-11s %-7d %-8s %-8s %-8s\n",$1,$3,$5,$6,$7,$8)}' $2.inp > new.pdb
cat new.pdb  >> $2.pdb
tail -n1 $2.inp | awk '{ print $2 "\n" 1 "\n" $5 "\n" "'"$2.pdb"'" "\n" "'"$2.out"'" "\n" 1.0}' > $2.txt
./pdbgepol <$2.txt
./gepol <$2.out > $2.output
rm -rf fort.15 $2.out $2.output $2.txt new.pdb $2.inpu 
sed '1,7d' fort.7 > $2.center 
rm -rf fort.7 fort.8 
sed '1d' $2.center > tmpfile
mv tmpfile $2.text
rm -rf $2.center tmpfile

cat $1.text $2.text >> $1+$2.centre 
rm -f $1.text $2.text 

cat $1 $2 >> comp.pdb
grep -w ATOM < comp.pdb | awk '$3 == "C" ' > pomc.pdb
rm -rf comp.pdb
perl arrangeR.pl pomc.pdb new.pdb
echo "REMARK PCI" > complex.pdb
awk '{printf ("%-10s %-11s %-7d %-8s %-8s %-8s\n",$1,$3,$5,$6,$7,$8)}' new.pdb > one.pdb
cat one.pdb >> complex.pdb
tail -n1 new.pdb | awk '{ print $2 "\n" 1 "\n" $5 "\n" "'"complex.pdb"'" "\n" "'"complex.out"'" "\n" 1.0}' > complex.txt
./pdbgepol <complex.txt
./gepol <complex.out > complex.output
rm -rf fort.15 complex.out complex.output new.inp complex.txt pomc.pdb complex.pdb one.pdb
sed '1,7d' fort.7 > new.center 
sed '1,5d' fort.8 > complex.normal
rm -rf fort.7 fort.8 
sed '1d' new.center > tmpfile
mv tmpfile new.text
rm -rf new.center tmpfile
awk '$9 == 0' $1+$2.centre > molecule.text
awk '$9 == 0' new.text > complex.text 
rm -rf new.text $1+$2.centre 