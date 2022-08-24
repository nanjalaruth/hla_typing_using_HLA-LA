#rename folders to remove hyphens for the hlatyping tool to work
for i in */*-*; do   mv "$i" `echo $i | sed -e 's/-//g'`; done

#State the input directory
indir=/scratch3/users/nanje/HLA-LA/output

#Create the HLA type list
for file in $(ls ${indir}); do new=${file%%.*}; cut -f 3 ${new}/hla/R1_bestguess_G.txt > test0; tr "\n" "\t" < test0 > test1; cut -f 2-17 test1 > test2; for i in {1..4}; do sed -i 's/^/0\t/' test2 ; done; for i in {1..2}; do sed -i "s/^/${new}\t/" test2 ; done; cat test2 >> H3A.hped; done

#rename files to their original state
for i in *SCD*; do   mv "$i" `echo $i | sed -e 's/SCD/SCD-/g'`; done
