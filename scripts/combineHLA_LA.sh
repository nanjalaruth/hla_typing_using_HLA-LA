#Extract column 3 
cut -f 3 R1_bestguess_G.txt > test

#Transpose column to row
tr "\n" "\t" < test > test1

#remove 1st column
cut -f 2-17 test1 > test2

#add 0 for the first 4 columns
for i in {1..4}; do sed -i 's/^/0\t/' test2 ; done

#add the sample name column twice
for i in {1..4}; do sed -i 's/^/"${dataset}"\t/' test3 ; done 

# Combine all the output
cat test3 >> H3A.hped
