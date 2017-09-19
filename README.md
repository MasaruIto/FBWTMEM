# FBWTMEM

FBWTMEM is an algorithm to compute maximal exact matches (MEMs) using FBWT

Please download sais from https://sites.google.com/site/yuta256/sais to ./sais

# Usage
start to make index  
./fbwtmem_seq -save index_path reference_file_path  
start to find MEM  
./fbwtmem_seq -load index_path query_file_path_1 ... [query_file_path_n]  

## Options related to making index
-save dirname      save the index  
-kmer num          set the hash key length. default is 10  
-k    num          set the size of FBWT  

## Options related to finding MEMs
-l      num        set the minimal MEM length. default is 50  
-load   dirname    load the index  
-directcompth  num set the threshold of interval size to switch computing left maximal length to directly compare with reference. default is 10  
-intervallenth num set the threshold of interval size to decide if more exact mathing in first step of algorithm. default is 10  
-skip   num        sparsify the MEM-finding algorithm. default is possible maximum value  
-print  num        print out the result when -print 1, not to do when -print 0. default is 1  
-fourcolumn     print out result by fourcolumn align  

## example
Make index
fbwtmem_seq -kmer 8 -save path ref.fa  
Finding MEMs from index  
fbwtmem_seq -l 20 -load index query.fa  


