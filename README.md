MODEST: motif outlier detection system (combined with IGLOSS)

Program takes five arguments:
proteome_file_name,  
query_file_name,  
scale,  
number_of_iterations,  
threshold 

Proteome file can be in fasta or text format.

Query file consists of one or more rows of amino acids symbols ARNDCQEGHILKMFPSTWYV. Symbols 'X' and 'x' standing for any amino acid symbol are allowed.
One extra row of 0 and 1 is allowed, where 1 stands for a conserved position. All rows in the query file must be of equal length.

Scale is a positive float. This parameter measures similarity between the input query and the response. The higher it is - higher is the similarity. Reasonable scale is from 3 to 15.

Number of iterations is selfexplanitory.

Threshold is used for discarding positives. Threshold 0 means no discarding. One tested threshold is 0.001. 

MODEST.cpp calls logistic.py to make a logistic fit of the data.

handlemore.py is python script that combines more motif searches into one tsv table.
