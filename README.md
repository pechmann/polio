
Hsp90 dictates viral sequence space balancing the evolutionary tradeoffs between protein stability, aggregation and translation rate

Ron Geller, Sebastian Pechmann, Ashley Acevedo, Raul Andino & Judith Frydman


Data + computer code for Hsp90 polio paper. 

Reproduce all analyses and generate all corresponding figures from the raw input data (codon count tables):

> python reproduce_all.py 


REQUIREMENTS: 	This code is optimized to run under Pyton2.7, and requires numpy, scipy, and rpy2 libraries installed. 

FORMAT: 	The codon count tables contain the number of sequenced codons at each codon position in the coding sequence (rows)
		for each possible codon (columns are codons in alphabetical order, i.e. 'AAA', 'AAC', 'AAG', 'AAT', 'ACA', etc.)
