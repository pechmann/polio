######################################################################################################################################################################
#       Copyright 2017, Sebastian Pechmann, Universite de Montreal (The MIT License)
#
#       Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"),
#       to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
#       and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#       The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#       THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#       FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#       WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
######################################################################################################################################################################

Data + computer code for Hsp90 polio paper. 

Reproduce all analyses and generate all corresponding figures from the raw input data (codon count tables):

> python reprduce_all.py 


REQUIREMENTS: 	This code is optimized to run under Pyton2.7, and requires numpy, scipy, and rpy2 libraries installed. 

FORMAT: 	The codon count tables contain the number of sequenced codons at each codon position in the coding sequence (rows)
		for each possible codon (columns are codons in alphabetical order, i.e. 'AAA', 'AAC', 'AAG', 'AAT', 'ACA', etc.)
