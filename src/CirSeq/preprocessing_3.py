import sys,gzip

workdir = sys.argv[1]

infiles = ["/9_alignment.sam.gz","/10_alignment.sam.gz"]
outfile = gzip.open(workdir + "/11_alignment.sam.gz","wb")

for file in infiles:
	infile = gzip.open(workdir + file,"rb")
	
	#running lists
	Alignment = []
	AlignmentScore = []
	
	#current read
	SequenceID = ""


	
	for line in infile:
		RawLine = line
		line = line.split()

		seqID = ""
		s1 = ""
		s2 = 0 

		seqID = line[0].split('_')

		if len(seqID) > 1:
			s1 = seqID[0]
			s2 = seqID[1]
		else:
			s1 = seqID[0]
			s2 = 0

		line2 = line[:]
		line2[0] = s1
		line2[7] = str(s2)
		NewLine = "\t".join(line2) + '\n'  

		#add alignments of the same read to the running list
		if line2[0] == SequenceID:
			#add only alignments lacking gaps and clipped bases
			if line[5].count("D") == 0 and line[5].count("I") == 0 and line[5].count("S") == 0:
				AS = line[11].split(":") 
				AlignmentScore.append(int(AS[2]))
				#Alignment.append(RawLine)
				Alignment.append(NewLine)
				
		else:
			if len(Alignment) > 0:
				#write alignment with the best alignment score
				outfile.write(Alignment[AlignmentScore.index(max(AlignmentScore))])
#			print s2, NewLine
			SequenceID = line2[0]
			Alignment = []
			AlignmentScore = []
