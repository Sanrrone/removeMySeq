#!/usr/bin/python3
from __future__ import with_statement
# ==============================================================================
# 						multiGenomicContext
#
# Author: Sandro Valenzuela (sandrolvalenzuead@gmail.com) 
#
# Please type "python multiGenomicContext.py -h" for usage help
#
# ==============================================================================


import sys, os, re, subprocess, csv, glob
from operator import itemgetter
from collections import deque
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, ProteinAlphabet

def which(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def isDNA(squence_str):

    for i in (squence_str):
        if str(i).upper() not in "ACTGN-RYSW":
            return False

    return True

def callblast(blastcmd):
	command=str(blastcmd+" -query "+inputseq+" -subject tmp.fasta -out tmp.out -evalue 0.001 -outfmt 10")
	#print("- running:",command)
	subprocess.call(command, shell=True)

def translateSeq(seq):
	if not isDNA(seq):
		return seq
	else:
		surplus=len(seq) % 3
		if surplus == 0:
			return Seq(seq).translate()
		else:
			return Seq(seq[0:(len(seq)-surplus)]).translate()



#######################################################################################################################
#######################################################################################################################
#######################################################################################################################
#######################################################################################################################


parser = OptionParser(usage = "Usage: python multiGenomicContext.py -q myseq.fasta -s myseqs.fasta")
parser.add_option("-q","--query",dest="inputseq",help="Your sequence you want to remove on the subject file")
parser.add_option("-s","--subject",dest="subjectseq",help="a multifasta file with sequences that contains your input sequence (or not)")
parser.add_option("-o","--output",dest="outputfile",help="output file name",default='clean.fasta')
parser.add_option("-g","--glue",dest="glueseq",help="default:5 number of genes to search upstream on the gbks",default='')
parser.add_option("-i","--identity",dest="identity",help="range 1-100 % of identity on sequence alignment to consider the gene/protein exists",default=85)
parser.add_option("-a","--alignmentLength",dest="alignL",help="range 1-100 % of aligment length to consider the gene/protein exists",default=85)
parser.add_option("-l","--length",dest="minLen",help="minimum length for resulting sequences after remove the query match sequence", default=0)
parser.add_option("-v","--inverse",dest="inverse",help="inverse the match to keep the query sequence or equivalent in subject sequences", default=False, action='store_true')
parser.add_option("-t","--translate",dest="translate",help="translate the output sequence if they are DNA", default=False, action='store_true')


(options,args) = parser.parse_args()

inputseq = options.inputseq
subjectseq = options.subjectseq
outputfile = options.outputfile
glueseq = options.glueseq
identity = int(options.identity)
alignL = int(options.alignL)
minLen = int(options.minLen)
inverse = options.inverse
translate = options.translate

#######################################################################################################################
#check variables
#######################################################################################################################
if not inputseq:
	print("ERROR: No input provided (-i), use -h for help")
	sys.exit()
else:
	if not os.path.isfile(inputseq):
		print(str("ERROR: "+inputseq+"doesn't exist, check the file directory"))
		sys.exit()

if not subjectseq:
	print("ERROR: No subject file provided (-s), use -h for help")
	sys.exit()
else:
	if not os.path.isfile(subjectseq):
		print(str("ERROR: "+subjectseq+"doesn't exist, check the file directory"))
		sys.exit()


#searching for blastp
if which('blastp') is None:
	print("ERROR: No blastp found, install it before continue")
	sys.exit()
if which('blastn') is None:
	print("ERROR: No blastn found, install it before continue")
	sys.exit()
if which('blastx') is None:
	print("ERROR: No blastx found, install it before continue")
	sys.exit()
if which('tblastn') is None:
	print("ERROR: No tblastn found, install it before continue")
	sys.exit()

#read query sequence
fasta = SeqIO.read(inputseq, "fasta")
name, sequence = fasta.id, str(fasta.seq)
seqsourcetype = "DNA" if isDNA(sequence) else "AA"
qlength=len(sequence)



#################################################################################################
########## defining strategy

outputFasta=open(outputfile,'w')
fasta_sequences = SeqIO.parse(open(subjectseq),'fasta')
for fasta in fasta_sequences:
	name, sequence = fasta.id, str(fasta.seq)
	tmp=open('tmp.fasta','w')
	tmp.write(">%s\n%s\n" % (name,sequence))
	tmp.close()
	seqtarget = "DNA" if isDNA(sequence) else "AA"
	
	##to generate tmp.out
	if seqsourcetype == "DNA":
		if seqtarget == 'DNA':
			callblast("blastn")
		else:
			callblast("blastx"),
	else:
		if seqtarget == 'DNA':
			callblast("tblastn")
		else:
			callblast("blastp")

	############
	if os.path.getsize("tmp.out")>0:
		tmp=open("tmp.out","r")
		
		
		for uniquerow in csv.reader(tmp,delimiter=','):
			#uniquerow[0] is our query protein
			#uniquerow[1] is the name of protein that match with our query (header of the fasta to be specific)
			#uniquerow[2] is identity
			#uniquerow[3] is alignment coverage (length)
			#uniquerow[8] subject sequence start
			#uniquerow[9] subject sequence end
			if float(uniquerow[2])>=identity and (float(uniquerow[3])/qlength)>=float(alignL/100.0):
				#if we are here, so, the protein exist in tmp.fasta, the next step is delete the sequence
				#call the function
				sstart=int(uniquerow[8])-1
				send=int(uniquerow[9])-1
				if inverse:
					newseq=str(sequence[sstart:send])

				else:
					newseq = str(sequence[0:sstart]+glueseq+sequence[send:len(sequence)])
					
				newseq = translateSeq(newseq) if translate else newseq
				if len(newseq) >= minLen:
					outputFasta.write(">%s\n" % (name))
					while len(newseq) > 0:
						outputFasta.write("%s\n" % (newseq[:70]))
						newseq = newseq[70:]

					outputFasta.write("\n")
				else:
					print("Warning: sequence",name,"was discarded because it is less than minimum length of "+str(minLen))
			else:
				print("Warning: sequence",name,"was discarded because it does not contain the query sequence for specified parameter of identity or coverage")
		
			
		tmp.close()
		os.remove("tmp.out")
		os.remove("tmp.fasta")


	else:
		print("Warning: No match found in subject file:",str(">"+name),"for",inputseq)
		if os.path.isfile("tmp.out"):
			os.remove("tmp.out")
			os.remove("tmp.fasta")


outputFasta.close()



