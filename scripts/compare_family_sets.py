#!/usr/bin/python

#Author: Akshay Yadav
#Version: 1.0

import re
import sys
import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description="Tool for comparing two family sets through family HMM sets")
parser._optionals.title="Arguments"
parser.add_argument('--hmm1',help="Concatenated HMM file for first set of families", required=True, dest="set1_hmm")
parser.add_argument('--hmm2',help="Concatenated HMM file for second set of families", required=True, dest="set2_hmm")
parser.add_argument('--fasta', help="Fasta file containing sequences that will be clustered into both set of families", required=True, dest="fasta_file")
parser.add_argument('--n1', help="Name for the first set of families", required="True", dest="set1_name")
parser.add_argument('--n2', help="Name for the second set of families", required="True", dest="set2_name")
args = parser.parse_args()


################################################################################
## function for hmmpress
def execute_hmmpress(famlist_hmm_fileName):
	sys.stdout.write("***** Executing hmmpress on {0} *******\n\n".format(famlist_hmm_fileName))	
	run_hmmpress=subprocess.Popen(["hmmpress",famlist_hmm_fileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmprocess_out = run_hmmpress.communicate()
	#print hmmprocess_out


def execute_hmmscan(famlist_hmmset_fileName, master_fasta_fileName, family_set_name):
	sys.stdout.write("***** Executing hmmscan on {0} vs {1} *******\n\n".format(master_fasta_fileName, famlist_hmmset_fileName))	
	run_hmmscan=subprocess.Popen(["hmmscan","--tblout",family_set_name+".hmmtblout","--noali",famlist_hmmset_fileName, master_fasta_fileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmscan_out = run_hmmscan.communicate()
	#print hmmscan_out


def print_famlist_file(set_hmmscan_tblout_fileName, famlist_set_fileName):
	assign_sequences_to_families_using_hmmscan_tblout(set_hmmscan_tblout_fileName, famlist_set_fileName)

## assigns sequences to best matching HMMs thus clustering them into corresponding family
def assign_sequences_to_families_using_hmmscan_tblout(hmmscan_tblout_fileName, famlist_set_fileName):
	hmmscan_tblout_file=open(hmmscan_tblout_fileName, "r")
	current_seq=""
	previous_seq=""
	famid_seqid_dict={}
	for line in hmmscan_tblout_file:
		line=line.rstrip()
		if(re.match(r'\#',line)):
			continue
		linearr=re.split(r'\s+', line)
		famid=linearr[0]
		current_seq=linearr[2]
		if(current_seq!=previous_seq):
			if(famid_seqid_dict.has_key(famid)):
				famid_seqid_dict[famid][current_seq]=1
			else:
				famid_seqid_dict[famid]={}
				famid_seqid_dict[famid][current_seq]=1
			previous_seq=current_seq

	print_hmmscan_predicted_families(famid_seqid_dict, famlist_set_fileName)
	hmmscan_tblout_file.close()

## prints famlist file format: <famid> <seqid>
def print_hmmscan_predicted_families(famid_seqid_dict, famlist_set_fileName):
	famlist_set_file=open(famlist_set_fileName, "w")
	for famid in famid_seqid_dict:
		for seqid in famid_seqid_dict[famid]:
			famlist_set_file.write(famid+" "+seqid+"\n")

	famlist_set_file.close()


## for reading the famlist file format: <famid> <seqid> with new sequence on each line
## reads the famlist file into a 2D dictionary(where 1st key is the family id and second key is the sequence id) and a 1D dictionary (where key is the sequence id and value is family id)
def read_famlist_file(famlist_fileName, famid_seqid_dict, seqid_dict):
	famlist_file=open(famlist_fileName,"r")
	for line in famlist_file:
		line=line.rstrip()
		linearr=re.split(r'\s+', line)
		famid=linearr[0]
		seqid=linearr[1]
		seqid_dict[seqid]=famid
		if(famid_seqid_dict.has_key(famid)):
			famid_seqid_dict[famid][seqid]=1
		else:
			famid_seqid_dict[famid]={}
			famid_seqid_dict[famid][seqid]=1
	
	
	famlist_file.close()

## family set in the first argument to the family set in the second argument
def compare_family_set_dicts(fam_seqid_dict, seqid_dict, comparison_result_outfileName):
	
	comparison_result_outfile=open(comparison_result_outfileName, "w")	
	
	for famid in fam_seqid_dict:
		set2_famid_counts={}
		other_set_famid_counts_dict={}
		get_other_set_famid_counts(fam_seqid_dict[famid], seqid_dict, other_set_famid_counts_dict)
		print_get_other_set_famid_counts_to_file(comparison_result_outfile, other_set_famid_counts_dict, famid)

	comparison_result_outfile.close()

## calculates how many different families for second argument the family in the first argument belongs to.
def get_other_set_famid_counts(seqid_dict_for_famid, seqid_dict, famid_counts_dict):
	for seqid in seqid_dict_for_famid:
		if not (seqid_dict.has_key(seqid)):
			continue
		other_set_famid=seqid_dict[seqid]
		if(famid_counts_dict.has_key(other_set_famid)):
			famid_counts_dict[other_set_famid]+=1
		else:
			famid_counts_dict[other_set_famid]=1

def print_get_other_set_famid_counts_to_file(comparison_result_outfile, other_set_famid_counts_dict, reference_famid):
	for famid in other_set_famid_counts_dict:
		comparison_result_outfile.write(reference_famid+" "+famid+" "+str(other_set_famid_counts_dict[famid])+"\n")

## wrapper for hmmpress
def hmmpress(family_set1_hmm, family_set2_hmm):
	execute_hmmpress(family_set1_hmm)
	execute_hmmpress(family_set2_hmm)

## wrapper for hmmscan
def hmmscan(family_set1_hmm, family_set2_hmm, master_fasta_fileName):
	execute_hmmscan(family_set1_hmm, master_fasta_fileName, family_set1_name)
	execute_hmmscan(family_set2_hmm, master_fasta_fileName, family_set2_name)

## wrapper for printing famlists
def print_famlists(family_set1_name, family_set2_name):
	family_set1_famlist_fileName=family_set1_name+".famlist"
	family_set2_famlist_fileName=family_set2_name+".famlist"


	set1_hmmscan_tblout_fileName=family_set1_name+".hmmtblout"
	set2_hmmscan_tblout_fileName=family_set2_name+".hmmtblout"

	print_famlist_file(set1_hmmscan_tblout_fileName, family_set1_famlist_fileName)
	print_famlist_file(set2_hmmscan_tblout_fileName, family_set2_famlist_fileName)

## wrapper for comparing families
def compare_family_sets_using_famlists(family_set1_name, family_set2_name):

	family_set1_famlist_fileName=family_set1_name+".famlist"
	family_set2_famlist_fileName=family_set2_name+".famlist"

		
	family_set1_famid_seqid_dict={}
	family_set2_famid_seqid_dict={}

	family_set1_seqid_dict={}
	family_set2_seqid_dict={}

	read_famlist_file(family_set1_famlist_fileName, family_set1_famid_seqid_dict, family_set1_seqid_dict)
	read_famlist_file(family_set2_famlist_fileName, family_set2_famid_seqid_dict, family_set2_seqid_dict)


	compare_family_set_dicts(family_set1_famid_seqid_dict, family_set2_seqid_dict, family_set1_name+"-"+family_set2_name)
	compare_family_set_dicts(family_set2_famid_seqid_dict, family_set1_seqid_dict, family_set2_name+"-"+family_set1_name)
#############################################################################################################################################

family_set1_name=args.set1_name
family_set2_name=args.set2_name

family_set1_hmm=args.set1_hmm
family_set2_hmm=args.set2_hmm

master_fasta_fileName=args.fasta_file

hmmpress(family_set1_hmm, family_set2_hmm)
hmmscan(family_set1_hmm, family_set2_hmm, master_fasta_fileName)
print_famlists(family_set1_name, family_set2_name)
compare_family_sets_using_famlists(family_set1_name, family_set2_name)

