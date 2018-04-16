## Tool for comparing two sets of gene families 

This tool compares genefamily sets built from two different sets of HMM databases to give family correspondences between the two sets. This tool can be used to compare two family sets built at different levels in the species tree to detect large scale duplications or deletions. This can also be used to compared families obtained from different family building tools.  

#### System requirments
* PYTHON 2.7
* _hmmpress_ and _hmmscan_ tools from HMMER package installed and add to PATH
* PYTHON modules _re_, _sys_, _os_, _subprocess_ and _argparse_


#### Tool options
```
usage: compare_family_sets.py [-h] --hmm1 SET1_HMM --hmm2 SET2_HMM --fasta
                              FASTA_FILE --n1 SET1_NAME --n2 SET2_NAME

Tool for comparing two family sets through family HMM sets

Arguments:
  -h, --help          show this help message and exit
  --hmm1 SET1_HMM     Concatenated HMM file for first set of families
  --hmm2 SET2_HMM     Concatenated HMM file for second set of families
  --fasta FASTA_FILE  Fasta file containing sequences that will be clustered
                      into both set of families
  --n1 SET1_NAME      Name for the first set of families
  --n2 SET2_NAME      Name for the second set of families
```


#### Output files

* The final output is written in files <n1>-<n2> and <n2>-<n1> that contain the correspondences between the set1 vs set2 and set2 vs set1 respectively.
* Format for: _<first-set-famid>_ _<second-set-famid>_ _number-sequences-common-between-two-families>_

#### Intermediate output files
* <n1>.famlist and <n2>.famlist contain the actual respective family sets. Format: _famid_ _sequenceid_
* <n1>.hmmtblout and <n2>.hmmtblout contain the respective hmmscan results
* <n1>.hmm.\* and <n1>.hmm.\* are the respective index files  


All the 3 input files must the present in the same working directory
	
