# Kraken2 and BlastN Ensemble Classification
The script used to detect intersecting Blast and Kraken2 reads in "Examining pathogen DNA recovery across the remains of a 14th century Italian monk (St. Brancorsini) infected with _Brucella melitensis_" (Hider _et al._,In Review).

# Setup
Place the three scripts in the same folder and run _ReadIntersection.sh_. This script assumes that the blast and Kraken runs have already been completed.

## Blast Notes
_BlastDefaultParser.awk_ is present to account for cases when the blastn output isn't tab deliminated. 

## Requirements
This script requires _taxonkit_ (https://github.com/shenwei356/taxonkit) and _R_ to run properly. The _R_ script will automatically download the required packages should they not be available.
