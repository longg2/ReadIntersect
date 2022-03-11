# Kraken2 and BlastN Ensemble Classification
The script used to detect intersecting Blast and Kraken2 reads in "Examining pathogen DNA recovery across the remains of a 14th century Italian friar (Giansante Brancorsini) infected with _Brucella melitensis_" (Hider _et al._,In Review).

# Setup
Place the three scripts in the same folder and run _ReadIntersection.sh_. This script assumes that the Blast (default output) and Kraken runs have already been completed. Please note that this script requires the Standard Kraken Output.

## Requirements
This script requires _taxonkit_ (https://github.com/shenwei356/taxonkit), _EUtils_ (https://dataguide.nlm.nih.gov/eutilities/how_eutilities_works.html) and _R_ to run properly. 
The _R_ script will automatically download the required packages if they are not already present.
