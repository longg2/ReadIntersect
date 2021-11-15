#! /usr/bin/env sh
# These are the files and variables that will be needed
usage() { printf 'Ensemble Classification Program V0.01
	-b\tBlast Output File, default being used
	-k\tKraken Output file.
	-t\tNumber of CPU Threads to be used
	-h\tShow this help message and exit\n' 1>&2; exit 1; }

while getopts "b:k:t:h" arg; do
	case $arg in
		b)
			blastTable=${OPTARG}
			echo "${OPTARG} is your Blast Output File"
			;;
		k)
			krakenOut=${OPTARG}
			echo "${OPTARG} is your Kraken output file"
			;;
		t)
			ncores=${OPTARG}
			echo "${OPTARG} Number of CPU for Blast LCA Search"
			;;
		h | *)
			usage
			exit 0
			;;

	esac
done
exit 0

################
# Blast Output #
################
echo "Let's start with Blast"
# First let's take care of the blast file.  This'll take the longest time...
# NOTE: Assuming no taxa ids are present
echo "Pulling out the Accessions"
./BlastDefaultParser.awk $blastTable > ModBlast.table # If the default blast setting was used
#cut -f 1,2 $blastTable | sed -e "s/\.[0-9]$//"> ModBlast.table # Only useful if tab 6 was used
cut -f 2 ModBlast.table | sort | uniq > accessions.list
#
## Now we need to get their respective Taxa IDs. Code from https://www.biostars.org/p/353500/
echo "Translating the Accessions to TaxaIDs"
cat accessions.list | epost -db nuccore | esummary -db nuccore |\
      	xtract -pattern DocumentSummary -element Caption,TaxId > AccessionTaxaIds.table

# Merging the files together prior to determining the lineages
sed -i -e "s/\.[0-9]$//g" ModBlast.table # Removes the version numbers from accession IDs
join -1 2 -2 1 <(sort -k2 ModBlast.table) <(sort -k1 AccessionTaxaIds.table) | cut -f 2,3 -d " " > SeqWithTaxa.table

# Determining the Lineages.  Using taxonkit
sort SeqWithTaxa.table -k 3 > SeqWithTaxaSorted.table
cut -f 3 SeqWithTaxaSorted.table | uniq | taxonkit lineage |\ 
	taxonkit reformat -t | cut -f 1,4 | sed "s/;/\t/g" > taxonKitBlast.out

join -1 3 -2 1 SeqWithTaxaSorted.table taxonKitBlast.out | cut -f 2- -d " " |\
       	sed -e "s/ /\t/g" | sort -k 2 > LineagedBlastData.delim

# After this, we need to figure out the LCAs for each Blast Hit.
BlastLCA.R LineagedBlastData.delim $ncores BlastLCAout.tab

#################
# Kraken Output #
#################
echo "Now we'll look at Kraken"

grep -v "0$" $krakenOut | sort -k 3 > Krak.tmp # Sorting the Kraken file based on the taxa and removing unclassified reads from the analysis

cut -f 3 Krak.tmp | uniq | taxonkit lineage |\
	taxonkit reformat -t | cut -f 1,4 | sed "s/;/\t/g" > taxonKitKraken.out

join -1 3 -2 1 Krak.tmp taxonKitKraken.out | cut -f 3- -d " " | sed -e "s/ /\t/g" > KrakenTable.out
