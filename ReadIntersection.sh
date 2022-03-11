#! /usr/bin/env bash
# These are the files and variables that will be needed
usage() { printf 'Ensemble Classification Program V0.02

	This pipline requires that taxonKit is installed and placed
	in your PATH. <https://github.com/shenwei356/taxonkit>
        -b\tBlast Output File
	-k\tKraken Output (not the report!)
	-o\tOutput Folder. Default = ReadIntersectionOut
        -t\tNumber of CPU Threads for the BlastLCA analysis. Default = 8
	-c\tMinimum Weight for BlastLCA [0,1]. Default = 0.75
	-d\t(D)o not delete the tmp folder
        -h\tShow this help message and exit\n' 1>&2; exit 1; }

log() {	printf "ReadIntersect settings for $(date):
	Log File:\t${log}
	Blast Input:\t${blastTable}
	Kraken Input:\t${krakenOut}
	Output folder:\t${out}
	CPU Threads:\t${ncores}
	-------------------------------------
	Minimum Blast Composition:\t${compo}
	-------------------------------------\n"; exit 0;
}

######################
### Default values ###
######################
compo="0.75"
ncores=8
out="ReadIntersectionOut"
scriptLocation=$(dirname $0)
delete=1
log="$(date +'%Y%m%d').log"

##############
### The UI ###
##############
while getopts "b:l:k:o:t:c:dh" arg; do
        case $arg in
                b)
                        blastTable=${OPTARG}
                        ;;
                k)
                        krakenOut=${OPTARG}
                        ;;
                o)
                        out=${OPTARG}
                        ;;
                t)
                        ncores=${OPTARG}
                        ;;
                c)
                        compo=${OPTARG}
                        ;;
                d)
                        delete=0
                        ;;
                l)
                        log=${OPTARG}
                        ;;
                h | *)
                        usage
                        exit 0
                        ;;

        esac
done
log | tee $log # The inital log file

mkdir -p $out # Since we'll be making subfolders, we'll make this first
mkdir -p $out/tmp # Intermediate files located here
mkdir -p $out/Blast

sample=$(basename $blastTable)
sample=${sample%.*}
################
# Blast Output #
################
echo "Analyzing Blast Results"
# First let's take care of the blast file.  This'll take the longest time...
echo "Pulling out the Accessions"
$scriptLocation/BlastDefaultParser.awk $blastTable > $out/tmp/ModBlast.table # If the default blast setting was used
cut -f 2 $out/tmp/ModBlast.table | sort -u --compress-program gzip > $out/tmp/accessions.list

## Now we need to get their respective Taxa IDs. Code from https://www.biostars.org/p/353500/
echo "Translating the Accessions to TaxaIDs"
cat $out/tmp/accessions.list | epost -db nuccore | esummary|\
        xtract -pattern DocumentSummary -element Caption,TaxId > $out/tmp/AccessionTaxaIds.table
 
# Merging the files together prior to determining the lineages
sed -i -e "s/\.[0-9]$//g" $out/tmp/ModBlast.table # Removes the version numbers from accession IDs

#### Sorting issues.  Can't do it in a single line...
sort -k2 --compress-program gzip $out/tmp/ModBlast.table > $out/tmp/ModBlastSorted.table
sort -k1 --compress-program gzip $out/tmp/AccessionTaxaIds.table > $out/tmp/AccessionTaxaIdsSorted.table

join -1 2 -2 1 $out/tmp/ModBlastSorted.table $out/tmp/AccessionTaxaIdsSorted.table | cut -f 2,3 -d " " | sort | uniq -c > $out/tmp/SeqWithTaxa.table
sed -i -e "s/^ *//g" -e "s/ /\t/g" $out/tmp/SeqWithTaxa.table

# Determining the Lineages.  Using taxonkit
echo "Getting the taxonomic lineages using taxonKit"
sort $out/tmp/SeqWithTaxa.table -k 3 > $out/tmp/SeqWithTaxaSorted.table
cut -f 3 $out/tmp/SeqWithTaxaSorted.table | uniq | taxonkit lineage |\
        taxonkit reformat -t | cut -f 1,4 | sed "s/;/\t/g" > $out/tmp/taxonKitBlast.out

join -1 3 -2 1 $out/tmp/SeqWithTaxaSorted.table $out/tmp/taxonKitBlast.out | cut -f 2- -d " " |\
        sed -e "s/ /\t/g" | sort -k 2 > $out/tmp/LineagedBlastData.delim

# After this, we need to figure out the LCAs for each Blast Hit.
./$scriptLocation/BlastLCA.R $out/tmp/LineagedBlastData.delim $ncores $out/Blast/${sample}Blast.tab $compo

######################
### Kraken Output ####
######################
mkdir -p $out/Kraken

echo "Now we'll look at Kraken"
grep -v "^U" $krakenOut | cut -f 2,3 | sort -k 2 > $out/tmp/Krak.tmp # Sorting the Kraken file based on the taxa and removing unclassified reads from the analysis

cut -f 2 $out/tmp/Krak.tmp | uniq | taxonkit lineage |\
	taxonkit reformat -t | cut -f 1,4 | sed "s/;/\t/g" > $out/tmp/taxonKitKraken.out

join -1 2 -2 1 $out/tmp/Krak.tmp $out/tmp/taxonKitKraken.out | cut -f 2- -d " " | sed -e "s/ /\t/g" > $out/tmp/tmpKrakenOut.tab
head -n 1 $out/Blast/${sample}Blast.tab > $out/tmp/tmpHeader.tmp
cat $out/tmp/tmpHeader.tmp $out/tmp/tmpKrakenOut.tab > $out/Kraken/${sample}Kraken.tab 

#######################################
#### Now to get the final Results. ####
#######################################
mkdir -p $out/Intersections

./$scriptLocation/FinalStep.R $out/Blast/${sample}Blast.tab $out/Kraken/${sample}Kraken.tab $out/Intersections/$sample

if [[ $delete == 1 ]]
then
	rm -rf $out/tmp/
fi

echo "ReadIntersect analysis Successful for $sample"
exit 0
