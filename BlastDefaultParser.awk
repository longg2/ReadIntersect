#!/usr/bin/awk -f 
# This needs to be able to parse the default blast file into something similar to the table.
BEGIN {ORS="\n"; OFS="\t"}
{
	if($0 ~ /Query= /){ # Getting the Query if we find it
		split($0, a, " ") # Splitting the line
		query=a[2] # Making sure we pull out the query sequence
		
	}else if($0 ~ />/){ # When there's a hit
		split($0, a, " ")
		subject=a[1]
		sub(">","", subject)
		print query,subject

	}
}
