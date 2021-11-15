#! /usr/bin/Rscript

# Loading the Libraries and installing if need be
listPackages <- c("dplyr", "pbapply", "tibble", "parallel")
newPackages <- listPackages[!(listPackages %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tibble))
#suppressPackageStartupMessages(library(taxonomizr))
suppressPackageStartupMessages(library(pbapply))

# Need to get those variables inserted
arguments <- commandArgs(trailingOnly = T)
blastTable <- arguments[1]
ncores <- arguments[2]
outName <- arguments[3]
composition <- 0.51

# Reading the blast output
cat("Reading the Data.\n")
blastTable <- as_tibble(read.delim(blastTable, header = T, stringsAsFactors = F, na.strings = ""))
#blastTable <- blastTable[!is.na(blastTable$genus),] # Less stringent than other options.  Also speeding things up!
invisible(gc())
#taxa <- read.table(taxa,head = F, stringsAsFactors = F)[,1]

cat("Splitting the dataframe into a list.  This should be quick.\n")
blastTable <- split(blastTable, f = blastTable$Sequence)

# This is where it splits into multiple cores
op <- pboptions(type = "timer")
cat("Performing a Weighted LCA Search: This can take a while \n")
WeightedLCATable <- pblapply(blastTable, cl = ncores, function(tmp){
#WeightedLCATable <- do.call(bind_rows,pblapply(blastTable, cl = ncores, function(tmp){
	if(nrow(tmp) == 1){ # If only one
		return(tmp %>% select(-c("Count", "Taxa")))
	}else if(nrow(tmp) == 0){ # If NOTHING
		return(NA)
	}else{
		# Making a weighted LCA table
		thresh <- ceiling(sum(tmp$Count) * composition)
		speciesTest <- as.data.frame(tmp %>% group_by(species) %>% summarize(Count = sum(Count), .groups = "drop_last"))
		if(any(speciesTest$Count >= thresh)){ # If species A OK
			index <- speciesTest %>% filter(Count >= thresh) %>% pull(species)
			tmpDat <- tmp %>% filter(species == index) %>% select(-c("Count", "Taxa")) %>% distinct()
			return(tmpDat)
		}else{ # We dig deeper
			genusTest <- as.data.frame(tmp %>% group_by(genus) %>% summarize(Count = sum(Count), .groups = "drop_last"))
			if(any(genusTest$Count >= thresh)){ # If genus A OK
				index <- genusTest %>% filter(Count >= thresh) %>% pull(genus)
				tmpDat <- tmp %>% filter(genus == index) %>% select(-c("Count", "Taxa", "species")) %>% distinct()
				tmpDat$species <- NA
				return(tmpDat)
			}else{ # Family?
				familyTest <- as.data.frame(tmp %>% group_by(family) %>% summarize(Count = sum(Count), .groups = "drop_last"))
				if(sum(familyTest$Count >= thresh)){
					index <- familyTest %>% filter(Count >= thresh) %>% pull(`family`)
					tmpDat <- tmp %>% filter(`family` == index) %>% select(-c("Count", "Taxa", "species","genus")) %>% distinct()
					tmpDat$species <- NA
					tmpDat$genus <- NA
					return(tmpDat)
				}else{
					return(NA)
				}
			}
		}
	}
})
#}))
cat("The Table is complete, now writing the results!\n")

WeightedLCATable <- WeightedLCATable[!(is.na(WeightedLCATable))]
final.df <- do.call(bind_rows,WeightedLCATable)

write.table(final.df, file = outName, row.names = F, col.names = T, quote = F, sep = "\t")

### Plotting Results
#library(ggplot2)
#library(ggpubr)
#library(tidyr)
#tmp <- WeightedLCATable %>% count(genus)
#tmp$genus <- sapply(tmp$genus, FUN = function(x){lastNotNa(getTaxonomy(x,"../NCBIdb.sql"))})
#tmp$genus[tmp$genus == "Unknown"] <- "Higher Level"
#p1 <- ggplot(data = tmp, aes(x = "", y = `n`, fill = genus)) + geom_col() + ylab("Contigs") + xlab("") + theme_classic() + theme(axis.ticks.x = element_blank())
#
#tmp <- WeightedLCATable %>% count(family)
#tmp$family <- sapply(tmp$family, FUN = function(x){lastNotNa(getTaxonomy(x,"../NCBIdb.sql"))})
#tmp$family[tmp$family == "Unknown"] <- "Higher Level"
#p2 <- ggplot(data = tmp, aes(x = "", y = `n`, fill = family)) + geom_col() + ylab("Contigs") + xlab("") + theme_classic() +
#	theme(axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank())
#
#ggarrange(p1,p2, ncol = 2, align = "hv")
#ggsave("~/ContigCounts.pdf", height = 4, width = 6)
#
## Getting the actual bp counts
#tmp <- WeightedLCATable %>% select(Sequence, genus) %>% separate(Sequence, into = c("Blah", "Blah2", "Blah3", "bp")) %>% select(bp, genus)
#tmp$bp <- as.numeric(tmp$bp)
#tmp <- tmp %>% group_by(genus) %>% summarize("TotalLength" = sum(bp))
#tmp$genus <- sapply(tmp$genus, FUN = function(x){lastNotNa(getTaxonomy(x,"../NCBIdb.sql"))})
#tmp$genus[tmp$genus == "Unknown"] <- "Higher Level"
#p1 <- ggplot(data = tmp, aes(x = genus, y = TotalLength, fill = genus)) + geom_col() + ylab("Total nt") + xlab("Genus") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
##ggsave("~/GenusNucleotideCount.pdf", height = 4, width =6)
#
#tmp <- WeightedLCATable %>% select(Sequence, family) %>% separate(Sequence, into = c("Blah", "Blah2", "Blah3", "bp")) %>% select(bp, family)
#tmp$bp <- as.numeric(tmp$bp)
#tmp <- tmp %>% group_by(family) %>% summarize("TotalLength" = sum(bp))
#tmp$family <- sapply(tmp$family, FUN = function(x){lastNotNa(getTaxonomy(x,"../NCBIdb.sql"))})
#tmp$family[tmp$family == "Unknown"] <- "Higher Level"
#p2 <- ggplot(data = tmp, aes(x = family, y = TotalLength, fill = family)) + geom_col() + ylab("Total nt") + xlab("Family") + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
#ggarrange(p1,p2, nrow = 2, align = "hv")
#ggsave("~/ContigNucleotideCounts.pdf", height = 6, width =8)
#
