####dart2splitstree function based on Jason's dart2snapp script
### i included a function that appends NSW number and meta pop info for the splitstree label
###

dart2splitstree <- function(dart_data, basedir, species, dataset, add_pop=FALSE, pop) {

require(ape)

# Step 1, get the genotypes ready
treatment <- dart_data$treatment 
if (dart_data$encoding == "altcount") {
  cat(" Dart data object for ", dataset, "in species", species, "\n")
  cat(" Dart data object found with altcount genotype encoding. Commencing conversion to genind. \n")
} else {
  cat(" Fatal Error: The dart data object does not appear to have altcount genotype encoding. \n"); stop()
}

# genotypes -- change char for missing
gt <- dart_data$gt
gt[ is.na(gt) ] <- "?"

# make directory, write files 
dir <- paste(basedir, species, "/popgen",sep="")
if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists... content might be overwritten. \n")
}

dir <- paste(basedir, species, "/popgen/",treatment,sep="")

if(!dir.exists(dir)) {
  cat("  Directory: ", dir, " does not exist and is being created. \n")
  dir.create(dir)
} else {
  cat("  Directory: ", dir, " already exists...  \n")
}

splitstree_dir    <- paste(basedir,species,"/popgen/",treatment,"/splitstree", sep="")

if(!dir.exists(splitstree_dir)) {
  cat("  RE directory: ", splitstree_dir, " does not exist and is being created. \n")
  dir.create(splitstree_dir)
} else {
  cat("  RE directory: ", splitstree_dir, " already exists, content will be overwritten. \n")
}

snapp_nexus_file   <- paste(basedir,species,"/",species,"_SPLITSTREE_",analysis,".nex",sep="")

if(!add_pop) {
  sample_names <- rownames(gt)
} else {
  cat(   "   adding taxon names to samples \n")
  pop <-gsub(" ", "", pop, fixed = TRUE)
  sample_names <- paste(pop, rownames(gt), sep="_")
} 

n_samples <- length(sample_names)

nexus.data <- capture.output({
  taxa.labels <- as.factor(sample_names)
  n.taxa <- n_samples
  
  # write the NEXUS header
  cat('#nexus\n\n')
  
  # write the Taxa block
  cat('BEGIN Taxa;\n')
  cat('DIMENSIONS ntax=', n.taxa, ';\n', sep='')
  cat('TAXLABELS\n')
  cat(paste0("  [", seq_along(taxa.labels), "] ", taxa.labels, ""), sep='\n')
  cat(';\n')
  cat('END;\n')
  
  # write the Distances block
  cat('BEGIN Characters;\n')
  cat('DIMENSIONS nchar=', length(colnames(gt)), ';\n', sep='')
  cat('FORMAT datatype=standard\n')
  cat(' missing=?\n')
  cat(' gap=-\n')
  cat(' symbols="012"\n')
  cat(' labels\n')
  cat(' interleave\n;')

  cat('\nMATRIX\n')
  rownames(gt) <- paste0(sample_names,"\t")
  write.table(gt, row.names = T, col.names=F,quote=FALSE,sep="")
  cat(';\n')
  cat('END;\n')
})

# save the nexus file
writeLines(nexus.data, snapp_nexus_file)

return(snapp_nexus_file)
print(paste0("Output has been diverted to ",basedir,species))

}