#This script takes an input of a folder titled "sam" containing .sam.gz files for each sample after mapping to the concatenated genome

YOURDIR <- "/scratch/Lasi/"
setwd(paste(YOURDIR, "/sam", sep = ""))
library(dplyr)
library(ggplot2)
library(reshape2)

TARGETSPECIES <- "Lasius"                                                           


# function to create summary map dataframe by pasting columns together
cbind.fill <- function(...) {
  nm <- list(...)
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow))
  do.call(cbind, lapply(nm, function(x) {
    rbind(x, matrix(, n - nrow(x), ncol(x)))
  }))
}

# create list of sample IDs based on file names
reads <- as.data.frame(list.files(path = (paste(YOURDIR, "/sam", sep = "")))) # reads is all files in the /sam directory

for (filename in (1:nrow(reads))) {
  filer <- ((reads[filename, 1]))
  reads$size[filename] <- file.size(filer)
}
reads <- reads[reads$size > 1000, ] # ONLY KEEP FILES OVER 1000 BYTES (empty files mess up loop)

reads <- as.data.frame(gsub(".sam.gz*", "", reads[, 1]))
reads <- as.data.frame(gsub("catgen_*", "", reads[, 1]))
reads <- reads[,1]


# initialize dataframe for summary plot 
spec <- (c(
  "Camponotus", "Formica", "Lasius", "Leptothorax",
  "Myrmica", "Tapinoma", "Temnothorax", "Tetramorium"
))
summap <- data.frame(matrix(nrow = length(spec), ncol = 0))
rownames(summap) <- spec

# for each sample id...
for (row in reads) {
  setwd(paste(YOURDIR, "/sam", sep = ""))
  print(row)
  sam <- read.table(paste("catgen_", row, ".sam.gz", sep = ""), fill = T, row.names = NULL) # load sam file by id name

  # create column names based on .sam file output
  colnames(sam) <- c(
    "qname", "flag", "scaf", "pos", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
  )
  sam <- sam[!is.na(sam$scaf), ] # remove na scaf name alignment rows
  sam <- subset(sam, select = c(
    "qname", "flag", "scaf", "pos", "mapq", "cigar",
    "mrnm", "mpos", "isize", "seq", "qual"
  )) # only keep the first few columns

  # filter out reads that had a mapping quality less than 30
  sam.filter <- subset(sam, mapq >= 30)

  # label each alignment with the species it mapped to 			
  sam.filter.sp <- sam.filter %>%
    mutate(spec = case_when(
      startsWith(sam.filter$scaf, "Camponotus_") ~ "Camponotus",
      startsWith(sam.filter$scaf, "Formica_") ~ "Formica",
      startsWith(sam.filter$scaf, "Leptothorax_") ~ "Leptothorax",
      startsWith(sam.filter$scaf, "Myrmica_") ~ "Myrmica",
      startsWith(sam.filter$scaf, "Tapinoma_") ~ "Tapinoma",
      startsWith(sam.filter$scaf, "Temnothorax_") ~ "Temnothorax",
      startsWith(sam.filter$scaf, "Tetramorium_") ~ "Tetramorium",
      startsWith(sam.filter$scaf, "Lasius_") ~ "Lasius"
    ))


  print(table(sam.filter.sp$spec)) # prints number of alignments per genome

  # create file to check that species names were added
  uniques <- (sam.filter.sp[!duplicated(sam.filter.sp$spec), ]) # create list of each genome name


  tempdf <- data.frame()
  for (specs in (spec)) { # for each genus ...
    specprop <- as.numeric((sum(sam.filter.sp$spec == specs)) / (as.numeric(nrow(sam.filter.sp))))
    # the proportion of that genus = the sum of all rows that match the iterated genus, divided by number of rows in sam file
    tempdf1 <- as.data.frame(cbind((specs), (as.numeric(specprop)))) # make table with specs and prop
    tempdf1$V2 <- as.numeric(tempdf1$V2)
    tempdf <- rbind(tempdf, tempdf1) # create a temporary dataframe with the prop of each genus for the current sample
  }

  colnames(tempdf) <- c("spec", row) # name the temp data frame columns with spec and the current sample
  summap <- cbind.fill(summap, tempdf) # add this temp data frame to the summap data frame

  # select only reads mapping to target species, in this case Lasius    							
  onlyForm <- subset(sam.filter.sp, spec == TARGETSPECIES) # subset to only target genome alignments			             
  Freads <- as.data.frame(onlyForm$qname)
  colnames(Freads) <- row
  if (nrow(Freads) == 0) { # let me know if any sample has no reads mapping to the target genome
    print(paste(row, "has no reads left"))
    next
  }


  Freads[, 1] <- paste(Freads[, 1], "/1", sep = "") # add /1 like in original fastq file
  setwd(paste(YOURDIR, "/IDS", sep = ""))
  write.table(Freads, paste(row, "_ID.txt", sep = ""), # write a file with a list of every read ID that mapped to the target genome
    append = FALSE, sep = " ", dec = ".",
    row.names = F, col.names = F, quote = F
  )
}

setwd(YOURDIR)
# remove dupilcate cols or every other
summap <- as.data.frame(summap[, !duplicated(colnames(summap))])

# Create summary plot
summap1 <- melt(summap, id.vars = "spec")
summap1[is.na(summap1)] <- 0 # change NAs to 0
summap1$value <- as.numeric(summap1$value)

summapP <- ggplot(summap1, aes(spec, value, fill = spec)) +
  geom_jitter(aes(color = spec), size = 2, position = position_jitter(seed = 1)) +
  labs(
    title = "Summary of Competitive Mapping",
    x = "species", y = "proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.0, hjust = 1)) #+
 # geom_text(aes(label = ifelse(value > 0.05, round(value, 4), "")), hjust = 0.5, vjust = 0, position = position_jitter(seed = 1))
#summapP

#output summary to .jpeg
jpeg(paste("mapplot_summary.jpeg", sep = ""), width = 480, height = 480, quality = 100)
print(summapP)
dev.off()

write.table(summap, "summap.txt", row.names = F)
