library(Biostrings)

window.size <- 1000 # configure the window size

# create some test yeast data, the coordinates must fit with the coverage ofc, 
# so this might not work directly:
mh_genome <- readDNAStringSet("Mh_assembly.fa")
chrom <- DNAString(mh_genome)

# compute gc content in a sliding window (as a fraction, if you want % : *100):
gc <- rowSums(letterFrequencyInSlidingView(chrom, window.size, c("G","C")))/window.size


# the coverage vector can be shorter than the sequence. 
# The length is equal to the largest end position of an alignment in the data:
pad.right <- function (rle, value=0, len) {
  if (length(rle) > len) stop("length mismatch: len must be >= length(rle)!")
  if (len > length(rle)) 
    c(rle, Rle(value, len-length(rle))) 
  else rle
}
cvg <- pad.right(cvg, len=length(gc))


# we have to process the coverage vector to have the same window size:
cvgw <- runmean(x=cvg, k=window.size, endrule = c("constant")) 
# endrule is important to get a vector of same size, see ?runmean