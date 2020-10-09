# Read the data into R
small_counts <- read.table("data/small_counts.txt", header = TRUE)
ResultsTable_small <- read.table("data/ResultsTable_small.txt", header=TRUE)
# Sum the counts for each sample
sample_sums = apply(small_counts, MARGIN = 2, sum)


