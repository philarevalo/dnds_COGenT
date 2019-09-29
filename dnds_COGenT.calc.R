library(seqinr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]
groupfile = args[2]
group_type = args[3]
outfile = args[4]

groups = read.csv(groupfile, row.names='Strain')

is_in_diff_group <- function(s1, s2){
  result <- rep(0L, length(s1))
  for (i in seq_along(s1)){
    result[i] = groups[as.character(s1[i]), group_type] != groups[as.character(s2[i]), group_type]
  }
  return(result)
}


align_nt = read.alignment(infile, format='fasta')
kaks(align_nt) -> result
omega = as.matrix(result$ka / result$ks)


# turn into upper triangular matrix
omega[upper.tri(omega)] <- NA

# turn into long format
omega = melt(omega)

# remove NA rows
omega = na.omit(omega)

colnames(omega) = c('Strain1', 'Strain2', 'dnds')

omega = transform(omega, diff_group=is_in_diff_group(Strain1, Strain2))

write.csv(omega, outfile, row.names=FALSE)