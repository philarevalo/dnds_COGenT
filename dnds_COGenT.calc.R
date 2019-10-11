library(seqinr)
library(reshape2)

args = commandArgs(trailingOnly=TRUE)

infile = args[1]
groupfile = args[2]
group_type = args[3]
outfile = args[4]
prot = strsplit(basename(infile), split='.', fixed=TRUE)[[1]][[1]]

groups = read.table(groupfile, sep='\t', row.names='Strain', header=TRUE)

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
ka = as.matrix(result$ka)
ks = as.matrix(result$ks)

# turn into upper triangular matrix
omega[upper.tri(omega)] <- NA
ka[upper.tri(ka)] <- NA
ks[upper.tri(ks)] <- NA

# turn into long format
omega = melt(omega)
ka = melt(ka)
ks = melt(ks)

# remove NA rows
#omega = na.omit(omega)
#ka = na.omit(ka)
#ks = na.omit(ks)

colnames(omega) = c('Strain1', 'Strain2', 'dnds')
colnames(ka) = c('Strain1', 'Strain2', 'ka')
colnames(ks) = c('Strain1', 'Strain2', 'ks')


omega = merge(omega, ka, by=c('Strain1', 'Strain2'))
omega = merge(omega, ks, by=c('Strain1', 'Strain2'))
omega = na.omit(omega)

omega = transform(omega, diff_group=is_in_diff_group(Strain1, Strain2))
omega = transform(omega, prot=prot)

write.csv(omega, outfile, row.names=FALSE)