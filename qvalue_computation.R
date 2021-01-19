library(qvalue)
args = commandArgs(trailingOnly=T)

dat <- read.table(args[1], header=F, sep='\t')

qdat <- qvalue(dat$V2)

write.qvalue(qdat, file=args[2])