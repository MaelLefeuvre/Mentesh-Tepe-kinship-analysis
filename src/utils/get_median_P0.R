args=commandArgs(TRUE)

meansP0=read.table(args[1], header=T, sep=" ")

cat(median(meansP0$NonNormalizedP0))
