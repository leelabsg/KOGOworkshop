library(qqman)

gwas<-read.table('finalresult.txt',h=T)

View(gwas)

colnames(gwas)[c(2,3,14)]<-c('BP','SNP','P')

manhattan(gwas, main='Chromosome 1 Manhattan Plot')

qq(gwas$P, main='Q-Q plot of chromosome 1')
