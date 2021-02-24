library(qqman)

gwas<-read.table('finalresult.txt',h=T)

colnames(gwas)[c(2,3,14)]<-c('BP','SNP','P')
png(file='manhattan_plot.png',width=1000, height=1000)
manhattan(gwas, main='Chromosome 1 Manhattan Plot')
dev.off()

png(file='manhattan_plot.png',width=1000, height=1000)
qq(gwas$P, main='Q-Q plot of chromosome 1')
dev.off()
