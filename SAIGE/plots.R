library(qqman)

gwas<-read.table('finalresult.txt',h=T)

colnames(gwas)[c(2,3,14)]<-c('BP','SNP','P')
jpeg(file='manhattan_plot.jpeg',width=1000, height=1000)
manhattan(gwas, main='Chromosome 1 Manhattan Plot')
dev.off()

jpeg(file='manhattan_plot.jpeg',width=1000, height=1000)
qq(gwas$P, main='Q-Q plot of chromosome 1')
dev.off()
