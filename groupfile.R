args <- commandArgs(trailingOnly = TRUE)
method=args[1]
filename=args[2]
func=args[3]

func_gene=args[4]
if(!is.na(func_gene)){
  func_gene=strsplit(func_gene,split=',')
  func_gene=func_gene[[1]]
}
exonic_func_gene=args[5]
if(!is.na(exonic_func_gene)){
  exonic_func_gene=strsplit(exonic_func_gene,split=',')
  exonic_func_gene=exonic_func_gene[[1]]
  for(k in 1:length(exonic_func_gene)){
	if(grepl('_',exonic_func_gene[k])){
		exonic_func_gene[k]=gsub('_',' ',exonic_func_gene[k])
	}	
  }
}

print(func_gene)
print(exonic_func_gene)

if(method=='SAIGE-GENE'){
	Annovar_output=read.csv(filename)
	group=''
	gene_list=unique(Annovar_output$Gene.refGene)

	for(i in 1:length(gene_list)){
 		if(i==1){
			if(is.na(func)){
			  idx=which(Annovar_output$Gene.refGene==gene_list[i] & Annovar_output$Func.refGene=='exonic' & Annovar_output$ExonicFunc.refGene %in% c('frameshift deletion', 'frameshift insertion','nonsynonymous SNV', 'startloss', 'stopgain', 'stoploss'))
			}else if(func=='All'){
			  idx=which(Annovar_output$Gene.refGene==gene_list[i])
			}else if(func=='manual'){
			  idx=which(Annovar_output$Gene.refGene==gene_list[i] & Annovar_output$Func.refGene %in% func_gene & Annovar_output$ExonicFunc.refGene %in% exonic_func_gene) 
			}else{
			  print('Fix the flag')
			  break
			}
			if(length(idx)!=0){
				group=paste0(group, gene_list[i])
			}
 		}else{
 		  if(is.na(func)){
 		    idx=which(Annovar_output$Gene.refGene==gene_list[i] & Annovar_output$Func.refGene=='exonic' & Annovar_output$ExonicFunc.refGene %in% c('frameshift deletion', 'frameshift insertion','nonsynonymous SNV', 'startloss', 'stopgain', 'stoploss'))
 		  }else if(func=='All'){
 		    idx=which(Annovar_output$Gene.refGene==gene_list[i])
 		  }else if(func=='manual'){
 		    idx=which(Annovar_output$Gene.refGene==gene_list[i] & Annovar_output$Func.refGene %in% func_gene & Annovar_output$ExonicFunc.refGene %in% exonic_func_gene) 
 		  }else{
 		    print('Fix the flag')
 		    break
 		  }
 		  
  		if(length(idx)!=0){
				group=paste(group, gene_list[i], sep='\n')
			}
 		}
 		if(length(idx)!=0){
			#Match the form of markers' names
			for(j in idx){
				marker_name=paste0(Annovar_output[j,1],':',Annovar_output[j,2],'_',Annovar_output[j,4],'/',Annovar_output[j,5])
				group=paste(group,marker_name)
			}
		}
	}
		#save result as group file
	write.table(group,file='groupfile.txt',row.names=F,col.names=F,quote=F)
}else if(method=='SKAT'){

	Annovar_output<-read.csv(filename)
	if(is.na(func)){
	  idx=which(Annovar_output$Func.refGene=='exonic' & Annovar_output$ExonicFunc.refGene %in% c('frameshift deletion', 'frameshift insertion','nonsynonymous SNV', 'startloss', 'stopgain', 'stoploss'))
	}else if(func=='All'){
	  idx=1:nrow(Annovar_output)
	}else if(func=='manual'){
	  idx=which(Annovar_output$Func.refGene %in% func_gene & Annovar_output$ExonicFunc.refGene %in% exonic_func_gene) 
	}else{
	  print('Fix the flag')
	  break
	}
	Annovar_output=Annovar_output[idx,1:ncol(Annovar_output)]
	if(length(idx)==0){
		print('No markers')
	}else{
		SetID = Annovar_output[c(7,16)]
		write.table(SetID,file='result.SetID',row.names=F,col.names=F,quote=F)
	}
}

