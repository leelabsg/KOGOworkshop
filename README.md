# KOGOworkshop

1. Codes for installation, step1, step2, drawing plots are provided in the SAIGE directory.

2. First, install SAIGE with devtools library in R. You may choose either a or b method.

0) SAIGE is now available in conda environment 
 <pre>
 <code>

(Admin)
export PATH=~/anaconda3/bin:$PATH
source ~/anaconda3/etc/profile.d/conda.sh
conda create -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige
conda activate saige
conda install -c conda-forge r-skat
conda install -c bioconda htslib

conda install -c compbiocore cget
conda install -c anaconda cmake
conda install -c bioconda savvy
(Done already)

(Users)
export PATH=../edu01/anaconda3/bin:$PATH
source ../edu01/anaconda3/etc/profile.d/conda.sh
conda activate saige

#Copy all the files needed 

cp -r ../edu01/sglee ~
cp ../edu01/plink2 ~ 
 </code>
 </pre>

a) Install dependencies on R with devtools library


R code

<pre>
<code>
install.packages('devtools')
install.packages('SKAT')
library(devtools)
devtools::install_github("leeshawn/MetaSKAT")
library(MetaSKAT); library(SKAT)
devtools::install_github("weizhouUMICH/SAIGE")

#try pip cget or cmake3 if not working 
library(SAIGE)
</code>
</pre>

b) Installation by conda virtual environment on linux shell 

Required yml file can be found in https://github.com/weizhouUMICH/SAIGE
or in /BiO/kogo/data

code on bash shell

<pre>
<code>
conda env create -f environment-RSAIGE.yml
conda activate RSAIGE

FLAGPATH=`which python | sed 's|/bin/python$||'`
export LDFLAGS="-L${FLAGPATH}/lib"
export CPPFLAGS="-I${FLAGPATH}/include"
pip install savvy

src_branch=master
repo_src_url=https://github.com/weizhouUMICH/SAIGE
git clone --depth 1 -b $src_branch $repo_src_url
R CMD INSTALL --library=path_to_final_SAIGE_library SAIGE

</code>
</pre>

path is ./SAIGE
When calling library in R, use 'library(SAIGE, lib.loc=path_to_final_SAIGE_library)'


3. Please copy all the example files in Example files in /BiO/kogo/data directory and the plink2 program in /BiO/kogo/apps/ to your directory. (cf .  cp -r /BiO/kogo/data/sglee/ ./ ,  cp /BiO/kogo/apps/plink2 ./)

4. vcf.gz file and tbi files are needed in order to make single variant association test.
If you want to make association tests on all of the variants in the plink files, use the following codes

<pre>
<code>
./plink2 --bfile saige_example --recode vcf id-paste=iid --out practice

bgzip practice.vcf

tabix -p vcf practice.vcf.gz  # create tbi index file
tabix -C practice.vcf.gz # create csi  # index file for step2 in SAIGE-GENE
</code>
</pre>

or if you want to convert vcf files to plink binary files, 

<pre>
<code>
./plink2 --vcf saige_example.vcf.gz --make-bed --out result 
</code>
</pre>

Installation of bgzip and tabix can be done with command 'conda install -c bioconda tabix' (conda environment) or 'apt-get install tabix' or follow the instructions in http://www.htslib.org/download/ (Download htslib)

5. Running SAIGE - fitting null GLMM
nohup command makes the process run in the background so that you can work in other processes.
<pre>
<code>
nohup Rscript step1_fitNULLGLMM.R \
--plinkFile=saige_example \
--phenoFile=saige_pheno.txt \
--phenoCol=y_binary \ 
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=./step1_result --nThreads=4 \
--LOCO=FALSE --IsOverwriteVarianceRatioFile=TRUE & 
</code>
</pre>

6. Running SAIGE- Step2 Single variant association test

<pre>
<code>
nohup Rscript step2_SPAtests.R \
--vcfFile=practice.vcf.gz \
--vcfFileIndex=practice.vcf.gz.tbi \
--vcfField=GT \
--chrom=1 --minMAF=0.0001 --minMAC=1 \
--sampleFile=sampleIDindosage.txt \
--GMMATmodelFile=step1_result.rda \
--varianceRatioFile=step1_result.varianceRatio.txt \
--SAIGEOutputFile=finalresult.txt --numLinesOutput=2 \
--IsOutputAFinCaseCtrl=TRUE --LOCO=FALSE &
</code>
</pre>


7. To make Manhattan plot and QQ plot, please use the plots.R code.
Rscript plots.R
More detailed information in https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE

# Code for file preparation

<pre>
<code>

#SetID file for SKAT 
SetID = Annovar_output[c(7,16)]
write.table(SetID,file='example.SetID',row.names=F,col.names=F,quote=F)

#7th column is the name of gene and 16th column is the name of the markers


#Group file for SAIGE-GENE
group=''
gene_list=unique(Annovar_output$Gene.refGene)

for(i in 1:length(gene_list)){
 if(i==1){
  group=paste0(group, gene_list[i])
 }else{
  group=paste(group, gene_list[i], sep='\n')
 }
 #index of the markers in the gene
 idx=which(Annovar_output$Gene.refGene==gene_list[i])
 
 #Match the form of markers' names
 for(j in idx){
  marker_name=paste0(Annovar_output[j,1],':',Annovar_output[j,2],'_',Annovar_output[j,4],'/',Annovar_output[j,5])
  group=paste(group,marker_name)
 }
}

#save result as group file
write.table(group,file='groupfile.txt',row.names=F,col.names=F,quote=F
 

</code>
</pre>


# SKAT (updated 2021-07-09)

For the binary phenotype, use **SKATBinary.SSD.All** function


<pre>
<code>
library(SKAT)

File.Bed<-'./Example1.bed'
File.Bim<-'./Example1.bim'
File.Fam<-'./Example1.fam'
File.SetID<-'./Example1.SetID'

#Calling covariate file (Order of the sample ID must be same with the fam file)
#When there is no covariate file, use FAM<-Read_Plink_FAM(File.Fam,Is.binary = F) 
#and object obj<-SKAT_Null_Model(y~X1+X2, out_type = 'C')
File.Cov<-'./Example1.Cov'
FAM_Cov<-Read_Plink_FAM_Cov(File.Fam,File.Cov,Is.binary = F)

#Object file for Null model
X1<-FAM_Cov$X1
X2<-FAM_Cov$X2
y<-FAM_Cov$Phenotype
obj<-SKAT_Null_Model(y~X1+X2, out_type = 'C')

# If the phenotype is binary one, out_type='D'
#When there is no covariate file, use FAM<-Read_Plink_FAM(File.Fam,Is.binary = F) 
#and object obj<-SKAT_Null_Model(y~1, out_type = 'C')

#Please set file location and name for SSD file and SSD.info file 
File.SSD<-'./Example1.SSD'
File.Info<-'./Example1.SSD.info'

#Generate and open SSD file for analysis
Generate_SSD_SetID(File.Bed,File.Bim,File.Fam,File.SetID,File.SSD,File.Info )
SSD.INFO<-Open_SSD(File.SSD,File.Info)
SSD.INFO$nSample
SSD.INFO$nSets

#Analysis
out<-SKAT.SSD.All(SSD.INFO,obj,method='SKATO')
#close SSD file 
Close_SSD()

</code>
</pre>

# SAIGE GENE (updated 2021-07-09)

**Step 0** : Creating Sparse GRM (Genomic Relationship Matrix)

<pre>
<code>
createSparseGRM.R --plinkFile=saige_gene_example \
--nThreads=4 \
--outputPrefix=step0_result \
--numRandomMarkerforSparseKin=2000 \
--relatednessCutoff=0.125
</code>
</pre>

**Step 1** : Fitting Null model 

<pre>
<code>
step1_fitNULLGLMM.R --plinkFile=saige_gene_example \
--phenoFile=saige_gene_pheno.txt \
--phenoCol=y_quantitative \
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=quantitative \
--invNormalize=TRUE \
--outputPrefix=step1_result \
--outputPrefix_varRatio=step1_result_ratio \
--sparseGRMFile=step0_result_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
--sparseGRMSampleIDFile=step0_result_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
--nThreads=4 \
--LOCO=FALSE \
--skipModelFitting=FALSE \
--IsSparseKin=TRUE \
--isCateVarianceRatio=TRUE
</code>
</pre>

**Step 2** : Gene based test

<pre>
<code>
step2_SPAtests.R --vcfFile=genotype_10markers.vcf.gz  \
--vcfFileIndex=genotype_10markers.vcf.gz.csi \
--chrom=1 \
--vcfField=GT \
--minMAF=0 \
--minMAC=0.5 \
--maxMAFforGroupTest=0.01 \
--sampleFile=samplelist.txt \
--GMMATmodelFile=step1_result.rda \
--varianceRatioFile=step1_result_ratio.varianceRatio.txt \
--SAIGEOutputFile=step2_result_gene.txt \
--numLinesOutput=1 \
--groupFile=groupFile_geneBasedtest_simulation.txt \
--sparseSigmaFile=step1_result_ratio.varianceRatio.txt_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseSigma.mtx \
--IsSingleVarinGroupTest=TRUE \
--method_to_CollapseUltraRare=absence_or_presence \
--MACCutoff_to_CollapseUltraRare=10  \
--LOCO=FALSE
</code>
</pre>


