# KOGOworkshop


#SAIGE single variant association test

1. Codes for installation, step1, step2, drawing plots are provided in the SAIGE directory.

2. First, install SAIGE with devtools library in R. You may choose one of a,b, and c method.

**a) SAIGE is now available in conda environment**
 <pre>
 <code>

(Admin)
export PATH=~/anaconda3/bin:$PATH
source ~/anaconda3/etc/profile.d/conda.sh
conda create -n saige -c conda-forge -c bioconda "r-base>=4.0" r-saige
(Version 0.44.5 : which does not have Cauchy combination function)

conda activate saige
conda install -c conda-forge r-skat
conda install -c bioconda htslib  
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

**b) Install dependencies on R with devtools library**


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

**c) Installation by conda virtual environment on linux shell** 

Required yml file can be found in https://github.com/weizhouUMICH/SAIGE


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


3. Please copy all the example files in Example files in 

4. vcf.gz file and tbi files are needed in order to make single variant association test.
If you want to make association tests on all of the variants in the plink files you have, use the following codes

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

nohup step1_fitNULLGLMM.R \
--plinkFile=saige_example \
--phenoFile=saige_pheno.txt \
--phenoCol=y_binary \
--covarColList=x1,x2 \
--sampleIDColinphenoFile=IID \
--traitType=binary \
--outputPrefix=./step1_result --nThreads=4 \
--LOCO=FALSE \
--IsOverwriteVarianceRatioFile=TRUE &

</code>
</pre>

6. Running SAIGE- Step2 Single variant association test

<pre>
<code>

nohup step2_SPAtests.R \
--vcfFile=saige_example.vcf.gz \
--vcfFileIndex=saige_example.vcf.gz.tbi \
--vcfField=GT \
--chrom=1 \
--minMAF=0.0001 \
--minMAC=1 \
--sampleFile=sampleIDindosage.txt \
--GMMATmodelFile=step1_result.rda \
--varianceRatioFile=step1_result.varianceRatio.txt \
--SAIGEOutputFile=finalresult.txt \
--numLinesOutput=2 \
--IsOutputAFinCaseCtrl=TRUE \
--LOCO=FALSE &

</code>
</pre>


7. To make Manhattan plot and QQ plot, please use the plots.R code.
Rscript plots.R
<pre>
<code>
library(qqman)

gwas<-read.table('finalresult.txt',h=T)

colnames(gwas)[c(2,3,14)]<-c('BP','SNP','P')
png(file='manhattan_plot.png',width=1000, height=1000)
manhattan(gwas, main='Chromosome 1 Manhattan Plot')
dev.off()

png(file='qq_plot.png',width=1000, height=1000)
qq(gwas$P, main='Q-Q plot of chromosome 1')
dev.off()

</code>
</pre>

More detailed information in https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE

# Code for file preparation (gene based test)

Use the groupfile.R code in the sglee directory.

**Usage**
Rscript groupfile.R <method> <annovar output filename> <mode> <gene function> <exonic function>

ex)
Rscript groupfile.R SAIGE-GENE Annovar_output.csv manual exonic,intronic synonymous_SNV

**Method** : SKAT / SAIGE-GENE 
 Makes setID file if SKAT, groupfile if SAIGE-GENE
 
 **Mode** : Default(blank) / All / manual
If you only type method and annovar output filename, the default mode selects only nonsynonymous SNV and loss of function (frameshift deletion, frameshift insertion, startloss, stopgain, stoploss) variant from each gene.
 
The mode 'All' selects every variant in the gene
 
If you set the mode as 'manual', you must type specific gene function and exonic function of the markers.
You type gene functions with ',' and do not put space between them. Then, put space and do the same with exonic functions. You can select among the lists below.

 gene_function_list = { downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_exonic;splicing, ncRNA_intronic, ncRNA_splicing, ncRNA_UTR5, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3 }

exonic_function_list = { frameshift_deletion, frameshift_insertion, nonframeshift_deletion, nonframeshift_insertion, nonsynonymous_SNV, startloss, stopgain, stoploss, synonymous_SNV }

ex) Rscript SKAT Annovar_output.csv manual exonic,intronic,splicing frameshift_deletion,nonsynonymous_SNV 

Code can be found in sglee directory or in this github.
 
 
 
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

obj<-SKAT_Null_Model(Phenotype~X1+X2, data=FAM_Cov, out_type = 'C')

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

# Appendix. Annotation - Annovar
ANNOVAR (ANNOtate VARiation) is a bioinformatics software tool for the interpretation and prioritization of single nucleotide variants (SNVs), insertions, deletions, and copy number variants (CNVs) of a given genome. It has the ability to annotate human genomes hg18, hg19, and hg38.

![annovar](https://user-images.githubusercontent.com/73377376/97069199-74801580-1609-11eb-8775-0b07cadf878d.png)


### Example of basic workflow

1. Download annovar, perl and annotation database

    Here you can download annovar and see the available databases : <http://annovar.openbioinformatics.org/en/latest/user-guide/download/>  
    Then, download perl here <https://www.perl.org/get.html/>

    You can download the reference genome database that suits your study (Ex. hg19(GRCh37), hg38(GRCh38) .. )

   * Download the hg38 database
<pre>
<code>
        ./perl annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
        # Or if perl is already installed,
        ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
</code>
</pre>

2. Run annotate-varation.pl and proceed with the annotation work using the vcf file

<pre>
<code>
        ./perl table_annovar.pl  [VCF_FILENAME.vcf]  humandb/ --outfile [FILENAME] --buildver hg38 
        --protocol refGene --operation g –vcfinput
        # Or if perl is already installed,
        ./table_annovar.pl  [VCF_FILENAME.vcf]  humandb/ --outfile [FILENAME] --buildver hg38 
        --protocol refGene --operation g –vcfinput
</code>
</pre>  
* VCF_FILENAME.vcf : input vcf filename  
* FILENAME : name of the file you want to save  

For more information, please refer to this tutorial : <http://annovar.openbioinformatics.org/en/latest/user-guide/startup/#annotate_variationpl>

If the process above is successfully performed, the following file will be created. > [FILENAME].hg38_multianno.txt.  
The following analysis can be carried out using the multianno.txt file above.

