# KOGOworkshop

1. Codes for installation, step1, step2, drawing plots are provided in the SAIGE directory.

2. First, install SAIGE with devtools library in R. You may choose either a or b method.

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


3. Please copy all the example files in Example files in /BiO/kogo/data directory and the plink2 program in /BiO/kogo/apps/ to your directory. (cf .  cp -r /BiO/kogo/data/sglee/ ./ cp /BiO/kogo/apps/plink2 ./)

4. vcf.gz file and tbi files are needed in order to make single variant association test.
If you want to make association tests on all of the variants in the plink files, use the following codes

<pre>
<code>
./plink2 --bfile saige_example --recode vcf id-paste=iid --out practice

bgzip practice.vcf

tabix -p vcf practice.vcf.gz 
</code>
</pre>

Installation of bgzip and tabix can be done with command 'conda install -c bioconda tabix' (conda environment) or 'apt-get install tabix' or follow the instructions in http://www.htslib.org/download/ (Download htslib)

5. Running SAIGE - fitting null GLMM
nohup command makes the process run in the background so that you can work in other processes.
<pre>
<code>
nohup /BiO/kogo/apps/R-4.0.2/bin/Rscript step1_fitNULLGLMM.R
--plinkFile=saige_example
--phenoFile=saige_pheno.txt
--phenoCol=y_binary
--covarColList=x1,x2
--sampleIDColinphenoFile=IID
--traitType=binary
--outputPrefix=./step1_result --nThreads=4
--LOCO=FALSE --IsOverwriteVarianceRatioFile=TRUE &
</code>
</pre>

6. Running SAIGE- Step2 Single variant association test

<pre>
<code>
nohup /BiO/kogo/apps/R-4.0.2/bin/Rscript step2_SPAtests.R
--vcfFile=practice.vcf.gz
--vcfFileIndex=practice.vcf.gz.tbi
--vcfField=GT
--chrom=1 --minMAF=0.0001 --minMAC=1
--sampleFile=sampleIDindosage.txt
--GMMATmodelFile=step1_result.rda
--varianceRatioFile=step1_result.varianceRatio.txt
--SAIGEOutputFile=finalresult.txt --numLinesOutput=2
--IsOutputAFinCaseCtrl=TRUE --LOCO=FALSE &
</code>
</pre>


7. To make Manhattan plot and QQ plot, please use the plots.R code.
/BiO/kogo/apps/R-4.0.2/bin/Rscript plots.R
More detailed information in https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE
