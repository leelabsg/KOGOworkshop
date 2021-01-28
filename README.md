# KOGOworkshop

1. Codes for installation, step1, step2, drawing plots are provided in the SAIGE directory.

2. vcf.gz file and tbi files are needed in order to make single variant association test.
If you want to make association tests on all of the variants in the plink files, use the following codes

<pre>
<code>
./plink2 --bfile saige_example --recode vcf id-paste=iid --out saige_example

bgzip saige_example.vcf

tabix -p vcf saige_example.vcf.gz 
</code>
</pre>

Installation of bgzip and tabix can be done with command 'conda install -c bioconda tabix' (conda environment) or 'apt-get install tabix'

3. Running SAIGE - fitting null GLMM

<pre>
<code>
nohup Rscript step1_fitNULLGLMM.R
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

4. Running SAIGE- Step2 Single variant association test

<pre>
<code>
nohup Rscript step2_SPAtests.R
--vcfFile=genotype_10markers.missingness.vcf.gz
--vcfFileIndex=genotype_10markers.missingness.vcf.gz.tbi
--vcfField=GT
--chrom=1 --minMAF=0.0001 --minMAC=1
--sampleFile=sampleIDindosage.txt
--GMMATmodelFile=step1_result.rda
--varianceRatioFile=step1_result.varianceRatio.txt
--SAIGEOutputFile=finalresult.txt --numLinesOutput=2
--IsOutputAFinCaseCtrl=TRUE --LOCO=FALSE &
</code>
</pre>


