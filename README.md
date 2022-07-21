# KOGO Workshop - Biobank Data Analysis Practice

## Connecting to the server

```
ssh edu01@59.26.46.120 -p 8517
```

## Session 1. Introduction to UK Biobank and SAIGE

### PC calculation using PLINK

```
plink2 \
--bfile ~/sglee/SAIGE/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
--pca 10 \
--out ~/sglee/SAIGE/output/pca
```

### Create Genetic Relationship Matrix (GRM)

```
plink2 \
--bfile ~/sglee/SAIGE/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
--make-king-table \
--out ~/sglee/SAIGE/output/GRM
```

### Running SAIGE using Docker

#### SAIGE

* Documentation: https://saigegit.github.io//SAIGE-doc/

SAIGE is an R package for genome-wide association tests in large-scale biobanks.

![SAIGE](https://saigegit.github.io/SAIGE-doc/assets/img/saige-siage-gene-outline.png)

#### Step 1 (Fitting the null model)

##### Binary trait

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE:/BiO/kogo/home/edu01/sglee/SAIGE wzhou88/saige:1.1.2 step1_fitNULLGLMM.R \
    --plinkFile=/BiO/kogo/home/edu01/sglee/SAIGE/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --phenoFile=/BiO/kogo/home/edu01/sglee/SAIGE/input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
    --phenoCol=y_binary \
    --covarColList=x1,x2 \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --outputPrefix=/BiO/kogo/home/edu01/sglee/SAIGE/output/example_binary \
    --nThreads=4 \
    --LOCO=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE
```

##### Quantitative trait

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE:/BiO/kogo/home/edu01/sglee/SAIGE wzhou88/saige:1.1.2 step1_fitNULLGLMM.R \
    --plinkFile=/BiO/kogo/home/edu01/sglee/SAIGE/input/nfam_100_nindep_0_step1_includeMoreRareVariants_poly \
    --phenoFile=/BiO/kogo/home/edu01/sglee/SAIGE/input/pheno_1000samples.txt_withdosages_withBothTraitTypes.txt \
    --phenoCol=y_quantitative \
    --covarColList=x1,x2 \
    --sampleIDColinphenoFile=IID \
    --traitType=quantitative \
    --invNormalize=TRUE \
    --outputPrefix=/BiO/kogo/home/edu01/sglee/SAIGE/output/example_quantitative \
    --nThreads=4 \
    --LOCO=FALSE \
    --IsOverwriteVarianceRatioFile=TRUE
```

#### Step 2 (Single-variant association test)

##### Using VCF file

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE:/BiO/kogo/home/edu01/sglee/SAIGE wzhou88/saige:1.1.2 step2_SPAtests.R \
    --vcfFile=/BiO/kogo/home/edu01/sglee/SAIGE/input/saige_example.vcf.gz \
    --vcfFileIndex=/BiO/kogo/home/edu01/sglee/SAIGE/input/saige_example.vcf.gz.csi \
    --vcfField=GT \
    --SAIGEOutputFile=/BiO/kogo/home/edu01/sglee/SAIGE/output/step2_result.txt \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=20 \
    --LOCO=FALSE \
    --GMMATmodelFile=/BiO/kogo/home/edu01/sglee/SAIGE/output/example_binary.rda \
    --varianceRatioFile=/BiO/kogo/home/edu01/sglee/SAIGE/output/example_binary.varianceRatio.txt \
    --is_output_moreDetails=TRUE
```


### Drawing Manhattan plot and Q-Q plot

Using `qqman` package in R, we can draw plots.

```
library(qqman)

# Read GWAS summary statistics
gwas <- read.table('~/sglee/SAIGE/output/step2_result.txt', header=T)

# Draw Manhattan plot
png(file='~/sglee/SAIGE/output/manhattan_plot.png', width=1200, height=600)
manhattan(gwas, bp='POS', snp='MarkerID', p='p.value', main='Manhattan Plot')
dev.off()

# Draw Q-Q plot
png(file='~/sglee/SAIGE/output/qq_plot.png', width=1000, height=1000)
qq(gwas$p.value, main='Q-Q plot')
dev.off()
```

Using `scp` command, you can move PNG files to your local machine.

```
scp -rP 8517 'edu01@59.26.46.120:~/sglee/SAIGE/output/*.png' .
```

### Resources

* UKB WES 200k PheWeb: https://ukb-200kexome.leelabsg.org/
* KoGES PheWeb: https://koges.leelabsg.org/

### FUMA

* FUMA (https://fuma.ctglab.nl/): platform that can be used to annotate, prioritize, visualize and interpret GWAS results.

## Session 2. Gene-based association test (SAIGE-GENE+)

### SKAT

Using `SKAT` package in R, we can conduct a gene-based association test.

```
library(SKAT)

File.Bed <- './Example1.bed'
File.Bim <- './Example1.bim'
File.Fam <- './Example1.fam'
File.SetID <- './Example1.SetID'

# Calling covariate file (Order of the sample ID must be same with the fam file)
# When there is no covariate file, use FAM<-Read_Plink_FAM(File.Fam, Is.binary = F) 
# and object obj<-SKAT_Null_Model(y~X1+X2, out_type = 'C')
File.Cov <- './Example1.Cov'
FAM_Cov <- Read_Plink_FAM_Cov(File.Fam, File.Cov, Is.binary = F)

# Object file for Null model
obj <- SKAT_Null_Model(Phenotype~X1+X2, data=FAM_Cov, out_type = 'C')

# If the phenotype is binary one, out_type='D'
# When there is no covariate file, use FAM<-Read_Plink_FAM(File.Fam, Is.binary = F) 
# and object obj <- SKAT_Null_Model(y~1, out_type = 'C')
# Please set file location and name for SSD file and SSD.info file 
File.SSD<-'./Example1.SSD'
File.Info<-'./Example1.SSD.info'

# Generate and open SSD file for analysis
Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID, File.SSD, File.Info)
SSD.INFO <- Open_SSD(File.SSD, File.Info)
SSD.INFO$nSample
SSD.INFO$nSets

# Analysis
out <- SKAT.SSD.All(SSD.INFO, obj, method='SKATO')
out$results

# close SSD file 
Close_SSD()
```

### SAIGE-GENE+ using Docker

#### Step 0 (Creating a sparse GRM)

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE-GENE:/BiO/kogo/home/edu01/sglee/SAIGE-GENE wzhou88/saige:1.1.2 createSparseGRM.R \
    --plinkFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/saige_gene_example \
    --nThreads=4 \
    --outputPrefix=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/sparseGRM \
    --numRandomMarkerforSparseKin=2000 \
    --relatednessCutoff=0.125
```

#### Step 1 (Fitting the null model)

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE-GENE:/BiO/kogo/home/edu01/sglee/SAIGE-GENE wzhou88/saige:1.1.2 step1_fitNULLGLMM.R \
    --sparseGRMFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --plinkFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/saige_gene_example \
    --useSparseGRMtoFitNULL=TRUE \
    --phenoFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/saige_gene_pheno.txt \
    --phenoCol=y_binary \
    --covarColList=x1,x2 \
    --qCovarColList=x2 \
    --sampleIDColinphenoFile=IID \
    --traitType=binary \
    --isCateVarianceRatio=TRUE \
    --outputPrefix=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/step1_result \
    --IsOverwriteVarianceRatioFile=TRUE
```

#### Step 2 (Gene-based association test)

```
docker run -v /BiO/kogo/home/edu01/sglee/SAIGE-GENE:/BiO/kogo/home/edu01/sglee/SAIGE-GENE wzhou88/saige:1.1.2 step2_SPAtests.R \
    --vcfFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/genotype_100markers.vcf.gz \
    --vcfFileIndex=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/genotype_100markers.vcf.gz.csi \
    --vcfField=GT \
    --SAIGEOutputFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/step2_result.txt \
    --LOCO=FALSE \
    --chrom=1 \
    --minMAF=0 \
    --minMAC=0.5 \
    --sampleFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/samplelist.txt \
    --GMMATmodelFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/step1_result.rda \
    --varianceRatioFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/step1_result.varianceRatio.txt \
    --sparseGRMFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx \
    --sparseGRMSampleIDFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/output/sparseGRM_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt \
    --groupFile=/BiO/kogo/home/edu01/sglee/SAIGE-GENE/input/groupfile.txt \
    --annotation_in_groupTest=lof,missense:lof,missense:lof:synonymous \
    --maxMAF_in_groupTest=0.0001,0.001,0.01
```

### DNAnexus

DNAnexus provides a cloud-based data analysis and management platform for DNA sequence data. (https://www.dnanexus.com/) \
You can also run SAIGE/SAIGE-GENE+ on DNAnexus.

More details can be found at: https://saigegit.github.io//SAIGE-doc/docs/UK_Biobank_WES_analysis.html

### Appendix 1 (Annotation using ANNOVAR)

#### ANNOVAR

ANNOVAR (ANNOtate VARiation) is a bioinformatics software tool for the interpretation and prioritization of single nucleotide variants (SNVs), insertions, deletions, and copy number variants (CNVs) of a given genome. It has the ability to annotate human genomes hg18, hg19, and hg38.

![annovar](https://user-images.githubusercontent.com/73377376/97069199-74801580-1609-11eb-8775-0b07cadf878d.png)

#### Basic Workflow

1. Download ANNOVAR, perl and annotation database
    * ANNOVAR: http://annovar.openbioinformatics.org/en/latest/user-guide/download/
    * Perl: https://www.perl.org/get.html
    * Download annotation database: 
        ```
        ./annotate_variation.pl -buildver hg38 -downdb -webfrom annovar refGene humandb/
        ```

2. Make a table with annotation using the VCF file
    ```
    ./table_annovar.pl [VCF_FILENAME.vcf] humandb/ --outfile [FILENAME] --buildver hg38 --protocol refGene --operation g –vcfinput
    ```

If the process above is successfully performed, a file named `[FILENAME].hg38_multianno.txt` will be created. 

### Appendix 2 (Making a group file)

You can manually generate a group file for SAIGE-GENE or SetID file for SKAT using `groupfile.R`.

**Usage**: `Rscript groupfile.R [method] [annotation_file] [mode] [gene_function_list] [exonic_function_list]`

#### Options

If you only type method and annovar output filename, the default mode selects only nonsynonymous SNV and loss of function (frameshift deletion, frameshift insertion, startloss, stopgain, stoploss) variant from each gene.

The mode `All` selects every variant in the gene

If you set the mode as `manual`, you must type specific gene function and exonic function of the markers. You type gene functions with ',' and do not put space between them. \
Then, put space and do the same with exonic functions. 

You can select among the lists below.

`mode = {, All, manual}`

`gene_function_list = {downstream, exonic, exonic;splicing, intergenic, intronic, ncRNA_exonic, ncRNA_exonic;splicing, ncRNA_intronic, ncRNA_splicing, ncRNA_UTR5, splicing, upstream, upstream;downstream, UTR3, UTR5, UTR5;UTR3}`

`exonic_function_list = {frameshift_deletion, frameshift_insertion, nonframeshift_deletion, nonframeshift_insertion, nonsynonymous_SNV, startloss, stopgain, stoploss, synonymous_SNV}`

#### Example usage

```
Rscript groupfile.R SAIGE-GENE Annovar_output.csv manual exonic,intronic synonymous_SNV
Rscript groupfile.R SKAT Annovar_output.csv manual exonic,intronic,splicing frameshift_deletion,nonsynonymous_SNV
```

