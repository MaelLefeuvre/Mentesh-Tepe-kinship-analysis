# 20220824 - Mentesh Tepe TKGWV2 analysis.

## Summary
**Date of analysis:** 2022-08-24
**Aim:** Apply pairwise kinship estimation using Fernandes DM's method [TKGWV2](https://github.com/danimfernandes/tkgwv2.git), and the provided 'EUR22M' SNP dataset. See the original publication [here](https://doi.org/10.1038/s41598-021-00581-3)

---

## Notes: 
- This analysis uses a *fork* of the TKGWV2 repository, which can be found [here](https://github.com/MaelLefeuvre/tkgwv2/tree/develop). This version of the source code is in no way different from the original repository, but merely *adds* upon it, by providing with a command line interface for the *helper scripts* initially provided by Daniel. D. Fernandes. More information regarding the changes and rationale behind this fork can be found on the corresponding [Pull Request](https://github.com/danimfernandes/tkgwv2/pull/5)

## Usage
### 01. Setup the work environment
```Shell
mamba env create -f envs/TKGWV2.yml
mamba activate TKGWV2
source envs/TKGWV2-post-deploy.sh
```
---
### 02. Setup the input raw data
#### A. Symlink our preprocessed MT individuals
```Shell
mkdir 00-raw

# ---- symlink the output of our reanalysis into 00-raw
```Shell
ln -srt 00-raw ../../20220315-ART_MT-reanalysis/out/3-mapdamage/MT*/MT*.bam
ln -srt 00-raw ../../20220315-ART_MT-reanalysis/out/3-mapdamage/MT*/MT*.bam.bai
```

#### B. Download and pre-process MTT001 bam files
**Study Accession**  : PRJEB37213
**Citation**: 

> E. Skourtanioti, Y.S. Erdal, M. Frangipane, F.B. Restelli, K.A. Yener, F. Pinnock, P. Matthiae, R. Özbal, U.-D. Schoop, F. GuliyevGenomic history of Neolithic to Bronze Age Anatolia, Northern Levant, and Southern Caucasus
> Cell, 181 (2020), pp. 1158-1175.e28
> **DOI**	: https://doi.org/10.1016/j.cell.2020.04.044

```Shell
pushd 00-raw

wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR402/ERR4027348/MTT001.1240K.PE.bam
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR402/ERR4027349/MTT001.1240K.SR.bam

# ---- Merge and sort
mamba activate samtools-1.15.0
samtools merge -@ 16 MTT001.1240K.merged.bam MTT001.1240K.PE.bam MTT001.1240K.SR.bam
samtools sort -@ 16 MTT001.1240K.merged.bam > MTT001.1240K.merged.srt.bam
mamba deactivate

# ---- Remove duplicates 
mamba activate picard-2.27.4
picard MarkDuplicates -I MTT001.1240K.merged.srt.bam -O MTT001.1240K.merged.srt.rmdup.bam -M MTT001.1240K.merged.srt.rmdup.metrics --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true
mamba deactivate

# ---- Rescale PMD
mamba activate mapdamage-2.2.1
mapDamage -i MTT001.1240K.merged.srt.rmdup.bam --reference /data/mlefeuvre/datasets/reference/hs37d5/hs37d5.fa.gz --rescale --rescale-out MTT001.1240K.merged.srt.rmdup.rescaled.bam --folder MTT001-MD --verbose
mamba deactivate

# ---- Index BAM for good measure
mamba activate samtools-1.15.0
samtools index MTT001.1240K.merged.srt.rmdup.rescaled.bam
mamba deactivate 

cd - 
```
---
### 03. Run SNP downsampling, as per D.Fernandes recommandations...
```Shell
mkdir -p 01-subsampled && cd 01-subsampled
TK-helpers.py downsampleBam --downsampleN 1800000 && find . -type l -exec rm {} \;
cd - 
```

---
### 04. Run TKGWV2 on all 4 downsampled bam files.
**Support Files URL:** https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4

```shell
# ---- Download D.Fernandes support files
TKGWV2_SUPPORT_FILES_URL="https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4"
TK_SUPPORTFILES="./TKGWV2_SupportFiles_Share"

mkdir -p "$TK_SUPPORTFILES"

mamba activate gdown-4.6.0
gdown "$TKGWV2_SUPPORT_FILES_URL" -O TKGWV2_SupportFiles_Share --folder

```

```Bash
mkdir -p 02-TKGW-1000G-22M-EUR

# ---- Symlink the output of previous step
ln -srt 02-TKGW-10000G-22M-EUR 01-subsampled/*_subsampled.bam

# ---- Run TKGWV2
TK_SUPPORTFILES="/data/mlefeuvre/datasets/TKGWV2_SupportFiles_Share"
REFERENCE="data/reference/hs37d5/hs37d5.fa"
TKGWV2.py bam2plink --referenceGenome ${REFERENCE} \
                    --gwvList ${TK_SUPPORTFILES}/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed \
                    --gwvPlink ${TK_SUPPORTFILES}/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed \
                    --bamExtension .bam \
                    --minMQ 30 \
                    --minBQ 30 \
          plink2tkrelated --freqFile ${TK_SUPPORTFILES}/genomeWideVariants_hg19/1000GP3_EUR_22M_noFixed.frq
```
**Output:**

```Text
 ################################################################################
 ### TKGWV2 - An ancient DNA relatedness pipeline for ultra-low coverage data ###
 ## Version 1.0b - Released 07/2022
 #
 # [2022-10-13 14:57:50] Running 'bam2plink' on folder /data/mlefeuvre/Projects/Mentesh-Tepe/20220824-MT_TKGWV2/02-TKGW
         # BAM >> Pileup >> Text-PLINK 
         # Files to be processed:
                 MT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.bam
                 MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.bam
                 MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.bam
                 MTT001.1240K.merged.srt.rmdup_subsampled.bam
         # Arguments used:
                 --referenceGenome = /data/mlefeuvre/datasets/reference/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa
                 --gwvList = /data/mlefeuvre/datasets/TKGWV2_SupportFiles_Share/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed
                 --gwvPlink = /data/mlefeuvre/datasets/TKGWV2_SupportFiles_Share/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed
                 --bamExtension = .bam
                 --minMQ = 30
                 --minBQ = 30
 # [2022-10-13 15:00:02] All BAM files processed

 # [2022-10-13 15:00:03] Running 'plink2tkrelated' on folder /data/mlefeuvre/Projects/Mentesh-Tepe/20220824-MT_TKGWV2/02-TKGW
         # Text-PLINK >> Pairwise transposed text-PLINK >> Relatedness estimates NULL
         # Arguments used:
                 --freqFile     /data/mlefeuvre/datasets/TKGWV2_SupportFiles_Share/genomeWideVariants_hg19/1000GP3_EUR_22M_noFixed.frq

         # Estimating coefficient of relatedness Rxy for   MT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled 78557 SNPs      0.241 HRC       (1/6)
         # Estimating coefficient of relatedness Rxy for   MT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled  59884 SNPs      0.04 HRC        (2/6)
         # Estimating coefficient of relatedness Rxy for   MT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MTT001.1240K.merged.srt.rmdup_subsampled       38441 SNPs  0.051 HRC        (3/6)
         # Estimating coefficient of relatedness Rxy for   MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled  55442 SNPs      0.053 HRC       (4/6)
         # Estimating coefficient of relatedness Rxy for   MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MTT001.1240K.merged.srt.rmdup_subsampled       31137 SNPs  0.042 HRC        (5/6)
         # Estimating coefficient of relatedness Rxy for   MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled   MTT001.1240K.merged.srt.rmdup_subsampled        25330 SNPs  0.062 HRC        (6/6)
 # [2022-10-13 15:01:51] All dyads processed
 # [2022-10-13 15:01:51] Results exported to TKGWV2_Results.txt
```

**Results:**

| Sample2 | Sample1 | Used_SNPs |    HRC | counts0 | counts4 | Relationship |
| ------- | :------ | --------: | -----: | ------: | ------: | ------------ |
| MT26    | MT23    |     78557 | 0.2415 |    5358 |   73199 | 1st degree   |
| MT7     | MT23    |     59884 | 0.0401 |    5197 |   54687 | Unrelated    |
| MTT001  | MT23    |     38441 | 0.0510 |    4813 |   33628 | Unrelated    |
| MT7     | MT26    |     55442 | 0.0532 |    4728 |   50714 | Unrelated    |
| MTT001  | MT26    |     31137 | 0.0418 |    4206 |   26931 | Unrelated    |
| MTT001  | MT7     |     25330 | 0.0618 |    3453 |   21877 | Unrelated    |

### 05. Run distSimulations

```Shell
TK-helpers.py distSimulations --sampleVec ./commMT*____MT*.frq --numSimPairs 6000
```

```Text
Running: /usr/bin/env Rscript --vanilla -e 'source("/home/mlefeuvre/miniconda3/envs/TKGWV2/bin/tkgwv2-develop/helpers/distSimulations.R"); distSimulations(sampleVec=c("./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq", "./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq", "./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq", "./commMT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq", "./commMT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq", "./commMT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq"), numSimPairs=6000, freqFileHeader=FALSE)'

[1] "Generating simulated individuals based on allele frequencies from ./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq"
[1] "Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities."
[1] "Generating simulated individuals based on allele frequencies from ./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq"
[1] "Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities."
[1] "Generating simulated individuals based on allele frequencies from ./commMT23.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq"
[1] "Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities."
[1] "Generating simulated individuals based on allele frequencies from ./commMT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled.frq"
[1] "Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities."
[1] "Generating simulated individuals based on allele frequencies from ./commMT26.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq"
[1] "Warning: Some of the curves do not seem to be overlapping. It is advised to re-run the script with a larger number of simulated individuals for more accurate probabilities."
[1] "Generating simulated individuals based on allele frequencies from ./commMT7.fastq.combined.fq.prefixed.mapped.mappedonly.sorted.qF.sorted.cleaned.merged_rmdup.rescaled_subsampled____MTT001.1240K.merged.srt.rmdup_subsampled.frq"
```

### 06 - Prettify Simulation figures

```Bash
# ---- Convert figures to png
for pdf in  *__6000.pdf;do
    output=$(echo $pdf | sed -E 's/(MT[0-9]+).*____(MT{1,2}[0-9]+)\..*/\1-\2_sims_6000/');
    echo $output;
    pdftoppm $pdf $output -png;
done
```

```Bash
# ---- Concatenate figures
montage *.png -mode concatenate -geometry +2+2 MT_TKsimulation-results_6000.png
```

```Bash
# ---- Add labels to figures
function add_label() {
   x="$1"
   y="$2"
   label="$3"
   convert -font helvetica -fill black -pointsize 60 -gravity center -draw "text $x,$y $label" - png:-
}

cat MT_TKsimulation-results_6000.png \
    | add_label -1450 -1150 'MT23-MT26' \
    | add_label 0 -1150 'MT23-MT7' \
    | add_label 1450 -1150 'MT23-MTT001' \
    | add_label -1450 80 'MT26-MT7' \
    | add_label 0 80 'MT26-MTT001' \
    | add_label 1450 80 'MT7-MTT001' \
    > MT-TKsimulation-results-6000-final.png
```
