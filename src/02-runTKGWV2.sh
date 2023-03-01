#!/usr/bin/env

# ---- 01. Setup the work environment w/ conda or mamba
# mamba env create -f envs/TKGWV2.yml
# mamba activate TKGWV2
# bash envs/TKGWV2-post-deploy.sh

# ---- 02. Setup the raw input data

# A. Symlink our preprocessed MT individuals

# A. symlink the output of our reanalysis into 00-raw
mkdir 00-raw
ln -srt 00-raw ./out/3-mapdamage/MT*/MT*.bam
ln -srt 00-raw ./out/3-mapdamage/MT*/MT*.bam.bai

# B. Download and preprocess MTT001 bam files
#    - Study accession: PRJEB37213
#    - Citation: https://doi.org/10.1016/j.cell.2020.04.044

pushd 00-raw
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR402/ERR4027348/MTT001.1240K.PE.bam
wget http://ftp.sra.ebi.ac.uk/vol1/run/ERR402/ERR4027349/MTT001.1240K.SR.bam

# ---- Merge and sort
# *mamba activate samtools-1.15.0
samtools merge -@ 16 MTT001.1240K.merged.bam MTT001.1240K.PE.bam MTT001.1240K.SR.bam
samtools sort -@ 16 MTT001.1240K.merged.bam > MTT001.1240K.merged.srt.bam

# ---- Remove duplicates
# *mamba activate picard-2.27.4
picard MarkDuplicates -I MTT001.1240K.merged.srt.bam -O MTT001.1240K.merged.srt.rmdup.bam -M MTT001.1240K.merged.srt.rmdup.metrics --VALIDATION_STRINGENCY LENIENT --REMOVE_DUPLICATES true

# ---- Rescale PMD
# *mamba activate mapdamage-2.2.1
mapDamage -i MTT001.1240K.merged.srt.rmdup.bam --reference /data/mlefeuvre/datasets/reference/hs37d5/hs37d5.fa.gz --rescale --rescale-out MTT001.1240K.merged.srt.rmdup.rescaled.bam --folder MTT001-MD --verbose

# ---- Index BAM for good measure
# *mamba activate samtools-1.15.0
samtools index MTT001.1240K.merged.srt.rmdup.rescaled.bam
mamba deactivate
popd

# ---- 03. Run SNP downsampling as per D. Fernandes recommendations
mkdir -p 01-subsampled && pushd 01-subsampled
TK-helpers.py downsampleBam --downsampleN 1800000 && find . -type l -exec rm {} \;
popd

# ---- 04. Run TKGWV2 on all downsampled bam files files
#  Support Files URL: https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4

# Download D.Fernandes support files.
TKGWV2_SUPPORT_FILES_URL="https://drive.google.com/drive/folders/1Aw-0v_7CUorHJOLpCJ0QVCwEdH43HlO4"
TK_SUPPORTFILES="./TKGWV2_SupportFiles_Share"

mkdir -p "$TK_SUPPORTFILES"

# *mamba activate gdown-4.6.0
gdown "$TKGWV2_SUPPORT_FILES_URL" -O TKGWV2_SupportFiles_Share --folder

mkdir -p 02-TKGW-1000G-22M-EUR
# ---- Symlink the output of previous step
ln -srt 02-TKGW-10000G-22M-EUR 01-subsampled/*_subsampled.bam

# ---- Run TKGWV2
REFERENCE="./data/reference/hs37d5/hs37d5.fa"
TKGWV2.py bam2plink --referenceGenome ${REFERENCE} \
                    --gwvList ${TK_SUPPORTFILES}/genomeWideVariants_hg19/1000GP3_22M_noFixed_noChr.bed \
                    --gwvPlink ${TK_SUPPORTFILES}/genomeWideVariants_hg19/DummyDataset_EUR_22M_noFixed \
                    --bamExtension .bam \
                    --minMQ 30 \
                    --minBQ 30 \
          plink2tkrelated --freqFile ${TK_SUPPORTFILES}/genomeWideVariants_hg19/1000GP3_EUR_22M_noFixed.frq

# ---- 06. Prettify simulation figures

# Convert figures to png
for pdf in  *__6000.pdf;do
    output=$(echo $pdf | sed -E 's/(MT[0-9]+).*____(MT{1,2}[0-9]+)\..*/\1-\2_sims_6000/');
    echo $output;
    pdftoppm $pdf $output -png;
done

# Concatenate figures
montage *.png -mode concatenate -geometry +2+2 MT_TKsimulation-results_6000.png

# Add labels to figures
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
