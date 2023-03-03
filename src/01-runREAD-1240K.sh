#!/usr/bin/env bash
# ------------------------------------------------------------------------------------------------ #
# 01-runREAD-1240K.sh                                                                              #
# ------------------------------------------------------------------------------------------------ #
# Description : Merges our fileset with the Reich 1240K compendium dataset, performs maf filering  #
#               and the re-extracts our samples, along with a set of 'proxy individuals' to run a  #
#               kinship analysis using READ.                                                       # 
#               Proxy individuals are used to compute the normalization value of READ and thus act #
#               as a 'surrogate population' for our Mentesh Tepe individuals.                      #
#               Note that the normalization value, and the **actual** kinship estimation of our    #
#               samples are done separately, in order to mitigate the risk of introducing biases   #
# ------------------------------------------------------------------------------------------------ #
# Dependencies: - Software:  READ | python2.7 (for READ) | plink-1.9 | eigensoft                   #
#               - Libraries: libgsl-dev libopenblas-dev liblapack-dev liblapacke-dev               #
# ------------------------------------------------------------------------------------------------ #
# Nodes       : - The script assumes you have an available conda v4.6+ in your path, and that the  #
#                 environments mentionned above (see the 'envs' directoryà are preinstalled        #
# 
#               - These dependencies will be installed in the './methods' directory and assumes    #
#                 you're working on a x86_64 architecture                                          #
# ------------------------------------------------------------------------------------------------ #
# Usage       : ./src/01-runREAD-1240K.sh                                                          #
# ------------------------------------------------------------------------------------------------ #


# -------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- PARAMETERS ------------------------------------------------------- #

MT_DATA="out/5-pileupCaller/AllMT_pileupCaller_RBq30Q30"
REICH_DATA="data/target_positions/Reich-dataset/v44.3_1240K_public"

REGEX="ART" # KRD / ART
MAF=0.05

PREPROCESS_DIR="./out/6-preprocess"

PYTHON=`which python2`
READ_HOME="`pwd`/methods/read"
PLINK_HOME="`pwd`/methods/plink1.9"
EIGENSOFT_HOME="`pwd`/methods/eigensoft"

PROCS=`nproc` 

mkdir -p "${PREPROCESS_DIR}"
mkdir -p "./in"

# -------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------- FUNCTIONS -------------------------------------------------------- #

function log(){
    # Basic logger, with 3 levels : INFO  - informational messages
    #                               WARN  - not fatal, but user should still get notified...
    #                               FATAL - Something went very wrong and should get fixed...
    local level message
    read -r level message <<< $(echo $1 $2)
    local reset='\e[97m'
    declare -A levels=([INFO]="\e[32m" [WARN]="\e[33m" [FATAL]="\e[31m")
    local color=${levels["${level}"]}    
    echo -e "${color}[${level}]: ${message}${reset}"
}

function abort(){
    # Sends a log with FATAL level, accompanied with the line number, exit code and an optional message. 
    # Then, exit.. 
    local exit_code=$?
    local message=$1
    log FATAL "Line ${LINENO}: $message (exit code = ${exit_code})"
    exit $exit_code
}

function drawline(){
    # Draws a horizontal line of dashes, spanning the entire console length.
    printf '%*s\n' "${COLUMNS:-$(tput cols)}" '' | tr ' ' -
}

function header(){
   # Calls drawline(), followed by a centered title
   local title=$1
   drawline
   title_len="${#title}"
   COLUMNS=$(tput cols)
   printf "%*s\n" "$(((title_len+COLUMNS)/2))" "$title"
}

function make_keep_file(){
    # Matches individuals within a .tfam/.fam/.peding file 
    # Writes a file containing samples for the --keep argument of plink. 
    local infile=$1
    local regex=$2
    local keep_file=$3
    local valid_suffixes="fam|tfam|pedind"
    log INFO  "Creating keep file..."
    CMD="grep -E \"${regex}\" \$(find \$(dirname ${infile})/ -name \"\$(basename ${infile}).*\" | grep -E \".*(${valid_suffixes})\") | awk '{print \$1,\$2}'  > ${keep_file}"
    echo " - [COMMAND]: ${CMD}"
    eval $CMD
}

function prime_READ(){
    # READ is unable to handle remote execution or I/O management...
    # We have to manually create the meansP0_AncientDNA_normalized file
    # and instantiate a symbolink link to READscript.R
    ln -sf ${READ_HOME}/READscript.R
    touch meansP0_AncientDNA_normalized
}

# -------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- DOWNLOAD DEPENDENCIES ----------------------------------------------- #
header "SETTING DEPENDENCIES" 


# ---- Download READ
if [ ! -f "${READ_HOME}/READ.py" ]; then
    log INFO "Installing read in ${READ_HOME}"  
    mkdir -p ${READ_HOME} && pushd ${READ_HOME}
    READ_REPO="https://bitbucket.org/tguenther/read.git"
    git clone ${READ_REPO} .
    popd
else
    log INFO "Using READ in ${READ_HOME}"
fi

# ---- Download PLINK1.9
if [ ! -d ${PLINK_HOME} ]; then 
    log INFO "Installing Plink1.9 in ${PLINK_HOME}"
    mkdir -p ${PLINK_HOME} && pushd ${PLINK_HOME} > /dev/null
    # NOTE :This is build is intended for X64_linux architectures.
    #       If you have a different architecture, go to https://www.cog-genomics.org/plink/1.9/
    #       and change the URL of \$PLINK_ZIP to one that is more appropriate
    log WARN "Assuming we're on X86_X64_linux architecture. Change the \$PLINK_ZIP variable value at (Line ${LINENO} to something more appropriate if this is not the case"
    PLINK_ZIP="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210606.zip"
    wget -qO tmp.zip $PLINK_ZIP && unzip -x tmp.zip && rm tmp.zip
    popd
else
    log INFO "Using Plink1.9 in ${PLINK_HOME}"
fi

# Download Eigensoft
if [ ! -x "${EIGENSOFT_HOME}/bin/convertf" ]; then
    # Check for and ask to install dependencies.
    log INFO "Installing eigensoft in ${EIGENSOFT_HOME}"
    log WARN "'sudo apt-get install libgsl-dev libopenblas-dev liblapack-dev liblapacke-dev' may be required "
    EIGENSOFT_REPO="https://github.com/DReichLab/EIG.git"
    mkdir -p ${EIGENSOFT_HOME} && pushd ${EIGENSOFT_HOME}
    git clone $EIGENSOFT_REPO .
   
    log INFO "Using ${PROCS} threads for compilation."
    sleep 2
    export LDLIBS="-llapacke" 
    cd ./src && make convertf --jobs ${PROCS}  && make install && cd - || abort "
    Failed to compile Eigensoft. Most probably related to this issue: https://github.com/argriffing/eigensoft/issues/2 . 
    Try running:\n\n
    awk -i inplace 'match(\$0,/ -lm/){print substr(\$0,1,RSTART),substr(\$0,RSTART+RLENGTH+1),substr(\$0,RSTART+1,RLENGTH-1); next}{print}' ${EIGENSOFT_HOME}/src/Makefile\n\n
    And rerun this script.
    "
    popd
else
    log INFO "Using convertf in ${EIGENSOFT_HOME}/bin/convertf"
fi
    

# -------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------- DOWNLOAD REICH DATASET ------------------------------------------------- #
header "DOWNLOADING REICH DATASET" 

# Download Reich Merge dataset.
if [ ! -f "${REICH_DATA}.geno" ]; then
    REICH_DIR=$(dirname $REICH_DATA)
    log INFO "Downloading Reich Merge dataset in ${REICH_DIR}"
    mkdir -p $REICH_DIR && pushd $REICH_DIR
    URL="https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V44/V44.3/SHARE/public.dir/v44.3_1240K_public.tar"
    wget $URL
    tar -xvf v44.3_1240K_public.tar || abort "Failed to untar ${REICH_DATA}.tar"  
    popd
fi

# Convert to Reich-merge to ped/map
if [ ! -f "${REICH_DATA}.ped" ]; then
    log INFO "Converting Reich Merge dataset from EIGENSTRAT to PED"
    REICH_DIR=$(dirname $REICH_DATA)
    REICH_FILE=$(basename $REICH_DATA) 
    CONVERTF_FILE="par.EIGENSTRAT.PED"

    pushd ${REICH_DIR}
    echo "" > $CONVERTF_FILE 
    echo "genotypename:           ${REICH_FILE}.geno"   >> ${CONVERTF_FILE} 
    echo "snpname:                ${REICH_FILE}.snp"    >> ${CONVERTF_FILE}
    echo "indivname:              ${REICH_FILE}.ind"    >> ${CONVERTF_FILE}
    echo "outputformat:           PED"                  >> ${CONVERTF_FILE}
    echo "genotypeoutname:        ${REICH_FILE}.ped"    >> ${CONVERTF_FILE}
    echo "snpoutname:             ${REICH_FILE}.map"    >> ${CONVERTF_FILE}
    echo "indivoutname:           ${REICH_FILE}.pedind" >> ${CONVERTF_FILE}

    cat ${CONVERTF_FILE}

    ${EIGENSOFT_HOME}/bin/convertf -p ${CONVERTF_FILE} || abort "Failed EIGENSTRAT -> PED conversion"
    popd
fi


# -------------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------- PREPROCESS FILES ---------------------------------------------------- #
header "PREPROCESSING DATA"

FILE_HEADER=$(basename ${MT_DATA})

# Merge MT data 
MT_ONLY_DATA="${PREPROCESS_DIR}/${FILE_HEADER}.1240K_maf_MT_only"
if [ ! -f "${MT_ONLY_DATA}.bed" ]; then
    log INFO "STEP I - Extracting MT_individuals from ${MT_DATA}"
    keep_file="./in/keep_MT_only.txt"
    make_keep_file ${MT_DATA} "MT[0-9]{1,2}" "${keep_file}"
    log INFO "Running Plink..."
    ${PLINK_HOME}/plink --threads ${PROCS} --file ${MT_DATA} --keep "${keep_file}" --out ${MT_ONLY_DATA} --make-bed --allow-no-sex || abort
else
    log WARN "Skipping STEP I - Found ${MT_ONLY_DATA} dataset"
fi

MERGED_REICH_MT="${PREPROCESS_DIR}/${FILE_HEADER}.v44.3_1240k_public_MTmerged"
if [ ! -f "${MERGED_REICH_MT}.bed" ]; then
    log INFO "STEP II - Merging MT individuals with Reich dataset"
    ${PLINK_HOME}/plink --threads ${PROCS} --file ${REICH_DATA} --bmerge ${MT_ONLY_DATA} --merge-mode 2 --make-bed --allow-no-sex --out ${MERGED_REICH_MT} || abort
else
    log WARN "Skipping STEP II - Found ${MERGED_REICH_MT} dataset"
fi
########

# Keep only autosomal positions and maf > 0.05
MT_MAF_DATA="${PREPROCESS_DIR}/${FILE_HEADER}.v44.3_1240k_public_MT_merged_maf_${MAF#0.}"
if [ ! -f "${MT_MAF_DATA}.bed" ]; then
    log INFO "STEP III - Filter SNP with MAF > ${MAF} from chr1-22"
    ${PLINK_HOME}/plink  --threads ${PROCS} --bfile ${MERGED_REICH_MT} --maf ${MAF} --chr 1-22 --make-bed --out "${MT_MAF_DATA}" || abort
else
    log WARN "Skipping STEP III - Found ${MT_MAF_DATA} dataset"
fi

ART_MT_DATA="${PREPROCESS_DIR}/${FILE_HEADER}.v44.3_1240k_public_MT_merged_maf_${MAF#0.}_${REGEX}_MT"
if [ ! -f "${ART_MT_DATA}.bed" ]; then
    log INFO "STEP IV - Extract ${REGEX} and MT individuals."
    keep_file="./in/v44.3_keep_inds.txt"
    make_keep_file ${MT_MAF_DATA} "M[T]{1,2}[0-9]|${REGEX}" "${keep_file}"
    ${PLINK_HOME}/plink --threads ${PROCS} --bfile ${MT_MAF_DATA} --keep ${keep_file} --make-bed --chr 1-22 --out "${ART_MT_DATA}" || abort
else
    log WARN " Skipping STEP IV - Found ${ART_MT_DATA}"
fi

ART_ONLY_DATA="out/READ/${REGEX}/${FILE_HEADER}.v44.3_1240k_public_MT_merged_maf_${MAF#0.}_${REGEX}_only"
if [ ! -f "${ART_ONLY_DATA}.tped" ]; then
    log INFO "STEP V - Extract $REGEX{} and transpose to tped/tfam"
    keep_file="./in/v44.3_keep_${REGEX}_only.txt"
    make_keep_file ${ART_MT_DATA} "${REGEX}" "${keep_file}"
    mkdir -p "./out/READ/${REGEX}"
    ${PLINK_HOME}/plink --threads ${PROCS} -bfile ${ART_MT_DATA} --keep ${keep_file} --recode transpose --out ${ART_ONLY_DATA} || abort
else
    log WARN "Skipping STEP V - Found ${ART_ONLY_DATA} dataset"
fi 

MT_ONLY_DATA="out/READ/MT/${FILE_HEADER}.v44.3_1240k_public_MT_merged_maf_${MAF#0.}_MT_only"
if [ ! -f "${MT_ONLY_DATA}.tped" ]; then
    log INFO "STEP VI - Extract MT and transpose to tped/tfam"
    keep_file="./in/v44.3_keep_MT_only.txt"
    make_keep_file ${ART_MT_DATA} "MT" "${keep_file}"
    mkdir -p "./out/READ/MT"
    ${PLINK_HOME}/plink --threads ${PROCS} -bfile ${ART_MT_DATA} --keep ${keep_file} --recode transpose --out ${MT_ONLY_DATA} || abort
else
    log WARN "Skipping STEP VI - Found ${MT_ONLY_DATA}"
fi 


# -------------------------------------------------------------------------------------------------------------------- #
# ------------------------------------------------ Running READ  ----------------------------------------------------- #
header "PERFORMING KINSHIP ANALYSIS USING READ"


log INFO "STEP VII -  Run READ on ART"
pushd "$(dirname ${ART_ONLY_DATA})"
prime_READ;
/usr/bin/python2 ${READ_HOME}/READ.py $(basename ${ART_ONLY_DATA}) median - || abort
popd
 
# ---- GET MEDIAN P0 value from our proxy individuals.
MEDIAN_P0=$(Rscript ./src/get_median_P0.R "$(dirname ${ART_ONLY_DATA})/meansP0_AncientDNA_normalized")
log WARN "Median P0 for ${REGEX}: ${MEDIAN_P0}"

log INFO "STEP VIII - Run READ on MT using ${MEDIAN_P0} as normalization value"

pushd "$(dirname ${MT_ONLY_DATA})"
prime_READ;
${PYTHON} ${READ_HOME}/READ.py $(basename ${MT_ONLY_DATA}) value "${MEDIAN_P0}" || abort
popd 


# --- Merge READ results.
log INFO "STEP IX - Merge ${ART_ONLY_DATA} and ${MT_ONLY_DATA} results for plotting."
MERGED_P0="results/MT_${REGEX}/meansP0_AncientDNA_normalized"
mkdir -p $(dirname ${MERGED_P0})
touch "${MERGED_P0}"
echo "PairIndividuals Normalized2AlleleDifference StandardError NonNormalizedP0 NonNormalizedStandardError" > ${MERGED_P0}
tail -q -n+2 `find "./out/READ/" -name "$(basename ${MERGED_P0})" | grep -E "${REGEX}|MT" | xargs` >> ${MERGED_P0}



MERGED_RESULTS="results/MT_${REGEX}/READ_results"
mkdir -p $(dirname ${MERGED_RESULTS})
touch "${MERGED_RESULTS}"
echo -e "PairIndividuals\tRelationship\tZ_upper\tZ_lower" > ${MERGED_RESULTS}
tail -q -n+2 `find "./out/READ" -name "$(basename ${MERGED_RESULTS})" | grep -E "${REGEX}|MT" | xargs` >> ${MERGED_RESULTS}


header "\e[32mAll done!\e[37m"

