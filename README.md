#### Summary
**Date of analysis:** 2022-03-15
**Aim:** Re-run READ Kinship analysis with MT-EAGER bam files WITHOUT running manual PMD damage trimming and interlacing, as in ([Skourtanioti,2019](https://doi.org/10.1016/j.cell.2020.04.044)).

#### Data 
###### Reference genome: 
 - hs37d5
 - From: `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence`

###### TARGET Positions and ancient genomes:
 - D.Reich '1240K Dataset', A.K.A allen-present-day-and-ancient-genotypes.
 - From: `https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V44/V44.3/SHARE/public.dir/v44.3_1240K_public.snp`


#### Instructions:
###### A. Dependencies
1. Install miniconda3. See: [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
2. Create conda environments in ./envs
   ```Bash
   user@desktop:~$ for env in ./envs; do conda env create -f $env; done
   ```

###### B. Prepare environment
1. Create symlink to the indexed hs37d5 ref. genome and 1240K target position (see [Data](#Data))
   ```Bash
   user@desktop:~$ mkdir -p ./data
   user@desktop:~$ ln -st ./data ../../MT_EAGER/data/reference/
   user@desktop:~$ ln -st ./data ../../MT_EAGER/data/target_positions/
   ```

   

2. Create symlink to EAGER output.
   ```Bash
   user@desktop:~$ mkdir -p ./out && ln -st ./out/ ../../MT_EAGER/out/2-dedup/
   ```
   
   - **Note:** These `.bam` files should be merged across runs + re-deduplicated

###### C. Run `00-call-genotypes`
Calling `./src/00-call-genotypes` from the root directory will:
  1. Perform PMD rescaling using mapDamage
  2. Perform mpileup `-RB -q30 -Q30`
  3. Perform Variant calling using pileupCaller `--randomHaploid --sampleNames`, targeting 'v44.3_1240K_public.snp' dataset.
  4. Convert the pileupCaller output to PLINK PED/MAP/PEDIND format.

**Expected output:** `out/5-pileupCaller/AllMT_pileupCaller_RBq30Q30.(ped/map/pedind)`


###### D. Run `01-runREAD-1240K.sh`
Calling `./src/01-runREAD-1240K.sh` from the root directory will:
  1. Install Software dependencies in `./methods` (plink-1.9, READ, eigensoft)
  2. Download `./v44.3_1240K_public.(snp|geno|ind)` dataset in `./data/Reich-dataset/`
  3. Convert Reich dataset to PED/MAP/PEDIND format in `./data/Reich-dataset/`
  4. Extract Mentesh-Tepe Individuals from `The output of `00-call-genotypes`
  5. Merge these individuals with the Reich-dataset.
  6. Filter out to only keep autosomal positions and MAF>0.05 from this merged dataset
  7. Extract surrogate population individuals from this dataset (ART or KRD) and run READ on them.  
     - This will extract an medianP0 which we'll use as a reference for our Mentesh Tepe Individuals
  8. Extract Mentesh Tepe individuals from the merged dataset and run READ on them, using the previously calculated 'medianP0' from step 7.
  9. Merge the READ_results and meansP0_AncientDNA_normalized output of Step 7 and Step 8, for plotting.

**Expected output:** `out/READ/MT_ART/READ_results`  and `out/READ/MT_ART/meansP0_AncientDNA_normalized`

###### E. Plot the results

```Bash
user@desktop:~$ pushd results/MT_ART/
user@desktop:~$ ../../src/plot_READ_results.R --regex "MT.*MT.*" --mainName "Mentesh Tepe" --proxyName "Arslan Tepe"
```

**Expected output:** `READ_plot.html` and `READ_lib_plot`

