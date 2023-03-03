# 20210723 - Mentesh Tepe EAGER preprocessing.

## Summary:
- **Date of analysis:** 2021-07-23  
- **Aim:** Replicate the data preprocessing pipeline found within (Skourtanioti et al. 2020)'s paper, since we plan to compare MT23, MT26 and MT7 with  with MTT001 and ArslanTepe / Tell-Kurdu individuals.
  - Latest release of EAGER-v1 is used. (Note that this software has not been maintained since Jan. 2018.)
  - READ analysis is then performed using a two-pass methods see [README-READ-analysis.md](README-READ-analysis.md):
    - 1. Perform Kinship analysis on a surrogate population (ART / KRD), and extract the median of meansP0
    - 2. Perform Kinship analysis on Mentesh Tepe individuals using the previously obtained median(meansP0)

---

## Pipeline description
### Adapters:

 - Adapter1: `GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG`
 - Adapter2: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT`

### EAGER preprocessing.
1. Sequencing adapters are clipped with `AdapterRemoval v2.2.0` and collapsed. fragments shorter than 30bp are discarded.
2. Mapping is performed with `bwa v0.7.12`", using the `hs37d5` + applying a quality filtering of `q30`.
3. PCR duplicates are removed using `dedup v0.12.` (this software considers the whole fragment composition to estimate potential PCR duplicates, whereas Picard will typically only match the beginning of the sequence, and leverage the sequence length to estimate duplication)
4. PMD misincorporations is evaluated with `MapDamage v2.0.6`

### Post-EAGER preprocessing
5. "Bam files are merged across libraries for the same individual"
6. "PCR duplicate removal is performed on the merged files a second time, using `dedup v0.12`.

### ***[LEGACY]*** Selective masking and interlacing of PMD misincorporations.
**Important notice:** This section was ***not*** included in our final data processing protocol, but was merely an attempt to replicate the approaches found in (Skourtanioti et al. 2020). This approach does not fair well with non UDG-treated data, and the more sensible approach for our case was to simply apply 'standard' PMD-rescaling using MapDamage. 

This section, and the corresponding scripts are thus only kept for archiving purposes. The rest of the analysis can be found in [README-READ-analysis.md](README-READ-analysis.md)

7. Using `trimBam` generated soft-clipped versions of the bam files, in which nucleotides are masked until the estimated nucleotide misincorporation frequency is `<1%`

8. "Generate pileup files with `-q 30` and `-Q 30` filters (on both masked and non-masked versions)."
9. "Perform random pseudo-haploid variant calling `pileupCaller` (on both masked and non-masked versions)"

10. "Interlace the genotype calling files, keeping only the *transitions* from the masked version, and the *transversions* from the non-masked version."

## Resources:
#### Skourtianoti:

 - **Citation:**
   ```Text
   Skourtanioti E, Erdal YS, Frangipane M, Balossi Restelli F, Yener KA, Pinnock F, Matthiae P, Özbal R, Schoop UD, Guliyev F, Akhundov T, Lyonnet B, Hammer EL, Nugent SE, Burri M, Neumann GU, Penske S, Ingman T, Akar M, Shafiq R, Palumbi G, Eisenmann S, D'Andrea M, Rohrlach AB, Warinner C, Jeong C, Stockhammer PW, Haak W, Krause J. Genomic History of Neolithic to Bronze Age Anatolia, Northern Levant, and Southern Caucasus. Cell. 2020 May 28;181(5):1158-1175.e28. doi: 10.1016/j.cell.2020.04.044. PMID: 32470401.
   ```
 - **DOI:** https://doi.org/10.1016/j.cell.2020.04.044

#### EAGER:
 - Source code: https://github.com/apeltzer/EAGER-GUI
 - Documentation: https://eager.readthedocs.io/en/latest/
 - Publication:
   ```Text
   A. Peltzer; G. Jäger; A. Herbig; S. Seitz; C. Kniep; J. Krause; K. Nieselt: EAGER: efficient ancient genome reconstruction (Genome Biology 2016, 17:60, doi:10.1186/s13059-016-0918-z)
   ```
#### READ:
 - Source code: https://bitbucket.org/tguenther/read/src/master/
 - Documentation: N/A (Github's README is all there is)
 - Publication: 
   ```Text  
   Monroy Kuhn JM, Jakobsson M, Günther T. Estimating genetic kin relationships in prehistoric populations. PLoS One. 2018 Apr 23;13(4):e0195491. doi: 10.1371/journal.pone.0195491. PMID: 29684051; PMCID: PMC5912749.
   ```