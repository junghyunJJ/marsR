# Leveraging allelic heterogeneity to increase power of association testing

The standard genome-wide association studies (GWAS) detect an association between a single variant and a phenotype of 
interest. Recently, several studies reported that at many risk loci, there may exist multiple causal variants. For a locus with multiple 
causal variants with small effect sizes, the standard association test is underpowered to detect the associations. Alternatively, an 
approach considering effects of multiple variants simultaneously may increase statistical power by leveraging effects of 
multiple causal variants. In this paper, we propose a new statistical method, Model-based Association test Reflecting causal 
Status (MARS), that tries to find an association between variants 
in risk loci and a phenotype, considering the causal status of the variants. One of the main advantages of MARS is that it only requires the 
existing summary statistics to detect associated risk loci. Thus, MARS is applicable to any association study with summary statistics, even 
though individual level data is not available for the study. Utilizing extensive simulated data sets, we show that MARS increases the 
power of detecting true associated risk loci compared to previous approaches that consider multiple variants, while robustly controls the type I error. 
Applied to data of 44 tissues provided by the Genotype-Tissue Expression (GTEx) consortium, we show that MARS identifies more eGenes compared to previous approaches in most of the tissues; e.g. MARS identified 16\% more eGenes than the ones reported by the GTEx consortium. Moreover, applied to Northern Finland Birth Cohort (NFBC) data, we demonstrate that MARS effectively identifies association loci with improved power (56\% of more loci found by MARS) in GWAS studies compared to the standard association test.


## Installation
- The C++ library for GNU [GSL](https://www.gnu.org/software/gsl/) is required.
- [makesigmasemiPDRcppGSL](https://github.com/junghyunJJ/makesigmasemiPDRcppGSL) is required.
- marsR works only on ***nix (Linux, Unix such as macOS) system**. please check **.Platform$OS.type** function.
- We currently only support R 3.5+.*

```
install.packages("Rcpp")
install.packages("RcppGSL")
install.packages("RcppArmadillo")
install.packages("devtools")
devtools::install_github("junghyunJJ/makesigmasemiPDRcppGSL")
devtools::install_github("junghyunJJ/marsR")
```

## Example

```
library(marsR)
> data(testdata)
> mars(testdata$stat, testdata$geno, threshold = 5e-08)
[2020-04-07 09:08:40] - calculating MARS alt LRT
[2020-04-07 09:08:40] - make null data
[2020-04-07 09:08:40] - calculating MARS null LRT
  |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=10s  
Total analysis time: 10.97 secs  
$alt
       LRT       uni
1 15985.79 1.968e-07

$results
  LRT_pvalue univariate_pvalue threshold_pvalue threshold_UNI significance_LRT significance_UNI
1          0        0.00099551       0.00099551         5e-08             TRUE            FALSE
```