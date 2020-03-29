\name{mars}
\alias{mars}
\title{Model-based Association test Reecting causal Status (MARS)}
\description{
To find an association between variants in risk loci 
and a phenotype, considering the causal status of the variants.}
\usage{
mars(stat, geno, 
     subsize = 50, 
     causalCount = 2, 
     NCP = 5.7, 
     gamma = 0.01, 
     simNum = 1000, 
     fast = TRUE, 
     threshold = 5e-6, 
     setseed = 1, 
     cores = 1)
}
\arguments{
  \item{stat}{a m x 1 numeric matrix}
  \item{geno}{a m x n numeric matrix}
  \item{subsize}{subsize of SNPs (default: 50)}
  \item{causalCount}{set the maximum number of causal SNPs (default: 2)}
  \item{NCP}{non-centrality parameter matrix (default: 5.7)}
  \item{gamma}{set the maximum number of causal SNPs (default: 0.01)}
  \item{simNum}{set the maximum number of causal SNPs (default: 10000)}
  \item{fast}{importance sampling (default: TRUE)}
  \item{threshold}{set the maximum number of causal SNPs (default: 5e-6)}
  \item{setseed}{seed for the random number generator (default: 1)}
  \item{cores}{he number of cores to use (default: 1)}
      
}
\details{
Over the past decade, genome-wide association studies (GWAS) have successfully identified many variants significantly associated with diseases and complex traits. Unfortunately, those variants explain extremely small proportion of phenotypic variation (Yang et al. Nature Genet. 2010; Wood et al. Nature Genet. 2014) and there are many more variants with even smaller effects that we are yet to identify (D.B. Goldstein NEJM 2009; Yang et al. Nature Genet. 2010; Loh et al. Nature Genet. 2015; Boyle et al. Cell 2017). Biological follow up studies are performed on loci that harbor at least one significant variant and detecting all loci that harbor associated risk loci can help elucidate the biological mechanisms of diseases. To detect association loci, most of the existing methods assume that each locus has only one causal variant. However, recently it has been reported that allelic heterogeneity (i.e., more than one causal variant exists in a locus) is widespread (Yang et al. Nature Genet. 2012, Flister et al. Genome Res. 2013).

}
\references{
  Leveraging allelic heterogeneity to increase power of association testing
(https://www.biorxiv.org/content/10.1101/498360v1)
}
\examples{
data(testdata)
mars(testdata$stat, testdata$geno, threshold = 5e-08)
}
\author{Junghyun Jung}


