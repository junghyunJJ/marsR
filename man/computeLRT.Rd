\name{computeLRT}
\alias{computeLRT}
\title{Likelihood Ratio Test (\emph{LRTscore})}
\description{
To test the association between a set of variants of a locus and a trait, 
MARS performs the Likelihood Ratio Test (LRT) 
to compute a test statistic, referred to as \emph{LRTscore}.
}
\usage{
computeLRT(stat, geno,
           subsize = 50, 
           causalCount = 2,
           NCP = 5.7, 
           gamma = 0.01)
}
\arguments{
  \item{stat}{a m x 1 numeric matrix}
  \item{geno}{a m x n numeric matrix}
  \item{subsize}{subsize of SNPs (default: 50)}
  \item{causalCount}{set the maximum number of causal SNPs (default: 2)}
  \item{NCP}{non-centrality parameter matrix (default: 5.7)}
  \item{gamma}{set the maximum number of causal SNPs (default: 0.01)}
}
\details{
MARS consider likelihoods of two models; the likelihood of the null model (L0) and the likelihood of the alternative model (L1). The null model assumes that there is no causal variant to the trait and the alternative model assumes that there is at least one causal variant to the trait. Then we compute the \emph{LRTscore} as L1/L0.
}
\references{
  Leveraging allelic heterogeneity to increase power of association testing
(https://www.biorxiv.org/content/10.1101/498360v1)
}
\examples{
data(testdata)
computeLRT(testdata$stat,testdata$geno)
}
\author{Junghyun Jung}


