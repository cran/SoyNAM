\name{Dataset}
\alias{soynam}
\alias{data.check}
\alias{data.line}
\alias{gen.check}
\alias{gen.line}
\alias{gen.IA3023}
\alias{soybase}
\alias{data.check.qa}
\alias{data.line.qa}
\alias{gen.check.qa}
\alias{gen.line.qa}
\alias{.Random.seed}
\title{
Dataset
}
\description{
Genotypes and phenotypes from quality assured (soybase) or original (soynam) datasets. 
}
\usage{
data(soybase)
data(soynam)
}
\details{
Datasets of the SoyNAM project, original and quality assured (QA) versions.

Phenotypes comprise two datasets, one with the lines ("data.line") and one with checks and parents ("data.checks"). Information on data objects include year, location, environment (combination of year and location), strain, NAMfamily, set (set in each environment), spot (combination of set and environment), height (in centimeters), R8 (number of days to maturity), planting date (501 represents may 1), flowering (701 represents july 1), maturity (901 represents september 1), lodging (score from 1 to 5), yield (in Kg/ha), moisture, protein (percentage in the seed), oil (percentage in the seed), fiber (percentage in the seed), seed size (in grams of 100 seeds).

Non-segregating markers were removed, then genotypes were imputed using a hidden Markov model implemented in fastPHASE (Scheet and Stephens 2006). Alleles with minor allele frequency lower than 0.1 were removed. The genotypic dataset comprise three objects: lines ("gen.line") and parents ("gen.check") having IA3023 as reference to the allele coding, where 2 means homozygous towards IA3032 and 0 means homozygous to the founder parent. The third object is the genotype of IA3023 having Williams82 as the reference, such that original allele coding can be recovered using the function 'reference' from the NAM package.

QA version was designed by William Beavis (wdbeavis@iastate.edu), James Specht, Brian Diers, Reka Howard and Alencar Xavier. The QA genotypic and phenotypic datasets were downloaded from \url{soybase.org} on September 18, 2015.

Some caveats were found on the raw dataset. Family NAM02 displays 38 repeated lines that were found to be selfs from the the parent TN05. Other families may also display redundant lines. The BLUP function implemented in this package control this problems with repeated lines. A second possible issue regards the identity of the founder parent of Family NAM46.
}
\references{
Scheet, P., & Stephens, M. (2006). A fast and flexible statistical model for large-scale population genotype data: applications to inferring missing genotypes and haplotypic phase. The American Journal of Human Genetics, 78(4), 629-644.
}
\author{
Alencar Xavier
}
\examples{
data(soybase)
data(soynam)
}