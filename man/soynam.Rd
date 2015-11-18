\name{Dataset}
\alias{soynam}
\alias{data.check}
\alias{data.line}
\alias{gen.raw}
\alias{soybase}
\alias{data.check.qa}
\alias{data.line.qa}
\alias{gen.qa}
\alias{.Random.seed}
\title{
Datasets
}
\description{
Genotypes and phenotypes from quality assured (soybase) or original (soynam) datasets. 
}
\usage{
data(soybase)
data(soynam)
}
\details{
Datasets of the SoyNAM project, original and quality assured (QA) versions. Data was downloaded on November 16th 2015 from \url{soynam.org}.

Genotypic matrices are named "gen.raw" and "gen.qa" for the raw and QA versions, respectively. In each dataset, phenotypes are allocated into two objects, one with the lines ("data.line") and one with checks and parents ("data.checks"). Information on data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set and environment), height (in centimeters), R8 (number of days to maturity), planting date (501 represents may 1), flowering (701 represents july 1), maturity (901 represents september 1), lodging (score from 1 to 5), yield (in Kg/ha), moisture, protein (percentage in the seed), oil (percentage in the seed), fiber (percentage in the seed), seed size (in grams of 100 seeds).
}

\author{
Alencar Xavier
}
\examples{
data(soybase)
data(soynam)
}