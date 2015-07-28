\name{Environment}
\alias{ENV}
\alias{.Random.seed}
\title{
Trait by Environment
}
\description{
Matrix of phenotypes of a given trait for each environment.
}
\usage{
ENV(trait="yield")
}
\arguments{
  \item{trait}{
Character. Trait of interest. The options are: "yield" (grain yield in Kg/ha), "maturity" (days to maturity), "height" (plant height in cm), "lodging" (lodging score from 1 to 5), "protein" (protein percentage in the grain), "oil" (oil percentage in the grain), "size" (seed size = mass of 100 seeds in grams) and "fiber" (fiber percentage in the grain).
}
}
\value{
Numeric matrix with genotypes as rows and environments as columns.
}
\author{
Alencar Xavier
}
\examples{
O = ENV("oil")
O = O[,c(1,3,4,6,13,14)]
pairs(O)
}