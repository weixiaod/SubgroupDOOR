% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_based_CIs.R
\name{getCI_clm}
\alias{getCI_clm}
\title{Calculate the pointwise confidence interval for DOOR probabilities for new dataset}
\usage{
getCI_clm(fit.clm, newdata, DOORlevel)
}
\arguments{
\item{fit.clm}{the model fit}

\item{newdata}{the new dataset contains arm indicator and covariate}

\item{DOORlevel}{Total DOOR level K}
}
\value{
a dataframe contains the covariate, DOOR probability, and confidence interval
}
\description{
Calculate the pointwise confidence interval for DOOR probabilities for new dataset
}
