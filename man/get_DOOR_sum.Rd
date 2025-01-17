% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DOOR_sum_v0.R
\name{get_DOOR_sum}
\alias{get_DOOR_sum}
\title{Summarize DOOR variable by treatment groups}
\usage{
get_DOOR_sum(data, DOORvar, trtvar, DOORlevel, trt = "1", ctrl = "0")
}
\arguments{
\item{data}{A dataframe contains DOOR variable and Treatment variable}

\item{DOORvar}{Variable name of DOOR}

\item{trtvar}{Variable name of treatment}

\item{DOORlevel}{number of DOOR level pre-specified}

\item{trt}{Char label for treatment indicator}

\item{ctrl}{char label for control indicator}
}
\description{
Summarize DOOR variable by treatment groups
}
