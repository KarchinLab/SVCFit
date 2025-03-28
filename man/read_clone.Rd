% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_clone.R
\name{read_clone}
\alias{read_clone}
\title{This function reads the ground truth designed for a simulation}
\usage{
read_clone(truth_path, mode = "inherited")
}
\arguments{
\item{truth_path}{an onject of class "Character". This variable is a path to
bed files storing true structural variants information with clonal assignment.
Each bed file should be named as "clone" + "number". Structural variants should
be saved in a seperate bed file if they belong to different (sub)clone.}

\item{mode}{an onject of class "Character". This variable describe how true
clonal information is saved. In "inherited" mode, bed files for all children
clone contains all ancestral structural variants of their parents. In "distinct"
mode, children clones don't contain any ancestral structural variants.}
}
\value{
a dataframe
}
\description{
This function reads the ground truth designed for a simulation
}
