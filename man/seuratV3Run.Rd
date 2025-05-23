% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seuratV3Run.R
\name{seuratV3Run}
\alias{seuratV3Run}
\title{Seurat V3 method}
\usage{
seuratV3Run(
  input,
  assay = NULL,
  reference = NULL,
  anchor.features = 2000,
  scale = TRUE,
  normalization.method = "LogNormalize",
  sct.clip.range = NULL,
  reduction = "cca",
  l2.norm = TRUE,
  dims = 1:30,
  k.anchor = 5,
  k.filter = 200,
  k.score = 30,
  max.features = 200,
  nn.method = "annoy",
  n.trees = 50,
  eps = 0,
  verbose = TRUE,
  new.assay.name = "integrated",
  features = NULL,
  features.to.integrate = NULL,
  k.weight = 100,
  weight.reduction = NULL,
  sd.weight = 1,
  sample.tree = NULL,
  preserve.order = FALSE
)
}
\arguments{
\item{input}{A list of Seurat objects.}

\item{assay}{A vector of assay names specifying which assay to use when constructing anchors. If NULL, the current default assay for each object is used.}

\item{reference}{A vector specifying the object/s to be used as a reference during integration.}

\item{anchor.features}{Number of features to be used in anchor finding.}

\item{scale}{A logical}

\item{normalization.method}{Name of normalization method used: LogNormalize (default) or SCT.}

\item{sct.clip.range}{sct.clip.range}

\item{reduction}{reduction}

\item{l2.norm}{l2.norm}

\item{dims}{dims}

\item{k.anchor}{k.anchor}

\item{k.filter}{k.filter}

\item{k.score}{k.score}

\item{max.features}{max.features}

\item{nn.method}{nn.method}

\item{n.trees}{n.trees}

\item{eps}{eps}

\item{verbose}{verbose}

\item{new.assay.name}{new.assay.name}

\item{features}{features}

\item{features.to.integrate}{features.to.integrate}

\item{k.weight}{k.weight}

\item{weight.reduction}{weight.reduction}

\item{sd.weight}{sd.weight}

\item{sample.tree}{sample.tree}

\item{preserve.order}{preserve.order}
}
\description{
Seurat V3 method
}
