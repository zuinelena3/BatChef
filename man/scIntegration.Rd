% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generic.R, R/scIntegration.R
\name{scIntegration}
\alias{scIntegration}
\alias{scIntegration,limmaMethod-method}
\alias{scIntegration,combatMethod-method}
\alias{scIntegration,seuratv3Method-method}
\alias{scIntegration,seuratv5Method-method}
\alias{scIntegration,fastMNNMethod-method}
\alias{scIntegration,harmonyMethod-method}
\alias{scIntegration,scanoramaMethod-method}
\alias{scIntegration,bbknnMethod-method}
\alias{scIntegration,scVIMethod-method}
\alias{scIntegration,scMergeMethod-method}
\title{scIntegration}
\usage{
scIntegration(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{limmaMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{combatMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{seuratv3Method}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{seuratv5Method}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{fastMNNMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{harmonyMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{scanoramaMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{bbknnMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{scVIMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)

\S4method{scIntegration}{scMergeMethod}(
  obj,
  batch = NULL,
  assay = NULL,
  hvgs = NULL,
  dims = NULL,
  reduction = NULL,
  anchor = NULL,
  k_anchor = NULL,
  genelist = NULL,
  cell_type = NULL,
  METHOD,
  alt_out = FALSE
)
}
\arguments{
\item{obj}{A SingleCellExperiment}

\item{batch}{a string specifying the batch. batch = NULL when obj is a list}

\item{assay}{a string specifying the assay to use for correction}

\item{hvgs}{number of highly variable genes}

\item{dims}{number of dimensions to use for dimension reduction}

\item{reduction}{a string specifying the dimension reduction to use for correction}

\item{anchor}{a string specifying the anchors finding type (CCAIntegration, RPCAIntegration, HarmonyIntegration, JointPCAIntegration)}

\item{k_anchor}{number of anchors (default: k_anchor = 5)}

\item{genelist}{negative controls}

\item{cell_type}{string specifying the cell-type labels}

\item{METHOD}{a \code{MethodParam} object specifying the batch correction method}

\item{alt_out}{alternative output: a \code{\link{AltOutput}} class}
}
\description{
Generic function for defining batch effect correction methods.

A common interface for single-cell batch correction methods.
}
