% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HomologySeeker.R
\name{HomoSeeker}
\alias{HomoSeeker}
\title{Seeker for highly varible homologous gene sets between species}
\usage{
HomoSeeker(
  RefSpec,
  QuySpec,
  RefSpec_mat,
  QuySpec_mat,
  homo_mat = NULL,
  RefSpec_gene = "Gene_sym",
  QuySpec_gene = "Gene_sym",
  HVGs_method = "seurat_vst",
  Integrated_mat = FALSE,
  Integrated_method = NULL,
  verbose = TRUE
)
}
\arguments{
\item{RefSpec, QuySpec}{Names of reference(ref) and query(quy) species in comparative analysis. Case is ignored.
See GetSpecNames() for detailed list.}

\item{RefSpec_mat, QuySpec_mat}{Single cell expression matrix of reference or query species with row as gene
symbol/ID and column as sample ID.}

\item{homo_mat}{Table of one to one homologous genes table between reference and query species. Default to use homologous
gene table collected by HomoSelector(). To use self constructed homologous gene table, please input table with:
#' \itemize{
 \item{1st Column} \strong{:} Gene symbel/ID of reference species.
 \item{2nd Column} \strong{:} Corresponding Gene symbel/ID of query species.
 }}

\item{RefSpec_gene, QuySpec_gene}{Type of gene name of Single cell expression matrix.
\itemize{
 \item{Gene_sym} \strong{:} Gene symbel (default).
 \item{Gene_id} \strong{:} Gene ID.
 }}

\item{HVGs_method}{HVG selection method to be used. Supported methods:
\itemize{
 \item{seurat_vst(default)}
 \item{seurat_sct}
 \item{scran}
 \item{scmap}
 \item{ROGUE}
 }
 See HVGSelector() for detailed description.}

\item{Integrated_mat}{Whether input matrics are returned by integration}

\item{Integrated_method}{If Integrated_mat=TRUE, please specify what integration method was used.
currently, only seurat_LogNorm and seurat_SCT intergation are supported}

\item{verbose}{Whether show calculation progress. Default is TRUE.}
}
\value{
Returns a HomologyHVG object with slot:
\itemize{
 \item{Species} \strong{:} Query species names.
 \item{HomologyHVG} \strong{:} Highly variable gene sets for species 1 and 2.
 \item{Matrix_orig} \strong{:} Original single cell matrices for species 1 and 2.
 \item{Matrix_homo} \strong{:} Single cell matrices with one-to-one homologous genes between species 1 and 2 as row.
 \item{HVG_feature} \strong{:} Returned HVG information of all homologous genes for species 1 and 2.
 \item{Table_homo} \strong{:} Returned homologous genes table between species 1 and 2.
 }
}
\description{
Get homologous information and variable gene sets between two species
}
