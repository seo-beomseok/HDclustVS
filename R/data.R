#' A single-cell RNA-sequencing data from primary mamalian cells
#'
#' A dataset consisting of 131 primary mammalian cells with 534 genes. This data set contains only two cell clusters: unstimulated and stimulated cells.
#'
#' @format A list containing four components
#'   \item{cell}{cell names}
#'   \item{cell.stage}{true labels of cells}
#'   \item{gnames}{gene names}
#'   \item{expr}{expression data}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/19729616}
"LPS"

#' A single-cell RNA-sequencing data from mouse cells.
#'
#' A dataset consisting of 92 single cells with 35 genes from two different mouse cell types: embryonic stem cells (ES R1) and embryonic fibroblasts (MEFs).
#'
#' @format A list containing four components
#'   \item{cell}{cell names}
#'   \item{cell.stage}{true labels of cells}
#'   \item{gnames}{gene names}
#'   \item{expr}{expression data}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/19729616}
"MEF"

#' A single-cell RNA-sequencing data from human preimplantation embryos and human embryonic stem cells.
#'
#' A dataset consisting of 124 individual cells with 3,840 genes from human preimplantation embryos and human embryonic stem cells. There are 7 cell clusters: Oocyte and Zygote, 2-cell, 4-cell, 8-cell, Morula, Late blastocyst, hESC.
#'
#' @format A list containing four components
#'   \item{cell}{cell names}
#'   \item{cell.stage}{true labels of cells}
#'   \item{gnames}{gene names}
#'   \item{expr}{expression data}
#' @source \url{https://www.ncbi.nlm.nih.gov/pubmed/23934149}
"YAN"

#' Simulated toy data
#'
#' A dataset containing 300 samples and 200 variables.
#'
#' @format A list containing four components
#'   \item{X_total}{total data}
#'   \item{X_relev}{relevant data}
#'   \item{X_irrel}{irrelevant data}
#'   \item{z}{true labels of clusters}
"sim1"


