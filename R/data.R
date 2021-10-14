#' Protein complexes.
#'
#' A dataset containing the protein complexes
#'
#' @format A data frame with 3512 rows and 6 variables:
#' \describe{
#'   \item{ComplexID}{ID of the protein complex}
#'   \item{ComplexName}{name of the protein complex}
#'   \item{Organism}{organism}
#'   \item{UniprotID}{Uniprot IDs of the proteins in the protein complex}
#'   \item{EntrezID}{Entrez IDs of the proteins in the protein complex}
#'   \item{GeneSymbol}{gene symbols of the proteins in the protein complex}
#' }
#' @source \url{https://mips.helmholtz-muenchen.de/corum/}
"PCGroups"


#' A dataset for classification
#'
#' @format A list containing a protein expression matrix and a response vector
#' \describe{
#'   \item{Exp}{a protein expression matrix}
#'   \item{Label}{a response vector}
#' }
"classData"

#' A dataset for prognostic model
#'
#' @format A list containing a protein expression matrix and survival data
#' \describe{
#'   \item{Exp}{a protein expression matrix}
#'   \item{survData}{Survival data. The first column is the time on study
#'   (follow up time); the second column is a binary variable with
#'   1 indicating that the event has occurred and 0 indicating right censoring.}
#' }
"survivalData"
