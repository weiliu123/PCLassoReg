#' get protein complexes
#'
#' @param Groups A data frame containing the protein complexes
#' @param Organism Organism. one of \code{Human}, \code{Mouse}, \code{Rat},
#' \code{Mammalia}, \code{Bovine}, \code{Dog}, or \code{Rabbit}.
#' @param Type The name type of the proteins in the protein complexes. One of
#' \code{GeneSymbol}, \code{EntrezID}, \code{UniprotID}
#'
#' @return A list of protein complexes
#' @export
#'
#' @examples
#' data(PCGroups)
#' PC.Human <- getPCGroups(Groups = PCGroups, Organism = "Human",
#' Type = "GeneSymbol")
getPCGroups <- function(Groups,
                        Organism = c("Human", "Mouse", "Rat", "Mammalia",
                                     "Bovine", "Dog", "Rabbit"),
                        Type = c("GeneSymbol", "EntrezID", "UniprotID")){
    Organism = match.arg(Organism)
    Type = match.arg(Type)

    proCom <- Groups[which(Groups$Organism == Organism),]

    # protein complex list
    Complexs <- vector(mode = "list", length = nrow(proCom))
    names(Complexs) <- paste0("Complex.", proCom$ComplexID)

    # protein set in all protein complexes
    for(i in 1:nrow(proCom)){
        prots.str <- as.character(proCom[i, Type])
        prots.i <- unlist(strsplit(prots.str, split = ";"))
        prots.i <- prots.i[which(prots.i != "None")]
        Complexs[[i]] <- prots.i
    }

    Complexs
}
