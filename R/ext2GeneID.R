ext2GeneID <- function(x, split = "_"){
    x.GeneSymbol <- rep(0, length = length(x))
    for(i in 1:length(x)){
        x.GeneSymbol[i] <- unlist(strsplit(x[i], split = split))[2]
    }
    unique(x.GeneSymbol)  # 去除重复
}
