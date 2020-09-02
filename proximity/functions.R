# Supporting functions for proximity score analysis

proximity <- function(graph, to, from) {
    # use igraph's shortest_paths function to get the shortest paths
    # just take the length of the first one as the proximity
    sp <- shortest_paths(g, from, to)
    l <- length(sp$vpath[[1]])
    # if there is no path b/w the two, return NaN
    if (l == 0) return(NaN) else return(l - 1)
}

proximity_matrix <- function(graph, gene_list_1, gene_list_2) {
    m <- sapply(gene_list_1, function(x) {
        sapply(gene_list_2, function(y) {
            proximity(graph, x, y)
        })
    })
}

filter_gene_lists <- function(x, universe) {
    tmp <- sapply(x, function(y) {
        intersect(y, universe)
    }, USE.NAMES = T)
    # remove any that now have zero genes
    return(tmp[lengths(tmp) > 0])
}