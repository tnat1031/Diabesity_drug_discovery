# Supporting functions for proximity score analysis

proximity <- function(graph, to, from) {
    # use igraph's shortest_paths function to get the shortest paths
    # just take the length of the first one as the proximity
    sp <- shortest_paths(graph, from, to)
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

# helper function to summarize a proximity matrix
summarize_prox_mat <- function(prox_mat) {
    if (is.matrix(prox_mat)) {
        nr <- nrow(prox_mat)
        nc <- ncol(prox_mat)
        mat_size <- nr * nc
    } else {
        nr <- NaN
        nc <- NaN
        mat_size <- length(prox_mat)
    }
    data.table(
        nr = nr,
        nc = nc,
        mat_size = mat_size,
        min_prox = min(prox_mat, na.rm=T),
        mean_prox = mean(prox_mat, na.rm=T),
        med_prox = median(prox_mat, na.rm=T),
        max_prox = max(prox_mat, na.rm=T)
    )
}