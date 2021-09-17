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

# run a few comparisons to test differences b/w distributions
compare_degree_distributions <- function(deg, gene_list) {
    idx <- which(names(deg) %in% gene_list)
    x <- deg[idx]
    y <- deg[-idx]
    ks <- ks.test(x, y)
    tt <- t.test(x, y)
    rs <- wilcox.test(x, y)
    data.table(
        ks_stat = ks$statistic,
        ks_pval = ks$p.value,
        t_stat = tt$statistic,
        t_pval = tt$p.value,
        rs_stat = rs$statistic,
        rs_pval = rs$p.value
    )
}

# sample from a background to reflect the degree distribution in a test gene list
sample_degree <- function(deg, deg_binned, gene_list) {
    n <- length(gene_list)
    idx <- which(names(deg) %in% gene_list)
    # how often does each bin occur in the gene list?
    bin_freq <- table(deg_binned[idx])
    bin_freq <- bin_freq[bin_freq > 0]
    # for each non-zero entry, sample that number of genes from
    # the corresponding bin
    Reduce(c, lapply(names(bin_freq), function(x) {
        freq <- bin_freq[[x]]
        idx <- which(deg_binned==x)
        sample(names(deg)[idx], freq, replace=F)
    }))
}

# for chunking a vector
chunk_vector <- function (x, chunk_size = 4) {
    stopifnot(is.vector(x))
    stopifnot(is.null(dim(x)))
    stopifnot(length(x) > 0)
    n <- length(x)
    l <- list()
    i <- 1
    while (i <= n) {
        l[[length(l) + 1]] <- x[i:min((i + chunk_size - 1), n)]
        i <- i + chunk_size
    }
    return(l)
}
