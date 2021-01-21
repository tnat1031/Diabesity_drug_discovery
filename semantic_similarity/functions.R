# Supporting functions for semantic similarity analyses

# compute the jaccard index between vectors x and y
jaccard <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
}

# define a function to compute similarity between a disease gene list
# and a collection of drug gene lists
# this will simply wrap the clusterSim function from GOSemSim
# also enable support for simple jaccard index, which we can use to compare
# with semantic similarity metrics
# rewrite this using future package, which has better support for parallel
# access to SQLite objects such as hsGO (AnnotationDbi object)
disease_drug_sim <- function(disease_gene_list, drug_gene_lists, measure,
                             mc.cores=4, ...) {
    # use future to run computations in parallel
    plan(multicore, workers = mc.cores)
    res <- future.apply::future_lapply(drug_gene_lists, function(x) {
        tryCatch({
            if (measure=="jaccard") {
                # just compute jaccard index
                jaccard(disease_gene_list, x)
            } else {
                clusterSim(disease_gene_list, x, measure=measure, ...)
                # catch any errors if they occur
            }}, error = function(e) {
                message(e, " ", x, "\n")
                return(NA)
            })
    })
    names(res) <- names(drug_gene_lists)
    # unlist so res becomes a vector
    return(unlist(res))
}

# function to generate null (permuted) gene lists
make_null_sets <- function(sets, universe=NULL, nperm=1e2) {
    if (is.null(universe)) {
        universe <- Reduce(union, sets)
    }
    lapply(1:nperm, function(i) {
        sapply(sets, function(x) {
            sample(universe, length(x))
        }, USE.NAMES = T, simplify = F)
    })
}