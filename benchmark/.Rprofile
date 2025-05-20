source("renv/activate.R")

# library(tools)
# source(file_path_as_absolute("renv/activate.R"))
# lib_path <- file_path_as_absolute("conda/lib")
Sys.setenv(LD_LIBRARY_PATH = paste("conda/lib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
library(reticulate)

use_python(paste0(getwd(),"/conda/bin/python"))
print(.libPaths())
print(py_config())
print(Sys.getenv("LD_LIBRARY_PATH"))

import("scanpy")
import("anndata")
# options(bitmapType = "cairo")


.convert <- function(obj) {
    builtins <- reticulate::import_builtins()
    Mapping <- reticulate::import("typing")$Mapping
    DataFrame <- reticulate::import("pandas")$DataFrame
    issparse <- reticulate::import("scipy.sparse")$issparse
    isinstance <- builtins$isinstance

    if (issparse(obj)) {
        obj <- obj$tocoo()
        return(Matrix::sparseMatrix(
            i = as.vector(reticulate::py_to_r(obj$row)) + 1,
            j = as.vector(reticulate::py_to_r(obj$col)) + 1,
            x = as.vector(reticulate::py_to_r(obj$data)),
            dims = unlist(as.vector(reticulate::py_to_r(obj$shape)))
        ))
    }
    if (isinstance(obj, DataFrame)) {
        obj <- reticulate::py_to_r(obj)
        attr(obj, "pandas.index") <- NULL
        return(obj)
    }
    if (!isinstance(obj, Mapping)) {
        return(reticulate::py_to_r(obj))
    }
    ret <- list()
    for (key in builtins$list(obj$keys())) {
        ret[[key]] <- .convert(obj[[key]])
    }
    return(ret)
}


convert.h5ad <- function(adata) {
    X <- .convert(adata$X)
    layers <- .convert(adata$layers)
    obs <- .convert(adata$obs)
    var <- .convert(adata$var)
    obsm <- .convert(adata$obsm)
    varm <- .convert(adata$varm)
    obsp <- .convert(adata$obsp)
    varp <- .convert(adata$varp)
    uns <- .convert(adata$uns)
    rownames(X) <- rownames(obs)
    colnames(X) <- rownames(var)
    return(list(
        X = X, layers = layers,
        obs = obs, var = var,
        obsm = obsm, varm = varm,
        obsp = obsp, varp = varp,
        uns = uns
    ))
}


read.h5ad <- function(filename) {
    ad <- reticulate::import("anndata", convert = FALSE)
    adata <- ad$read_h5ad(filename)
    return(convert.h5ad(adata))
}
