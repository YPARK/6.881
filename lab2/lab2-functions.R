library(tidyr)   # generally useful library
library(dplyr)   # generally useful library
library(ggplot2) # much prettier plottting
library(ggrepel) # avoid cluttered labels

## Remove NAs replacing them with zeros
rm.na.zero <- function(mat) replace(mat, is.na(mat), 0)

## Generate a matrix with random N(0,1)
rnorm.mat <- function(n, m) matrix(rnorm(n * m), nrow = n, ncol = m)

###########################
## simulation parameters ##
###########################

## n.causal <- 3  # between [1, n.snps]
## h2 <- 0.17     # between [0, 1] this is assumed heritability

simulate.phenotype <- function(X, n.causal, h2, n.pheno = 1) {

    Y <- matrix(NA, nrow = nrow(X), n.pheno)
    causal <- NULL

    for(j in 1:n.pheno) {
        causal.pos <- sample(n.snps, n.causal) # causal variant
        effect <- rnorm.mat(n.causal, 1) * sqrt(h2 / n.causal) # effect size
        y.hat <- X[, causal.pos, drop = FALSE] %*% effect
        y.err <- rnorm.mat(n.ind, 1) * sqrt(1 - h2)
        y <- y.hat + y.err
        
        Y[, j] <- y
        causal <- rbind(causal, data.frame(x.pos = causal.pos, y.col = j))
    }
    
    ret <- list(y = Y, causal = causal)
}

################################################################
## Calculation of SNP-by-SNP QTL statistics

## convert z-score to p-values (two-sided test)
zscore.pvalue <- function(z) {
    2 * pnorm(abs(z), lower.tail = FALSE)
}

## calculate univariate effect sizes and p-values
calc.qtl.stat <- function(xx, yy) {

    require(dplyr)
    require(tidyr)

    .xx <- scale(xx)
    .yy <- scale(yy)

    ## cross-product is much faster than covariance function
    n.obs <- crossprod(!is.na(.xx), !is.na(.yy))
    beta.mat <- crossprod(.xx %>% rm.na.zero(), .yy %>% rm.na.zero()) / n.obs

    ## residual standard deviation
    resid.se.mat <- matrix(NA, ncol(.xx), ncol(.yy))

    for(k in 1:ncol(.yy)) {

        beta.k <- beta.mat[, k]
        yy.k <- .yy[, k]
        err.k <- sweep(sweep(.xx, 2, beta.k, `*`), 1, .yy, `-`)
        se.k <- apply(err.k, 2, sd, na.rm = TRUE)

        resid.se.mat[, k] <- se.k
    }

    ## organize as consolidated table
    y.cols <- 1:ncol(yy)
    colnames(beta.mat) <- y.cols
    colnames(n.obs) <- y.cols
    colnames(resid.se.mat) <- y.cols

    beta.tab <- beta.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'beta', y.cols)
    
    resid.se.tab <- resid.se.mat %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'resid.se', y.cols)
    
    nobs.tab <- n.obs %>%
        as.data.frame() %>%
            dplyr::mutate(x.col = 1:n()) %>%
                tidyr::gather(key = 'y.col', value = 'n', y.cols)
    
    out.tab <- beta.tab %>%
        left_join(nobs.tab) %>%
            left_join(resid.se.tab) %>%
                dplyr::mutate(se = resid.se/sqrt(n)) %>%
                    dplyr::mutate(p.val = zscore.pvalue(beta/se))
    
    return(out.tab)
}

################################################################
## simply overlap with SNP information
annotate.snps <- function(qtl.stat.tab, plink.bim) {
    require(dplyr)

    ret <- plink.bim %>%
        dplyr::mutate(x.col = 1:n()) %>%
        right_join(qtl.stat.tab)
}

