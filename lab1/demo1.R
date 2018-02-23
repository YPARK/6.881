
## Show LD blocks, comparing different populations
##
## LD blocks downloaded from Joe Pickrell's repository
## https://bitbucket.org/nygcresearch/ldetect-data/get/ac125e47bf7f.zip
##

library(dplyr)
library(util6881)
library(readr)
library(Matrix)

`%&&%` <- function(a, b) paste(a, b, sep = '')

take.ld.file <- function(pop, chr) {
    'ldblocks/' %&&% pop %&&% '/fourier_ls-chr' %&&% chr %&&% '.bed'
}

.chr <- 21
eur.ld.tab <- take.ld.file('EUR', chr = .chr) %>% read_tsv()

eas.plink <- '1KG_phase3/eas/chr' %&&% .chr %>% read.plink()
eur.plink <- '1KG_phase3/eur/chr' %&&% .chr %>% read.plink()

################################################################
## Take a subset of SNPs
subset.plink <- function(.plink, ld.start, ld.stop) {
    .idx <- which(.plink$BIM[, 4] >= ld.start & .plink$BIM[, 4] <= ld.stop)
    ret <- .plink
    ret$BIM <- ret$BIM[.idx, , drop = FALSE]
    ret$BED <- ret$BED[, .idx, drop = FALSE]
    colnames(ret$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')
    return(ret)
}

################################################################
## match plink.other to plink.ref
match.plink <- function(plink.ref, plink.other) {

    if(is.null(plink.ref)) return(NULL)
    if(is.null(plink.other)) return(NULL)

    ret.ref <- plink.ref
    ret.other <- plink.other

    gwas.bim <- plink.ref$BIM %>%
        mutate(gwas.x.pos = 1:n()) %>%
            rename(gwas.plink.a1 = plink.a1,
                   gwas.plink.a2 = plink.a2) %>%
                       select(-missing)

    other.bim <- plink.other$BIM %>%
        mutate(other.x.pos = 1:n()) %>%
            rename(other.plink.a1 = plink.a1,
                   other.plink.a2 = plink.a2,
                   other.rs = rs) %>%
                       select(-missing)

    bim.matched <- gwas.bim %>%
        left_join(other.bim) %>%
            na.omit()

    if(nrow(bim.matched) < 1) return(NULL)

    bim.matched <- bim.matched %>%
        dplyr::filter(((gwas.plink.a1 == other.plink.a1) & (gwas.plink.a2 == other.plink.a2)) |
                          ((gwas.plink.a2 == other.plink.a1) & (gwas.plink.a1 == other.plink.a2))) %>%
                              arrange(chr, snp.loc)

    if(nrow(bim.matched) < 1) return(NULL)

    ret.ref$BIM <- ret.ref$BIM[bim.matched$gwas.x.pos, , drop = FALSE]
    ret.ref$BED <- ret.ref$BED[ , bim.matched$gwas.x.pos, drop = FALSE]

    ret.other$BIM <- ret.other$BIM[bim.matched$other.x.pos, , drop = FALSE]
    ret.other$BED <- ret.other$BED[ , bim.matched$other.x.pos, drop = FALSE]

    flip.tab <- ret.ref$BIM %>% mutate(gwas.x.pos = 1:n()) %>%
        left_join(ret.other$BIM %>% mutate(other.x.pos = 1:n()),
                  by = c('chr', 'snp.loc'),
                  suffix = c('.ref', '.other')) %>%                      
                      filter(plink.a1.ref != plink.a1.other)

    ret.other$BIM[flip.tab$other.x.pos, ] <- ret.ref$BIM[flip.tab$gwas.x.pos, ]

    flip.bed <- ret.other$BED[, flip.tab$other.x.pos]
    zero.idx <- flip.bed <= 0.5
    two.idx <- flip.bed >= 1.5
    flip.bed[two.idx] <- 0
    flip.bed[zero.idx] <- 2
    ret.other$BED[, flip.tab$other.x.pos] <- flip.bed

    return(list(ref = ret.ref, other = ret.other))
}

rm.na.zero <- function(mat) replace(mat, is.na(mat), 0)

fast.cov <- function(x, y) {
    n.obs <- crossprod(!is.na(x), !is.na(y))
    ret <- crossprod(replace(x, is.na(x), 0),
                     replace(y, is.na(y), 0)) / n.obs
    return(ret)
}

nrow(eur.ld.tab)

ld.idx <- 13

ld.start <- eur.ld.tab[ld.idx, 'start'] %>% as.integer()
ld.stop <- eur.ld.tab[ld.idx, 'stop'] %>% as.integer()

eas.sub.plink <- eas.plink %>% subset.plink(ld.start, ld.stop)
eur.sub.plink <- eur.plink %>% subset.plink(ld.start, ld.stop)

## match two plinks (making EUR as reference A1/A2)
matched <- match.plink(eur.sub.plink, eas.sub.plink)
eur.sub.plink <- matched$ref
eas.sub.plink <- matched$other

## 1. compare allele frequency

eur.mean <- apply(eur.sub.plink$BED, 2, mean, na.rm = TRUE)
eas.mean <- apply(eas.sub.plink$BED, 2, mean, na.rm = TRUE)

pdf(file = 'Fig_MAF_EUR_EAS_1.pdf')
par(mfrow = c(2, 1))
.plot <- function(...) plot(..., pch = 19, col = 'gray50', cex = .5)
.title <- .chr %&&% ':' %&&% ld.start %&&% '-' %&&% ld.stop
.plot(eur.mean, main = 'European frequency : '%&&% .title)
.plot(eas.mean, main = 'East Asian frequency : '%&&% .title)
dev.off()

pdf(file = 'Fig_MAF_EUR_EAS_2.pdf')
par(mfrow = c(1, 1))
.plot(eas.mean, eur.mean, xlab = 'East Asian frequency', ylab = 'Eurpean frequency')
abline(a=0, b=1, col = 2)
dev.off()

## 2. show correlation structure

eas.mat <- eas.sub.plink$BED[, 1:1000, drop = FALSE] %>% scale()
eur.mat <- eur.sub.plink$BED[, 1:1000, drop = FALSE] %>% scale()

eas.cov <- fast.cov(eas.mat, eas.mat)
eur.cov <- fast.cov(eur.mat, eur.mat)

png(filename = 'Fig_LD_EAS.png')
image(eas.cov)
dev.off()

png(filename = 'Fig_LD_EUR.png')
image(eur.cov)
dev.off()

