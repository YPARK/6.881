################################################################
## This is how I prepared data:
## library(util6881)
## plink <- read.plink('../../1KG_phase3/eur-low/chr22')
## X <- plink$BED
## maf <- apply(X/2, 2, mean, na.rm = TRUE)   # average frequency
## maf[maf > .5] <- 1 - maf[maf > .5]         # Why do we flip?
## ## Remove MAF < 0.05%
## rm.snps <- which(maf < 5e-3 | is.na(maf))
## plink$BIM <- plink$BIM[-rm.snps, , drop = FALSE]
## plink$BED <- plink$BED[, -rm.snps, drop = FALSE]
## save(plink, file = 'chr22_lowfreq.Rd', compression_level = 9)
################################################################

load('chr22_lowfreq.Rd')

library(SKAT)    # install SKAT by CRAN
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

source('../lab2/lab2-functions.R') # use functions in the previous lab

log.msg <- function(...) {
    ss <- as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

################################################################
## a. Calculate MAF from 0/1/2 doage matrix
## b. visualize histogram
take.maf <- function(X) {
    ret <- apply(X/2, 2, mean, na.rm = TRUE)
    .flip <- which(ret > .5)
    ret[.flip] <- 1 - ret[.flip]
    return(ret)
}

MAF <- take.maf(plink$BED)

colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'a1', 'a2')

snp.tab <- plink$BIM %>% mutate(maf = MAF) %>%
    mutate(pos = 1:n())

pdf(file = 'figure_lab3_maf.pdf')
hist(MAF, 20)
dev.off()


################################################################
## Take subset of ExAC annotation
exac.types <- 'ccciiiiiiidddddddd'

exac.subset <- read_tsv('ExAC_data.txt', col_types = exac.types) %>%
    filter(chr == 22) %>%
    mutate(syn.weights = n_syn / exp_syn,
           mis.weights = n_mis / exp_mis,
           lof.weights = n_lof / exp_lof) %>%
    as.data.frame()

set.seed(1227)
n.causal.genes <- 2
n.causal.snps.per.gene <- 10
h2 <- 0.25

################################################################
## Select causal genes with higher burden of synonymous mutations
causal.gene.tab <- exac.subset %>%
    sample_n(size = n.causal.genes, weight = syn.weights)

## Select causal SNPs inversely proportional to MAF
sample.snps.g <- function(g, .gene.tab, .snp.tab, n.snps) {
    g.start <- .gene.tab[g, 'tx_start'] %>% as.integer()
    g.end <- .gene.tab[g, 'tx_end'] %>% as.integer()

    ret <- .snp.tab %>%
        filter(snp.loc >= g.start, snp.loc <= g.end)

    if(nrow(ret) < 1) return(NULL)

    .n.snps <- min(nrow(ret), n.snps)

    ret <- ret %>%
        sample_n(size = .n.snps, weight = 1/maf)
}

causal.snps.tab <- lapply(1:n.causal.genes, sample.snps.g,
                          .gene.tab = causal.gene.tab,
                          .snp.tab = snp.tab,
                          n.snps = n.causal.snps.per.gene) %>%
    bind_rows()

################################################################
## Simulate genetic associations on the low frequency variants
## Note: for low MAF, scaling is unstable, just do centering
X.causal <- plink$BED[ , causal.snps.tab$pos] %>%
    scale(center = TRUE, scale = FALSE)

## Note: a simple linear model
y.causal <- X.causal %*% rnorm.mat(nrow(causal.snps.tab), 1)

## Properly scale noise by heritability
var.genetic <- var(y.causal, na.rm = TRUE) %>% as.numeric()
var.noise <- var.genetic * (1/h2 - 1)
## should have h2 = var.genetic / (var.genetic + var.noise)

y.obs <- y.causal + rnorm.mat(nrow(y.causal), 1) * sqrt(var.noise)

################################################################
## Do GWAS SNP by SNP
gwas.stat <- calc.qtl.stat(plink$BED, y.obs) ## will take some time

gwas.stat <- gwas.stat %>% select(-y.col) %>% ## only one phenotype
    rename(pos = x.col) %>%                   ## match by position
    right_join(snp.tab)

## Show results
plt <-
    ggplot(gwas.stat, aes(x = snp.loc, y = -log10(p.val))) +
    theme_bw() +
    geom_point(alpha = 0.5) +
    ylab('-log10 p-value') + xlab('genomic location (kb)') +
    scale_x_continuous(labels = function(x) format(round(x/1e3), big.mark = ','))

causal.tab <- gwas.stat %>% filter(snp.loc %in% causal.snps.tab$snp.loc)

plt <- plt + geom_point(data = causal.tab, aes(x = snp.loc, y = -log10(p.val)),
                        color = 'red', pch = 22, size = 4)

ggsave(filename = 'figure_lab3_gwas.png', plot = plt, width = 8, height = 4)

################################################################
## Do gene by gene rare variant burden test
run.skat.g <- function(g, xx, .gene.tab, .snp.tab, .skat.null) {
    require(SKAT)

    g.start <- .gene.tab[g, 'tx_start'] %>% as.integer()
    g.end <- .gene.tab[g, 'tx_end'] %>% as.integer()
    x.pos <- .snp.tab %>% filter(snp.loc >= g.start, snp.loc <= g.end) %>%
        select(pos) %>% unlist(use.names = FALSE) %>% as.integer()

    if(length(x.pos) < 1) return(NULL)

    skat.out <- SKAT(xx[ , x.pos, drop = FALSE], obj = .skat.null,
                     kernel = 'linear.weighted', method = 'SKATO',
                     is_dosage = TRUE)

    ret <- .gene.tab[g, ] %>% select(transcript, gene) %>%
        mutate(n.snps = skat.out$param$n.marker,
               p.val = skat.out$p.value)

    log.msg('Finished: %d, p-value = %.2e', g, ret$p.val)
    return(ret)
}

genes <- 1:nrow(exac.subset)

skat.null <- SKAT_Null_Model(y.obs ~ 1, out_type = 'C')

skat.result <- lapply(genes, run.skat.g, xx = plink$BED,
                      .gene.tab = exac.subset, .snp.tab = snp.tab,
                      .skat.null = skat.null) %>% bind_rows()

skat.result <- skat.result %>% left_join(exac.subset)

skat.causal <- skat.result %>% filter(gene %in% causal.gene.tab$gene)

## Show gene-based results
plt <-
    ggplot(skat.result, aes(x = tx_start, y = -log10(p.val))) +
    theme_bw() +
    geom_point(alpha = 0.5) +
    geom_point(data = skat.causal, color = 'red', pch = 22, size = 4) +
    geom_label_repel(data = skat.causal, aes(label = gene)) +
    ylab('-log10 p-value') + xlab('genomic location (kb)') +
    scale_x_continuous(labels = function(x) format(round(x/1e3), big.mark = ','))

ggsave(filename = 'figure_lab3_skat.png', plot = plt, width = 8, height = 4)

