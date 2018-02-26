
options(stringsAsFactors = FALSE) # Don't like R treat string as a factor
source('lab2-functions.R')  # See the separate file

## NOTE: We use "%>%" pipeline operator a lot!  It helps write a code
## with less parenthesis (provided by dplyr)

################################################################
## 0. For convenience, we saved chr22 data in an R-friendly way
##
## library(util6881)
## plink <- read.plink('../lab1/1KG_phase3/eur/chr22')
## save(plink, file = 'chr22.Rd')
##

load('chr22.Rd') # simply load the R data file (plink)

colnames(plink$BIM) <- c('chr', 'rs', 'missing', 'snp.loc', 'plink.a1', 'plink.a2')

################################################################
## 1. Simulate quantitative traits
## Standardize genotype matrix and replace NAs with zeros
## check dim(X): individual x SNP matrix
X <- plink$BED %>% scale() %>% rm.na.zero()
n.snps <- ncol(X)
n.ind <- nrow(X)

## Also have genotype matrix wihtout LD
X.noLD <- apply(X, 2, function(x) x[sample(n.ind)])

set.seed(1331) # some random seed to make it reproducible (try many different values)
y.data <- simulate.phenotype(X, n.causal = 2, h2 = 0.17, n.pheno = 1)

y <- y.data$y

## Calculate observed QTL association statistics
qtl.stat <- calc.qtl.stat(X, y) %>%
    annotate.snps(plink$BIM)

## Just perform sample permutation
y.null <- y[sample(n.ind), , drop = FALSE]

## Calculate observed QTL null statistics
qtl.stat.null <- calc.qtl.stat(X, y.null) %>%
    annotate.snps(plink$BIM)

################################################################
## 2. show QQ-plot in log10 scale

plt.df <- qtl.stat %>%
    filter(y.col == 1) %>% # just do it on the first phenotype
    arrange(p.val) %>%
    rename(observed = p.val)

n.tot <- nrow(plt.df)

plt.df <- plt.df %>%
    mutate(expected = (1:n()) / (n.tot + 1)) # Why is this expected under the null?

## highlight the causal ones
causal.pos <- y.data$causal %>%
    filter(y.col == 1) %>%
    select(x.pos) %>%
    unlist(use.names = FALSE)

plt.df.causal <- plt.df %>% filter(x.col %in% causal.pos)

plt1 <-
    ggplot(plt.df, aes(x = -log10(expected), y = -log10(observed))) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_point() +
    geom_point(data = plt.df.causal, pch = 21, color = 'red', size = 5) +
    geom_text_repel(data = plt.df.causal, aes(label = rs), nudge_x = .5)

ggsave(filename = 'figure_lab2_qqplot.png', plot = plt1)

################################################################
## 3. show Manhattan plot

plt2 <-
    ggplot(plt.df, aes(x = snp.loc, y = -log10(observed))) +
    theme_bw() +
    geom_point(alpha = 0.5) +
    geom_point(data = plt.df.causal, pch = 21, color = 'red', size = 5) +
    geom_label_repel(data = plt.df.causal, aes(label = rs),
                     segment.color = 'red', nudge_x = .5, nudge_y = .5) +
    scale_x_continuous(labels = function(x) format(round(x/1e3), big.mark = ',')) +
    xlab('genomic location (kb)') +
    ylab('-log10 p-value')

ggsave(filename = 'figure_lab2_manhattan.png', plot = plt2, width = 8, height = 4)
