# Lab 1 : exploring the 1000 genomes (1KG) data

## Data access

- Download 1KG VCF files from `ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20130502/`

- Download `plink 1.9` from `https://www.cog-genomics.org/plink2`

- Convert VCF files to plink binary filesets (`.bed`, `.bim`, `.fam`). For instance:

```
./bin/plink --make-bed --hwe 1e-4 --maf 0.05 --snps-only --vcf vcf/chr21.vcf.gz --out plink/chr21
```

- You can get a rough idea of LD blocks using `https://bitbucket.org/nygcresearch/ldetect-data/get/ac125e47bf7f.zip`

---

## Exercise

1. Play with `plink` and `vcftools` (optional): e.g., reading and writing a subset of data

2. Load the genotype data in your favorite environment and play with it. I included a simple `R` script `demo1.R` as an example.  To run this script, you need to install several packages, including the one I made `util68881`.  You will find more instructions under `../util/`



