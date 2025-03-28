# Bradley March 2025
library(tidyverse)

# Get paths
repo.dir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/code/IBDVerse-sc-eQTL-code/'
sumstats.all.basedir <- '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_03_20-multi_tissue_base_results/TensorQTL_eQTLS/'
out.dir <- 'input'
data.dir <- paste0(repo.dir,'/data/')
coloc.dir <- "/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/bradley_analysis/scripts/scRNAseq/snakemake_colocalisation/results/2025_03_IBDverse_coloc_all_gwas/collapsed/"
source(paste0(repo.dir,'qtl_plot/helper_functions.R'))

# Make outdirs
if(!file.exists(out.dir)){
    dir.create(out.dir, recursive=T)
}

##################
# Read in eQTLs and prep
##################
sumstat.df <- read_eqtls(sumstats.all.basedir)
sig.sumstat.df <- sumstat.df %>% 
  filter(qval < 0.05) 

# Get sig gene x condition pairs
sig.gene.condition = sig.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation)) %>% 
    pull(pheno_annotation) %>% 
    unique()

length(sig.gene.condition) # 330504

##################
# Read in conditional eQTLs and prep
##################
cond.sumstat.df <- read_conditional_eqtls(sumstats.all.basedir)
cond.sumstat.df <- cond.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation))

# get a unique list of phenotype IDs and variants
qtl_gene_leads = cond.sumstat.df %>% 
    filter(pheno_annotation %in% sig.gene.condition) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(pval_nominal) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="eQTL"
    )

##################
# Read in colocs and prep
##################

colocs_all = get_colocs(coloc.dir, known_ibd_only = F) %>% 
    filter(
        !(gwas_trait %in% c("ibd", "cd", "uc")),
        PP.H4.abf > 0.75
    )

# get a unique list of phenotype IDs and variants
coloc_gene_leads = colocs_all %>% 
    mutate(variant_id = gsub("\\_", "\\:", snp_id)) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(qtl_pval) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="coloc"
    )


##################
# Combine these together, saveas well as a summary of the genes
##################
coqtl = coloc_gene_leads %>% 
    bind_rows(qtl_gene_leads) %>% 
    arrange(phenotype_id, variant_id, P) %>% 
    rename(SNP = variant_id) %>% 
    select(phenotype_id, SNP, P, type)

nrow(coqtl) # pre add coloc=206708

# summary_genes
write.table(coqtl %>% pull(phenotype_id) %>% unique(), paste0(out.dir, "/gene_list.txt"), col.names=F, row.names=F, quote=F)
# variants per gene
write.table(coqtl, paste0(out.dir, "/variants_per_gene_list.txt"), col.names=T, row.names=F, quote=F, sep = "\t")


