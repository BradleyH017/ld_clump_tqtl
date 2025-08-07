# Bradley March 2025
library(tidyverse)

# Get paths
repo.dir <- '../IBDVerse-sc-eQTL-code/'
sumstats.all.basedir <- "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_11-multi_tissue_base_results/TensorQTL_eQTLS/" # Specify target path for the pipeline
sumstats.interaction.basedir <- '/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/repos/nf-hgi_qtlight/2025_06_12-multi_tissue_interaction_results/TensorQTL_eQTLS/'
out.dir <- 'input'
data.dir <- paste0(repo.dir,'/data/')
coloc.dir <- "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/bradley_analysis/IBDverse/snakemake_coloc/results/2025_06_11_IBDverse_coloc_all_gwas/collapsed/"
int.coloc.dir <- "/lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/analysis/tobi_qtl_analysis/results/coloc/IBDverse-multi_tissue_interaction_2025/collapsed"
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

length(sig.gene.condition) # 329996

##################
# Read in conditional eQTLs and prep
##################
cond.sumstat.df <- read_conditional_eqtls(sumstats.all.basedir)
cond.sumstat.df <- cond.sumstat.df %>% 
    mutate(pheno_annotation = paste0(phenotype_id, "-", annotation))

# Column names described in https://github.com/broadinstitute/tensorqtl/blob/0c4db65a0cdc47f3b824ae530b89d270ef5e0096/docs/outputs.md?plain=1#L51
# get a unique list of phenotype IDs and variants
qtl_gene_leads = cond.sumstat.df %>% 
    filter(
        pheno_annotation %in% sig.gene.condition,
        pval_beta < 0.05 
    ) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(pval_nominal) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="eQTL"
    ) # 188357

# Load in interaction eQTLs
interactions = read_ieqtls(sumstats.interaction.basedir) %>% 
    filter(
        pval_adj_bh < 0.05
    ) %>% 
    mutate(
        variant_id = factor(variant_id, levels=unique(variant_id)),
        interaction_new = interaction.mapping[interaction]
    )

int_gene_leads = interactions %>%
    group_by(variant_id, phenotype_id) %>%
    filter(pval_adj_bh < 0.05) %>%  
    summarise(
        P=min(pval_adj_bh)
    ) %>% 
    mutate(
        type="interaction_eQTL"
    )

# Add this
qtl_gene_leads = qtl_gene_leads %>% bind_rows(int_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() # 189,422


##################
# Read in colocs and prep
##################

colocs_all = get_colocs(coloc.dir, known_ibd_only = T)  %>% # All are IBD/CD/UC
    filter(PP.H4.abf > 0.75) # 6051

# get a unique list of phenotype IDs and variants
coloc_gene_leads = colocs_all %>% 
    mutate(variant_id = gsub("\\_", "\\:", snp_id)) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(qtl_pval) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="coloc"
    ) # 1,714

# Interactions: 
int_colocs_all = get_colocs(int.coloc.dir, known_ibd_only = T) %>% 
    filter(
        PP.H4.abf > 0.75
    )

int_coloc_gene_leads = int_colocs_all %>% 
    mutate(variant_id = gsub("\\_", "\\:", snp_id)) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(qtl_pval) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="interaction_coloc"
    ) # 6

# Combine
coloc_gene_leads = coloc_gene_leads %>% bind_rows(int_coloc_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() # 1720

##################
# Combine these together, saveas well as a summary of the genes
##################
coqtl = coloc_gene_leads %>% 
    bind_rows(qtl_gene_leads) %>% 
    group_by(variant_id, phenotype_id) %>% 
    slice_min(P) %>% 
    distinct() %>% 
    arrange(phenotype_id, variant_id, P) %>% 
    rename(SNP = variant_id) %>% 
    select(phenotype_id, SNP, P, type)

nrow(coqtl) # 189946
table(coqtl$type)
#            coloc              eQTL interaction_coloc  interaction_eQTL 
#             763            187621                6              1059

# summary_genes
write.table(coqtl %>% pull(phenotype_id) %>% unique(), paste0(out.dir, "/gene_list.txt"), col.names=F, row.names=F, quote=F)
# variants per gene
write.table(coqtl, paste0(out.dir, "/variants_per_gene_list.txt"), col.names=T, row.names=F, quote=F, sep = "\t")