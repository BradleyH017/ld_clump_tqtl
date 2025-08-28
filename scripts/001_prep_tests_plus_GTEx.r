# Bradley March 2025
# module load HGI/softpack/users/wl2/fine_mapping/1

library(tidyverse)
library(qvalue)

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

basedir = sumstats.all.basedir
base_dir = basedir
sumstat.files <- list.files(sumstats.all.basedir, pattern="*Cis_eqtls_qval.tsv",
                              full.names=FALSE, recursive=TRUE)
  
list.of.sumstat.dfs <- lapply(sumstat.files, function(file.name) {
  # Get info from filename
  split.file.name <- str_split(file.name,'/')[[1]]
  name <- split.file.name[1]

  # Get info from name
  split.name <- str_split(name,'__')[[1]]
  annotation <- str_remove(split.name[2], '_all')
  method <- split.name[1]

  # Read file
  df <- read_tsv(paste0(base_dir, file.name))

  # Add new columns
  mutate(df,
          name = name,
          annotation = annotation,
          method = method,
          bonf_pvalue = pval_beta * n())
})
  
sumstat.df <- bind_rows(list.of.sumstat.dfs) %>% 
left_join(dplyr::select(grch38, ensgene, symbol, biotype),
            by = c("phenotype_id" = "ensgene"),multiple='first') %>% 
mutate(symbol = case_when(symbol == '' ~ phenotype_id,
                            TRUE ~ symbol)) %>% 
left_join(annot.mapping, by=c('annotation' = 'label_machine'))

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
        pheno_annotation %in% sig.gene.condition # No need to apply additional filter to pval_perm or pval_beta. 
    ) %>% 
    group_by(variant_id, phenotype_id) %>% 
    summarise(
        P=min(pval_nominal) # Not interested in whether the same gene-snp pair is more/less significant in different conditions
    ) %>% 
    mutate(
        type="eQTL"
    ) # 194233

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
    distinct() # 195,296


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
#             763            193495                6              1057


##################
# Specific to this script - group with GTEx too
##################
eqtlcat_meta = "/lustre/scratch125/humgen/resources_v2/eQTL_catalogue/dataset_metadata.tsv"
eqtldir = "/lustre/scratch125/humgen/resources_v2/eQTL_catalogue/sumstats"
eqtl_cat = read.delim(eqtlcat_meta)
gtex_ge_datasets = eqtl_cat %>%
    filter(
        study_id == "QTS000015",
        quant_method == "ge"
    ) %>% 
    pull(dataset_id)

# load permuted results per dataset, save an annotation of these
gtex = do.call(rbind, lapply(gtex_ge_datasets, function(x){
    # Get tissue
    tissue = eqtl_cat %>%
        filter(dataset_id == !!x) %>% 
        pull(sample_group)

    print(paste0("..Loading GTEx: ", tissue))

    # Load
    path = paste0(eqtldir, "/QTS000015/", x, "/", x, ".permuted.tsv.gz")
    temp = read.delim(path) %>%
        mutate(
            qvalue = qvalue(p_beta)$qvalues, # Apply qvalue correction to beta-regressed pvalue
            type = paste0("GTEx-", tissue)
        ) %>% 
        filter(
            qvalue < 0.05, # Filter for significant
            chromosome != "X", # Remove sex genes
            chromosome != "Y"
        ) %>% 
        select(
            molecular_trait_id,
            variant, 
            pvalue,
            type
        ) %>% 
        rowwise() %>% 
        mutate(
            variant = gsub("\\_", "\\:", variant)
        ) %>% 
        ungroup() %>% 
        rename(
            phenotype_id = molecular_trait_id,
            SNP = variant,
            P = pvalue
        )

    return(temp)
}))

# Merge with rest
ours_gtex = coqtl %>% 
    bind_rows(gtex) %>% 
    distinct()

# summary_genes
write.table(ours_gtex %>% pull(phenotype_id) %>% unique(), paste0(out.dir, "/gene_list.txt"), col.names=F, row.names=F, quote=F)
# variants per gene
write.table(ours_gtex, paste0(out.dir, "/variants_per_gene_list.txt"), col.names=T, row.names=F, quote=F, sep = "\t")