# Bradley March 2025
# Check LD
# Commands to compute (bash)
module load HGI/softpack/groups/macromapsqtl/macromapsqtlR4/10
gene=ENSG00000162613
gene_name=FUBP1
ld_thresh=0.5

mkdir -p temp_plot
awk -v g=$gene '$1 == g' input/variants_per_gene_list.txt | awk '{print $2}' | uniq > temp_plot/${gene}_variants.txt
num_chr=$(awk -F':' 'NR==1 {gsub("chr", "", $1); print $1}' temp_plot/${gene}_variants.txt)
plink --bfile /lustre/scratch127/humgen/projects_v2/sc-eqtl-ibd/core_analysis_output/IBDverse_multi-tissue_eQTL_project/IBDverse_genotypes/2024_07_11-genotype_plate12345/plink/imputed_chr${num_chr} --extract temp_plot/${gene}_variants.txt --make-bed --out temp_plot/${gene}
awk -v g="${gene}" 'NR==1 || $1 == g' input/variants_per_gene_list.txt >> temp_plot/${gene}_assoc.txt
plink --bfile temp_plot/${gene} --clump temp_plot/${gene}_assoc.txt --clump-p1 1 --clump-p2 1 --clump-r2 $ld_thresh --clump-kb 1000 --out temp_plot/${gene}
plink --bfile temp_plot/${gene} --ld-snp-list temp_plot/${gene}_variants.txt --ld-window-kb 50000000 --ld-window-r2 0 --r2 --out temp_plot/clump_index_ld 

module load HGI/softpack/groups/macromapsqtl/macromapsqtlR4/10
Rscript scripts/plot_ld.r
mkdir -p temp_plot/$gene_name
mv temp_plot/*.* temp_plot/$gene_name