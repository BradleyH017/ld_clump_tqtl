# Bradley March 2025
# Check LD
# Commands to compute (bash)
module load HGI/softpack/groups/macromapsqtl/macromapsqtlR4/10
gene=ENSG00000143801
gene_name=PSEN2
awk -v g=$gene '$1 == g' input/variants_per_gene_list.txt | awk '{print $2}' | uniq > temp/${gene}_variants.txt
num_chr=$(awk -F':' 'NR==1 {gsub("chr", "", $1); print $1}' temp/${gene}_variants.txt)
plink --bfile /lustre/scratch126/humgen/projects/sc-eqtl-ibd/data/genotypes/eQTL_genotypes_march_2024/plink_march_2025/plink_geno_chr${num_chr}_plink --extract temp/${gene}_variants.txt --make-bed --out temp/${gene}
awk -v g="${gene}" 'NR==1 || $1 == g' input/variants_per_gene_list.txt >> temp/${gene}_assoc.txt
plink --bfile temp/${gene} --clump temp/${gene}_assoc.txt --clump-p1 1 --clump-p2 1 --clump-r2 0.5 --clump-kb 1000 --out temp/${gene}
plink --bfile temp/${gene} --ld-snp-list temp/${gene}_variants.txt --ld-window-kb 50000000 --ld-window-r2 0 --r2 --out temp/clump_index_ld 

module load HGI/softpack/groups/macromapsqtl/macromapsqtlR4/10
Rscript scripts/plot_ld.r
mkdir -p temp/$gene_name
mv temp/*.* temp/$gene_name