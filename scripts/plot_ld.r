# Read in results
library(dplyr)
library(ComplexHeatmap)
library(reshape2)
library(circlize)
ld = read.delim("temp_plot/clump_index_ld.ld", sep = "") %>% 
    filter(SNP_A != SNP_B)

# Prep for heatmap
snp_levels <- unique(c(ld$SNP_A, ld$SNP_B))  # Get unique SNPs
ld_matrix <- matrix(NA, nrow = length(snp_levels), ncol = length(snp_levels), dimnames = list(snp_levels, snp_levels))
for (i in 1:nrow(ld)) {
  ld_matrix[ld$SNP_A[i], ld$SNP_B[i]] <- ld$R2[i]
  ld_matrix[ld$SNP_B[i], ld$SNP_A[i]] <- ld$R2[i]  # Ensure symmetry
}
diag(ld_matrix) <- 1  # Set diagonal to 1 (self-LD)
ld_matrix <- as.matrix(ld_matrix)
ld_matrix[is.na(ld_matrix)] <- 0

col_fun <- colorRamp2(c(0, 1), c("white", "red"))

# Create heatmap
ppi=300
png("temp_plot/heatmap_ld.png", res=ppi, width=8*ppi, height=8*ppi)
Heatmap(ld_matrix,
        name = "R²",
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", ld_matrix[i, j]), x, y, gp = gpar(fontsize = 10))
        },
        cluster_rows = T,
        cluster_columns = T,
        clustering_distance_rows = "euclidean",  # Can be "correlation", "manhattan", etc.
        clustering_method_rows = "complete",    # Can be "average", "ward.D2", etc.
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "complete",
        show_row_names = TRUE,
        show_column_names = TRUE)
dev.off()