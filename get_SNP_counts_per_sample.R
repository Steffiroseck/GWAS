# Load necessary library
library(data.table)

# Set the path to your PED file
ped_file <- "C:/Users/afbi-roses/SHEEP_GENOTYPES/GWAS_meat_set1_84/Plink_files/Meatset_84.ped"

# Read the ped file; note that the PED file doesn't have a header.
# The first six columns are sample info; the rest are genotype calls (two columns per SNP).
ped_data <- fread(ped_file, header = FALSE, data.table = FALSE)

# Number of samples
num_samples <- nrow(ped_data)
print(paste0("The number of samples in the .ped file are: ", num_samples)) 

# Total number of genotype columns (excluding the first 6 sample info columns)
num_genotype_columns <- ncol(ped_data) - 6
print(paste0("The number of genotype columns in the .ped file are: ", num_genotype_columns))

# The number of SNPs (assuming 2 columns per SNP)
num_snps <- num_genotype_columns / 2
print(paste0("Total number of SNPs in the .ped file are: ", num_snps))

# Function to count non-missing genotype calls for one sample row
count_genotypes <- function(geno_row) {
  # Genotypes are stored in pairs; here we treat each SNP as one observation.
  # We'll consider a SNP genotype as missing if both alleles are "0".
  geno_matrix <- matrix(geno_row, ncol = 2, byrow = TRUE)
  # Count SNPs where at least one allele is not missing.
  non_missing <- sum(!(geno_matrix[, 1] == "0" & geno_matrix[, 2] == "0"))
  return(non_missing)
}

# Initialize a vector to store counts for each sample
snp_counts <- numeric(num_samples)

# Loop over each sample (row) and calculate the count
for (i in 1:num_samples) {
  # Extract genotype columns for the sample
  geno_row <- as.character(ped_data[i, -(1:6)])
  snp_counts[i] <- count_genotypes(geno_row)
}

# Combine sample IDs (from column 2 of ped file, for example) with SNP counts
results <- data.frame(SampleID = ped_data[, 2],
                      NonMissingSNPs = snp_counts,
                      TotalSNPs = num_snps,
                      stringsAsFactors = FALSE)

# Print the results
print(results)

