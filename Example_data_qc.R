# This requires package caravagn/CNAqc
# devtools::install_github("caravagn/CNAqc")

library(CNAqc)
library(dplyr)

# We use some data released with the tool, this is WGS from a CRC biopsy of a primary tumour.
# We have generated for this sample somatic mutation, copy number segments and estimated the
# sample purity.

# -> somatic mutations
snvs = CNAqc::example_dataset_CNAqc$snvs

# -> copy number segment (CNA), precisely absolute calls for clonal segments
cna = CNAqc::example_dataset_CNAqc$cna

# -> a scalar reporting sample purity (0 <= x <= 1)
purity = CNAqc::example_dataset_CNAqc$purity

# -> reference genome used
reference = CNAqc::example_dataset_CNAqc$reference

# Check what we have in this data. To map mutations on top of CNA we need to have genome
# locations, and to run a deconvolution with MOBSTER we need the VAF of the somatic mutation.
# These calls are obtained by substracting the set of somatic mutations from 
snvs

# For instance we see:
# 
# chr - chromosome name (format chrXX)
# from/ to - nucleotide position (relative to the chromosome beginning)
# ref/ alt - reference and alternative allele
# 
# These mutations span all the genome and are just SNVs (it is good to check those)
snvs %>% group_by(chr) %>% summarise(n = n()) %>% arrange(desc(n)) # Mapping of the muts per chr
snvs %>% group_by(ref, alt) %>% summarise(n = n()) %>% arrange(desc(n)) # Type of mutations

# We also have read counts data obtained from the caller that was used to generate. You
# can verify that VAF = NV/DP. Eg, 6/60 = 0.1
snvs %>% select(DP, NV, VAF)

# CNAqc data object
x = CNAqc::init(snvs,
                cna,
                purity,
                reference)

# S3 print
print(x)

# Some basic plotting functions
plot_data_histogram(x)
plot_data_histogram(x, which = 'DP')

# You can play around with the CNAqc functions following examples available at https://caravagn.github.io/CNAqc/
# 
# Here we just check that the quality of the calls is GOOD to study this tumour. We use a method based on a peak
# detection algorithm that is derived from a combinatorial arguments that links tumour ploidy, mutation multiplicity,
# sample purity, and VAF.
x = analyze_peaks(x, matching_strategy = 'closest')

print(x)

# We can visualise the result of this analysis
plot_peaks_analysis(x)

# We can proceed further and compute CCF for this sample
x = compute_CCF(x)

plot_CCF(x)
plot_data_histogram(x, which = 'CCF')


