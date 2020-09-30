# This requires package caravagn/mobster
# devtools::install_github("caravagn/mobster")

library(mobster)
library(dplyr)

# We use some data released with the tool, this is WGS from the paper
# 
# Nik-Zainal, Serena, et al. "The life history of 21 breast cancers." Cell 149.5 (2012): 994-1007.
# 
# https://www.cell.com/fulltext/S0092-8674(12)00527-2
# 
# Since the package wraps the result of the analysis with MOBSTER, we first extract the data
# and use it again to re-compute the analysis from scratch.

# -> somatic mutations data
input_calls = mobster::PD4120a_breast_sample$best$data %>% 
  select(-cluster) %>% 
  mutate(DP = round(NV/VAF)) %>% 
  select(chr, from, to, ref, alt, DP, NV, VAF, everything())

print(input_calls %>% as_tibble())

# Note what data we are using, SNVs that map on chromosome 3. We have selected this particular
# chromosome because it contains enough mutations, is largely diploid (so no need for CCF computation)
# and the calls seem very reliable (assessed elsewhere).
input_calls %>% group_by(chr) %>% summarise(n = n()) %>% arrange(desc(n)) # Mapping of the muts per chr
input_calls %>% group_by(ref, alt) %>% summarise(n = n()) %>% arrange(desc(n)) # Type of mutations

# We can visualise the data using ggplot
input_calls %>% 
  ggplot(aes(VAF)) +
  geom_histogram(binwidth = 0.01)

# We use one of the default parametrisations of the tool. This called FAST because it tests a few
# model parameters, and is usually fast to run. It is not necessarily the best analysis but in general
# is good enough to get a quick grasp of what is in the data.
mobster_analysis = mobster_fit(input_calls, seed = 1234, parallel = FALSE)

# Mobster returns a composite object (which we will eventually wrap up in an S3 object). So far you
# can see the fields returned by the tool in a named list.
names(mobster_analysis)

# A table reporting the score for each possible model returned by the tool
print(mobster_analysis$fits.table)

# The scoring used (defined in the paper): default reICL (~BIC + entropy of the latent variables)
print(mobster_analysis$model.selection)

# The best model fit -- what you will consider most of the times for downstream analysis
print(mobster_analysis$best)

# Plot the model fit (important! Check always your data fits!)
plot(mobster_analysis$best)

# Probability of assignment of a mutation to a cluster (check if some clusters are more difficult to define than others)
plot_latent_variables(mobster_analysis$best)

# How big is each cluster (normalised for total mutational burden)
plot_mixing_proportions(mobster_analysis$best)

# The likelihood function that we minimise
plot_NLL(mobster_analysis$best)

# Inspect alternative models - important if you do not like the model fit and you think that
# the model-selection procedure did not work well.
plot_model_selection(mobster_analysis)

# Evolutionary analyses
evolutionary_parameters(mobster_analysis)

# Clone trees -- will not work because we do not have drivers annotated
mobster::get_clone_trees(mobster_analysis$best)

# ... let's annotate some fake drivers ...
fake_annotations = mobster_analysis$best$data %>% 
  mutate(
    is_driver = FALSE,
    driver_label = paste(chr, from, to, ref, alt, sep = ':')
  ) %>% 
  as_tibble()

fake_annotations$is_driver[c(1, 2, 11)] = TRUE

# See where we annotated the drivers
print(fake_annotations %>%
        select(cluster, is_driver, driver_label) %>%
        filter(is_driver))

fake_fit = mobster_analysis
fake_fit$best$data = fake_annotations

# Compute the trees with the ctree package (called by mobster)
trees = mobster::get_clone_trees(fake_fit$best)

# Assemble all plots available for this sample
ggpubr::ggarrange(
  plot(trees[[1]]), 
  plot(trees[[2]]), 
  ncol = 2)

# Conclusions: with these data we cannot be 100% confident about the 
# lineage relation (linear vs branched) for C2/C3. We can dissect that
# by looking at other mutations (see the original Cell paper).
