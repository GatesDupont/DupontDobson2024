# Load necessary libraries
library(ape)
library(auk)
library(parallel)
library(MuMIn)
library(tidyverse)
library(MRFtools)
library(readxl)
library(mgcv)

df_bcr <- list()
for(i in 1:6){
  
  # BCR names
  bird_conservation_region <- c("Atlantic Northern Forest","Northern Pacific Rainforest",
                                "Appalachian Mountains", "Eastern Tallgrass Prairie",
                                "Northern Rockies", "Shortgrass Prairie")[i]
  
  # ---- Load in model results ----
  
  # Read the combined list as an RDS file in your directory
  BCR_file_suffix <- gsub(" ", "_", toupper(bird_conservation_region))
  BCR_output_file <- paste("data/BCR-results/results_bcr_", BCR_file_suffix, ".rds", sep = "")
  combined_results <- readRDS(file = BCR_output_file)
  
  # ---- Broad Groups categories -----
  file_name <- "data/BCR-species/complete/species_bcr_APPALACHIAN_MOUNTAINS_categories.xlsx"
  categories_dataset <- paste("data/BCR-species/complete/species_bcr_", BCR_file_suffix, "_categories.xlsx", sep = "")
  broad_groups <- read_excel(categories_dataset)
  
  if("Sagebrush_obligates" %in% colnames(broad_groups)){
    broad_groups <- broad_groups %>% select(-Sagebrush_obligates)
  }
  
  
  # ---- Assign habitats to species ----
  
  # Pull out important data
  df_bcr[[i]] <- combined_results$SpeciesEstimates %>%
    # Append variance of estimates and make weights
    left_join(combined_results$SpeciesVariances, by = "species") %>%
    rename(Variance = var_growth_rate) %>%
    mutate(weight = 1 / Variance) %>% # Calculate weights
    mutate(weight = round(weight/min(weight), 0)) %>%
    # Keep only required data
    select(species, growth_rate = r_mean, Initial_abundance = N0_mean, weight) %>%
    # Append my species groups
    mutate(species_clean = sub(" \\(all forms\\)", "", species)) %>%
    left_join(broad_groups, by = "species") %>%
    mutate(Migratory_species = if_else(Residents == 1, 0, 1)) %>%
    na.omit()
  
}

# Combine the datasets                                                                                                                 
df <- do.call(rbind, df_bcr) %>%
  mutate(bcr = as.factor(bcr),
         species = as.factor(species))

# This is from the 03_drivers_combinedBCRs script
my_species <- df %>%
  select(species) %>%
  distinct() %>%
  mutate(species_clean = gsub(" \\(all forms\\)", "", species)) %>%
  mutate(sci_name = ebird_species(.$species_clean)) %>%
  mutate(sci_name_dash = gsub(" ", "_", sci_name))

# Read the phylogenetic tree - I think this is Clements 4 taxonomy
phylo_tree <- read.tree("/Users/gatesdupont/Desktop/Stage2_Hackett_MCC_no_neg.tre")

# Data frame of tree tip names (seems to be from Clements 4)
phylo_tree_spp <- data.frame(
  sci_name_dash = phylo_tree$tip.label,
  exists = 1)

# Converting 2023 eBird taxonomy to Clements 4 taxonomy, manually
species_names_hash <- my_species %>%
  mutate(sci_name_dash = gsub("Setophaga_americana", "Parula_americana", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Setophaga_citrina", "Wilsonia_citrina", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Setophaga", "Dendroica", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Dendroica_ruticilla", "Setophaga_ruticilla", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Poecile", "Parus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Haemorhous", "Carpodacus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Antrostomus", "Caprimulgus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Cardellina", "Wilsonia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Spinus", "Carduelis", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Dryobates", "Picoides", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Leiothlypis", "Vermivora", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Parkesia", "Seiurus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Melozone_fusca", "Pipilo_fuscus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Melozone", "Pipilo", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Rhynchophanes", "Calcarius", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Corthylio", "Regulus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Centronyx", "Ammodramus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_philadelphia", "Oporornis_philadelphia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Vermivora_cyanoptera", "Vermivora_pinus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Troglodytes_pacificus", "Troglodytes_troglodytes", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Troglodytes_hiemalis", "Troglodytes_troglodytes", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Spatula_cyanoptera", "Anas_cyanoptera", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Circus_hudsonius", "Circus_cyaneus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Artemisiospiza_nevadensis", "Amphispiza_belli", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Peucaea_cassinii", "Aimophila_cassinii", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Ardea_alba", "Casmerodius_albus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Cistothorus_stellaris", "Cistothorus_platensis", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_formosa", "Oporornis_formosus", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Geothlypis_tolmiei", "Oporornis_tolmiei", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Ixoreus_naevius", "Zoothera_naevia", sci_name_dash)) %>%
  mutate(sci_name_dash = gsub("Selasphorus_calliope", "Stellula_calliope", sci_name_dash)) %>%
  left_join(phylo_tree_spp) %>%
  select(species, tree_sci_name = sci_name_dash)

# Final data
df2 <- df %>%
  left_join(species_names_hash) %>%
  mutate(tree_sci_name = as.factor(tree_sci_name)) %>%
  mutate(Initial_abundance = scale(Initial_abundance)[,1])


# Pruned tree
pruned <- drop.tip(phylo_tree, setdiff(phylo_tree$tip.label, df2$tree_sci_name))

# Penalty matrix from phylogeny
pmat <- mrf_penalty(pruned, type = "individual")

# Make sure the scientific species names are in the same order as the tip labels
spp_names_ordered_ref <- df2 %>%
  mutate(tree_sci_name = factor(tree_sci_name, levels = pruned$tip.label)) %>%
  select(tree_sci_name, species) %>%
  distinct() %>%
  arrange(tree_sci_name) %>%
  pull(species)

# Make sure the scientific species names are in the same order as the tip labels
df2 <- df2 %>%
  mutate(tree_sci_name = factor(tree_sci_name, levels = pruned$tip.label)) %>%
  arrange(tree_sci_name) %>%
  mutate(species = factor(species, levels = spp_names_ordered_ref))


# ---- Global model ----

# Global model
m.global <- bam(growth_rate ~ Introduced_species + Aerial_insectivores + Initial_abundance +
                  Secondary_growth_specialists + Urban_adapted_species +
                  Forest_obligates + Grassland_specialists + Wetland_obligates + 
                  Obligate_migrants + Facultative_migrants + Residents + Migratory_species +
                  Tree_cavity_nesters + 
                  s(bcr, bs = "re") +
                  s(species, bs = "re") +
                  s(tree_sci_name, bs = "mrf", xt = list(penalty = pmat)),
                discrete = T, select = F, gamma = 1.4,
                weights = df2$weight, data = df2, na.action = "na.fail")

summary(m.global)



# ----- Use global model for model selection ----

# Set up the cluster and fit models
clusterType <- if(length(find.package("snow", quiet = TRUE))) "SOCK" else "PSOCK"
clust <- try(makeCluster(getOption("cl.cores", 2), type = clusterType))
clusterEvalQ(clust, library(mgcv))
options(na.action = "na.fail")
clusterExport(clust, c("pmat", "df2"))
model_set <- dredge(m.global, cluster = clust, trace = 2, 
                    fixed = ~ s(bcr, bs = "re") + s(species, bs = "re") + s(sci_name, bs = "mrf", xt = list(penalty = pmat)),
                    subset = ~( !("Migratory_species" & "Residents") &
                                  !("Migratory_species" & "Facultative_migrants") &
                                  !("Migratory_species" & "Obligate_migrants") ))
stopCluster(clust)

# ---- Work with the top model ----

m <- get.models(model_set, 1)[[1]]
summary(m)
m_coefs <- summary(m, re.test = F)$p.coeff

# Calculate confidence intervals for plotting
results <- data.frame(
  term = names(m_coefs),
  estimate = as.numeric(m_coefs),
  stderr = summary(m, re.test = F)$se[1:length(m_coefs)]) %>%
  mutate(lower = estimate - 1.96 * stderr,
         upper = estimate + 1.96 * stderr) %>%
  `rownames<-`(NULL) %>%
  mutate(term = if_else(term == "(Intercept)", "Baseline", term)) %>%
  mutate(term = gsub("_", " ", term)) %>%
  mutate(significant = if_else(lower <= 0 & upper >= 0, "not significant", "significant")) %>%
  mutate(direction = if_else(estimate > 0, "positive", "negative")) %>%
  mutate(direction = case_when(
    direction == "positive" & significant == "not significant" ~ "Low certainty",
    direction == "negative" & significant == "not significant" ~ "Low certainty",
    direction == "positive" & significant == "significant" ~ "Positive",
    direction == "negative" & significant == "significant" ~ "Negative",
    TRUE ~ NA)) %>%
  mutate(direction = factor(direction, levels = c("Positive", "Low certainty", "Negative"))) %>%
  mutate(term = if_else(term == "Urban adapted species", "Urban-adapted species", term)) %>%
  filter(term != "Baseline")


# ---- Plot ----

# Create the plot
p <- ggplot(results, aes(y = reorder(term, estimate), x = estimate, color = direction)) +
  
  # Vertical line at zero
  geom_vline(xintercept = 0) +
  
  # Data geoms
  geom_point(size = 2.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0, lwd = 0.7) +
  
  # Scales
  scale_color_manual(values = c(
    "Negative" = "firebrick2", "Positive" = "royalblue3", "Low certainty" = "gray60")) +
  scale_x_continuous(#labels = scales::percent_format(),
    labels = scales::percent_format(scale=100),
    breaks = scales::pretty_breaks()) +
  
  # Labels
  labs(y = NULL, x = "Absolute effect on percent annual population growth rate") +
  
  # Theme
  theme_light(11) +
  theme(aspect.ratio = 7/16,
        legend.title = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text = element_text(color = 1),
        text = element_text(family = "Roboto", color = 1)); p

ggsave(p, filename = "Figure2.pdf", path = "/Users/gatesdupont/Desktop", 
       width = 7.5, height = 2.5,
       device = cairo_pdf)
