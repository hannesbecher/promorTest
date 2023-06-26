# Following promor tutrial for experiment withough technical reps: https://caranathunge.github.io/promor/articles/promor_no_techreps.html


setwd("~/git_repos/promorTest/")

#install.packages("promor")
library(promor)

# Create a raw_df object with default settings.
raw <- create_df(
  prot_groups = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
  exp_design = "https://raw.githubusercontent.com/caranathunge/promor_example_data/main/ed1.txt"
)

#?create_df

# downloads folder is in .gitignore
download.file("https://raw.githubusercontent.com/caranathunge/promor_example_data/main/pg1.txt",
              destfile = "downloads/pg1.txt")
# 1 Protein IDs     Majority protein IDs    Peptide counts (all)    Peptide counts (razor+unique)   Peptide counts (unique) Protein names   Gene names      Fasta headers   Number of proteins      Peptides        Razor >
# 2 A0AV96;B7Z8Z7;A0AV96-2;D6R9D6;D6RBS9;D6REZ6;D6RA49;D6RFL5;D6R9M7;D6RBP6;D6RCT1;Q5T0W7;B3KWU8;B4DZ27;Q9NQ94-5;Q8TBY0;Q9NQ94-6;Q9NQ94-3;C9JGD3;Q9NQ94-2;F8W9F8;Q9NQ94-4;Q9NQ94    A0AV96;B7Z8Z7;A0AV96-2;D6R9D6   13;12;>
# 3 A0AVT1;A0AVT1-2;H0Y8S8;A0AVT1-4;A0AVT1-3        A0AVT1;A0AVT1-2 14;11;4;3;3     14;11;4;3;3     14;11;4;3;3     Ubiquitin-like modifier-activating enzyme 6     UBA6    >sp|A0AVT1|UBA6_HUMAN Ubiquitin-like modifier->
# 4 H7BXI1;A0FGR8-6;A0FGR8-2;A0FGR8;C9JGI7;A0FGR8-4;F2Z3K9;A0FGR8-5;A6NFV7  H7BXI1;A0FGR8-6;A0FGR8-2;A0FGR8;C9JGI7;A0FGR8-4 11;11;10;10;7;7;4;4;2   11;11;10;10;7;7;4;4;2   11;11;10;10;7;7;4;4;2   Extended synaptotagmin>
 

# example of standard input:
download.file("https://raw.githubusercontent.com/caranathunge/promor_example_data/main/st.txt",
              destfile = "downloads/st.txt")
# 1 Protein H1      H2      H3      L1      L2      L3
# 2 A0AV96;B7Z8Z7;A0AV96-2;D6R9D6   166790000       176680000       140720000       202980000       245050000       259580000
# 3 A0AVT1;A0AVT1-2 34679000        46278000        24750000        72394000        71990000        58362000
# 4 H7BXI1;A0FGR8-6;A0FGR8-2;A0FGR8;C9JGI7;A0FGR8-4 56390000        55263000        58800000        85986000        101540000       91829000


# Filter out proteins with high levels of missing data in either condition
raw_filtered <- filterbygroup_na(raw, set_na = 0.4)
#> 224 proteins with higher than 40% NAs in at least one group removed.

dim(raw_filtered)
#> [1] 4360    6


# Visualize missing data in a subset of proteins.
heatmap_na(raw_filtered, palette = "mako")
?heatmap_na

# Order proteins by mean intensity.
heatmap_na(raw_filtered, reorder_y = TRUE, palette = "mako")


# Visualize missing data in a subset of proteins.
heatmap_na(raw_filtered, protein_range = 40:70, label_proteins = TRUE, palette = "mako")


# Impute missing data with minProb method. Don't forget to fix the random seed for reproducibility.
imp_df_mp <- impute_na(raw_filtered, seed = 327)


#heatmap_na(imp_df_mp)

# Visualize the imputed data with sample-wise density plots.
impute_plot(original = raw_filtered, imputed = imp_df_mp, n_row = 3, n_col = 3, palette = "mako")


# by sample
# Visualize the imputed data with sample-wise density plots.
impute_plot(original = raw_filtered, imputed = imp_df_mp, global = FALSE, n_row = 3, n_col = 3, palette = "mako")



# MAxLFQ data are normalised already. So, this should not be needed;
norm_df <- normalize_data(imp_df_mp)

norm_plot(original = imp_df_mp, normalized = norm_df, palette = "mako")

norm_plot(original = imp_df_mp, normalized = norm_df, type = "density", palette = "mako")


# differential expression
fit_df <- find_dep(imp_df_mp)
#> 1294 siginificantly differentially expressed proteins found.

#fit_df <- find_dep(imp_df_mp, file_path = getwd(), save_tophits = TRUE, n_top = 10)

# visualisation
volcano_plot(fit_df,
             text_size = 5,
             palette = "mako"
)

heatmap_de(fit_df, imp_df_mp, palette = "mako")
heatmap_de(fit_df, imp_df_mp, palette = "mako", n_top = 200)


# Seems to allow for simple designs only?
