# Here we are going to look at different diversity metrics
## We are still using MCtoolsR here but again, any of these can be used in vegan.

First, lets filter down to the subset of samples we are focusing on and make a heat map to visualize who is the most abundant within treatments. At this point we have decided to drop MT and stick to the more extreame treatments, CT and ST. We will check each taxonomic level to see where any significance may appear.

## Heatmaps
```
# Filtering down to surface, block 1, north
input_S = filter_data(input_rar, 'Depth', keep_vals = "Surface") # 135 samples left
input_S_B1 = filter_data(input_S, "Block", keep_vals = "B1") #67 samples left
input_S_B1_N = filter_data(input_S_B1, "Direction", keep_vals = "North") # 34 samples
input_S_B1_N_CTST = filter_data(input_S_B1_N, "Treatment", filter_vals = "MT") # 24 samples left


# summerize taxonomy by level
# levels: 1=K, 2=P, 3=C, 4=O, 5=F, 6=G, 7=S
tax_sum_families = summarize_taxonomy(input_S_B1_N_CTST, level = 2, report_higher_tax = FALSE)
# Plot heatmap
# Red is a higher abundance, Blue is the lowr abundance.
Heatmap <- plot_ts_heatmap(tax_sum_families, input_S_B1_N_CTST$map_loaded, 
  0.01, 'Treatment') + ylab("Treatment") + ggtitle("Heatmap ITS Phylum")
Heatmap
```
Heat maps for surface soils:
![Screen Shot 2023-07-20 at 12 19 28 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/01321bc6-6cd4-475b-8e31-cf1ee9df9ba0)
![Screen Shot 2023-07-20 at 12 21 00 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/52f62db2-729c-4b9c-b63e-84dfd2320f84)
![Screen Shot 2023-07-20 at 12 22 58 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/091febb1-63bf-4e66-a3aa-ec609da9bca2)
![Screen Shot 2023-07-20 at 12 23 46 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/8d312eb6-77a8-4944-b953-17126fe6fac5)
![Screen Shot 2023-07-20 at 12 24 57 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/4bf9055d-80da-402a-9f44-34d9300e6bff)
![Screen Shot 2023-07-20 at 12 25 49 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/e0a0ee87-bc6c-40bd-9516-9de68c793491)
Heat maps for rhizosphere soils:
![Screen Shot 2023-07-20 at 2 10 15 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/0d6375cb-9a62-4dd3-95e1-50db7d424c24)
![Screen Shot 2023-07-20 at 2 11 01 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/7c4a669c-2382-4cab-861a-c6920f0e23af)
![Screen Shot 2023-07-20 at 2 11 34 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/28441175-6a6c-4197-a093-42fdc99c85d2)
![Screen Shot 2023-07-20 at 2 12 07 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/454ed50b-ae87-436c-91a1-ec24707a112e)
![Screen Shot 2023-07-20 at 2 12 38 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/030e5a52-867b-40a8-9537-0bdb802b044c)
![Screen Shot 2023-07-20 at 2 13 14 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/da7bd127-86c4-4ed4-8b40-fdafa6b880a2)

## Stacked taxa bars
Now let's look at the top 10 taxa in a stacked taxa bar plot. Note that you need to change the "tax_sum_families" level for each new graph you generate. Also note that you need to change your "input..." to surface or rhizosphere. 
```
# levels: 1=K, 2=P, 3=C, 4=O, 5=F, 6=G, 7=S
tax_sum_families = summarize_taxonomy(input_S_B1_N_CTST, level = 2, report_higher_tax = FALSE)
# Plot taxa bars
taxa_bar_P <- plot_taxa_bars(tax_sum_families, input_S_B1_N_CTST$map_loaded, "Treatment",
  num_taxa = 10, data_only = FALSE) +
  scale_fill_brewer(palette="Spectral") + theme_classic() +
  ggtitle("ITS Phylum")
taxa_bar_P
```
Below is the code to arrange your graphs.
```
ggarrange(taxa_bar_P, taxa_bar_C, taxa_bar_O, taxa_bar_F, taxa_bar_G, taxa_bar_S,
          ncol = 3, nrow = 2)
```
Stacked taxa bars for surface soils:
![Screen Shot 2023-07-20 at 2 00 18 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/12516b74-74c6-4949-938c-a1d18428d196)
Stacked taxa bars for rhizosphere:
![Screen Shot 2023-07-20 at 2 21 18 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/4b18ce87-46b8-42b4-b44a-792e81b82555)

## Richness
Below are some richness plots that take the richness caculations from MCtoolsR and add some extra details using ggplot2 add-ons.








