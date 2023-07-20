Here we will make 2 NMDS graphs. The first one is just a quick one to check clusturing while the second is a bit more cleaned up. 

```
setwd("~/Desktop/Agbiome/Kerbel/02_ITS/R_analysis/Kerbel-ITS/")

# Load the needed libraries
library(mctoolsr) # microbial analysis
library(devtools)
library(vegan) # for stats
library(tidyverse) # ggplot2 + dplyr
library(ggpubr) # graph order


#The taxon table itself: "data_loaded"
tax_table_fp = "ASV_Table_ITS.txt"
# The metadata: "map_loaded"
map_fp = "ITS_Metadata_R.txt"
# Merge your mapping file and your OTU table
input = load_taxa_table(tax_table_fp, map_fp) #311 samples loaded
# for merging tables later in ggplot later
map = input$map_loaded %>% 
  rownames_to_column("Sample_ID")

# Rarefying sample reads
set.seed(1990)
input_rar = single_rarefy(input, 2500) #282 Samples remaining

# NMDS MCtoolsR way for a quick check
dm = calc_dm(input_rar$data_loaded)
ord = calc_ordination(dm, 'nmds')
plot_ordination(input_rar, ord, 'Treatment', hulls = TRUE)
```

![Screen Shot 2023-07-19 at 5 10 28 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/fec78b09-7e3c-4906-98fe-dd36dceb7242)

We can pull the ordination points (stored in "ord" and pipe them into ggplot in order to make this look better.

```
# using ggplot to make this look better
# Add MDS1 and MDS2 points to mapping file for ggplot
ord_map<- ord %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(map)
# Plot
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(shape = Depth, fill = Treatment, color = Treatment)) +
  scale_fill_manual(values = c("#CA054D", "#3B1C32", "#6FA37F")) +
  scale_color_manual(values = c("#CA054D", "#3B1C32", "#6FA37F")) +
  scale_shape_manual(values = c(16, 17, 16)) +
  stat_ellipse(aes(color = Treatment), level = 0.95)+
  ggtitle("Kerbel Bacterial full site") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```

![Screen Shot 2023-07-19 at 5 13 01 PM](https://github.com/LadyGrant/Kerbel/assets/95941680/e4f3f884-1d94-4f81-9a6f-d1b87a548576)

Now we need to run a PERMANOVA using adonis2

```
# run a PERMANOVA on your data (you need the DM to run a permanova)
permanova <- adonis2(dm ~ Treatment * Time * Depth * Block * Direction,
  data = input_rar$map_loaded, permutations= 999)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dm ~ Treatment * Time * Depth * Block * Direction, data = input_rar$map_loaded, permutations = 999)
                                      Df SumOfSqs      R2      F Pr(>F)
Treatment                              2    1.877 0.01530 2.2648  0.001
Time                                   2    1.085 0.00885 1.3096  0.001
Depth                                  2    1.973 0.01609 2.3812  0.001
Block                                  1    0.567 0.00462 1.3687  0.002
Direction                              1    0.850 0.00693 2.0504  0.001
Treatment:Time                         4    1.886 0.01538 1.1377  0.001
Treatment:Depth                        4    2.321 0.01893 1.4006  0.001
Time:Depth                             3    1.358 0.01107 1.0925  0.010
Treatment:Block                        2    1.096 0.00894 1.3231  0.001
Time:Block                             2    0.920 0.00750 1.1104  0.014
Depth:Block                            1    0.444 0.00362 1.0705  0.099
Treatment:Direction                    2    0.968 0.00790 1.1687  0.005
Time:Direction                         2    0.816 0.00665 0.9842  0.602
Depth:Direction                        1    0.504 0.00411 1.2164  0.004
Block:Direction                        1    0.581 0.00474 1.4025  0.001
Treatment:Time:Depth                   5    2.331 0.01901 1.1251  0.001
Treatment:Time:Block                   4    1.740 0.01419 1.0497  0.046
Treatment:Depth:Block                  2    0.888 0.00724 1.0718  0.043
Time:Depth:Block                       2    0.903 0.00737 1.0901  0.025
Treatment:Time:Direction               4    1.689 0.01377 1.0192  0.222
Treatment:Depth:Direction              2    0.862 0.00703 1.0404  0.143
Time:Depth:Direction                   2    0.856 0.00698 1.0329  0.192
Treatment:Block:Direction              2    1.327 0.01082 1.6018  0.001
Time:Block:Direction                   2    0.933 0.00761 1.1258  0.010
Depth:Block:Direction                  1    0.437 0.00356 1.0546  0.134
Treatment:Time:Depth:Block             4    1.712 0.01396 1.0329  0.124
Treatment:Time:Depth:Direction         4    1.644 0.01341 0.9919  0.597
Treatment:Time:Block:Direction         4    1.762 0.01437 1.0630  0.024
Treatment:Depth:Block:Direction        2    0.854 0.00697 1.0310  0.226
Time:Depth:Block:Direction             2    0.827 0.00674 0.9977  0.494
Treatment:Time:Depth:Block:Direction   3    1.264 0.01031 1.0172  0.268
Residual                             206   85.353 0.69603              
Total                                281  122.628 1.00000
```
Below, we changed some variable colors to view some patterns better. These are all the same NMDS!! just colored by different patterns.

```
# Colored by depth
# Deep: "#3395FF", Rhizosphere: "#44097A", Surface: "#FF33FB"
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(shape = Depth, fill = Depth, color = Depth)) +
  scale_fill_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_color_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_shape_manual(values = c(16, 17, 16)) +
  stat_ellipse(aes(color = Depth), level = 0.95)+
  ggtitle("Kerbel Fungal Full Site Depth") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```
![Screen Shot 2023-07-20 at 11 22 48 AM](https://github.com/LadyGrant/Kerbel/assets/95941680/9fa38129-4d96-451a-aa92-7a9530b55696)

```
# Colored by time
# July: "#3395FF", May: "#44097A", September: "#FF33FB"
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(shape = Depth, fill = Time, color = Time)) +
  scale_fill_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_color_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_shape_manual(values = c(16, 17, 16)) +
  stat_ellipse(aes(color = Time), level = 0.95)+
  ggtitle("Kerbel Fungal Full Site Time") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```
![Screen Shot 2023-07-20 at 11 24 39 AM](https://github.com/LadyGrant/Kerbel/assets/95941680/f6b8849f-1f00-459c-a1f8-89511be7fe66)

```
# colored by direction
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(fill = Direction, color = Direction, shape = Depth)) +
  scale_fill_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_color_manual(values = c("#3395FF","#44097A", "#FF33FB")) +
  scale_shape_manual(values = c(16, 17, 16)) +
  stat_ellipse(aes(color = Direction), level = 0.95)+
  ggtitle("Kerbel Fungal Full Site Direction") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```
![Screen Shot 2023-07-20 at 11 26 25 AM](https://github.com/LadyGrant/Kerbel/assets/95941680/97ffe287-0be4-435d-94db-f2bd24c79c17)

This project has a lot a variables, so lets's filter some out to see if we can get a better grasp at what is going on. Here, we are using the mctoolsr way of filtering, however a vegan or tidyverse filter will work as well. 

```
# Filtering down to surface, block 1, north
input_S = filter_data(input_rar, 'Depth', keep_vals = "S") # 135 samples left
input_S_B1 = filter_data(input_S, "Block", keep_vals = "B1") #67 samples left
input_S_B1_N = filter_data(input_S_B1, "Direction", keep_vals = "N") # 34 samples

# Set a new seed and calcualte new distance matrix
set.seed(1234)
dm = calc_dm(input_S_B1_N$data_loaded)
ord = calc_ordination(dm, 'nmds')

# Add MDS1 and MDS2 points to mapping file for ggplot
ord_map <- ord %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(map)
# Colored by Treatment
#  CT: #CA054D, MT: #3B1C32", ST: #6FA37F
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(color = Treatment)) +
  scale_color_manual(values = c("#CA054D", "#3B1C32", "#6FA37F")) +
  scale_shape_manual(values = c(16, 16, 16)) +
  stat_ellipse(aes(color = Treatment), level = 0.95)+
  ggtitle("Kerbel Fungal Surface") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```

![Screen Shot 2023-07-20 at 11 38 09 AM](https://github.com/LadyGrant/Kerbel/assets/95941680/0511ad9f-cb47-4a0e-9680-4aab017d15cb)

Now let's run a PERMANOVA for the two varibles left, treatment and time.

```
# run a PERMANOVA on your data
permanova <- adonis2(dm ~ Treatment * Time,
  data = input_S_B1_N$map_loaded, permutations= 999)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dm ~ Treatment * Time, data = input_S_B1_N$map_loaded, permutations = 999)
               Df SumOfSqs      R2      F Pr(>F)    
Treatment       2   1.3404 0.09425 1.6273  0.001 ***
Time            2   0.8854 0.06226 1.0749  0.074 .  
Treatment:Time  4   1.6999 0.11953 1.0319  0.169    
Residual       25  10.2961 0.72396                  
Total          33  14.2218 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
Now let's check the rhizosphere.

```
# Filtering down to surface, block 1, north
Rhiz = filter_data(input_rar, 'Depth', keep_vals = "rhizo") # 16 samples left

# Set a new seed and calcualte new distance matrix
set.seed(444)
dm = calc_dm(Rhiz$data_loaded)
ord = calc_ordination(dm, 'nmds')

# Add MDS1 and MDS2 points to mapping file for ggplot
ord_map <- ord %>% 
  rownames_to_column("Sample_ID") %>%
  left_join(map)

# Colored by Treatment
#  CT: #CA054D, MT: #3B1C32", ST: #6FA37F
NMDS_all <- ord_map %>% 
  ggplot(aes(MDS1, MDS2)) +
  geom_point(aes(color = Treatment)) +
  scale_color_manual(values = c("#CA054D", "#3B1C32", "#6FA37F")) +
  scale_shape_manual(values = c(16, 16, 16)) +
  stat_ellipse(aes(color = Treatment), level = 0.95)+
  ggtitle("Kerbel Fungal Rhizosphere") + theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))
NMDS_all
```
![Screen Shot 2023-07-20 at 11 46 47 AM](https://github.com/LadyGrant/Kerbel/assets/95941680/5355e1e6-7733-46bb-9b9c-e32d723ebc19)

Run the PERMANOVA on Rhizosphere and time.

```
# run a PERMANOVA on your data
permanova <- adonis2(dm ~ Treatment * Time,
  data = Rhiz$map_loaded, permutations= 999)
```
```
Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 999

adonis2(formula = dm ~ Treatment * Time, data = Rhiz$map_loaded, permutations = 999)
               Df SumOfSqs      R2      F Pr(>F)    
Treatment       2   1.0409 0.15888 1.2511  0.001 ***
Time            1   0.4485 0.06846 1.0782  0.134    
Treatment:Time  1   0.4862 0.07421 1.1688  0.016 *  
Residual       11   4.5758 0.69845                  
Total          15   6.5514 1.00000                  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```





























