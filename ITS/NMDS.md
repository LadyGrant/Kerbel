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


