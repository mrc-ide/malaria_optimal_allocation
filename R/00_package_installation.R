# Install required packages

# Mathematical model
devtools::install_github("mrc-ide/odin")
devtools::install_github("mrc-ide/deterministic-malaria-model@feat/vivax")

# Optimisation
devtools::install_version("akima", version = "0.6.2.2", repos = "http://cran.us.r-project.org")
install.packages("GenSA")

# Spatial data 
install.packages("sf")
install.packages("raster") 
install.packages("terra")
install.packages("exactextractr")
install.packages("rgdal")
install.packages("spatial")

# Data manipulation and plotting
install.packages("tidyverse")
install.packages("here")
install.packages("fuzzyjoin")
install.packages("ggplot2")
install.packages("patchwork")
install.packages("ggrepel")
install.packages("gridExtra")
install.packages("plotly")
install.packages("ggpubr")
