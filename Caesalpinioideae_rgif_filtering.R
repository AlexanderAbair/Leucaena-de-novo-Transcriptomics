library(CoordinateCleaner)
library(rnaturalearthdata)
library(dplyr)
library(scrubr)
library(raster)
library(sp)
library(countrycode)
library(ggplot2)

#world <- getData("countries")
#names(world)[names(world) == "ISO"] <- "iso_a3_eh"

Caesalpinioid_occs <- read.csv('Caesalpinioideae_gbif_occurrences_unfiltered.csv')
Caesalpinioid_occs <- filter(Caesalpinioid_occs,
                             decimalLongitude != 0 | decimalLatitude != 0,
                             basisOfRecord == "PRESERVED_SPECIMEN",
                             order == "Fabales")

Caesalpinioideae_gbif_occurrences_filtered <- cc_cap(Caesalpinioid_occs, lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_cen(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_dupl(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_zero(lat = "decimalLatitude", lon = "decimalLongitude") %>% 
  cc_equ(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_gbif(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_inst(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  cc_sea(lat = "decimalLatitude", lon = "decimalLongitude") %>%
  filter(decimalLatitude >= -60) %>% 
  filter(decimalLatitude <= 70)
  #cc_coun(lat = "decimalLatitude", lon = "decimalLongitude", iso3 =  "countryCode") # Crashes the code 

write.csv(Caesalpinioideae_gbif_occurrences_filtered, 
          "Caesalpinioideae_gbif_occurrences_filtered.csv")

# Plotting occurrences on a world map

world_map <- borders("world", colour="gray50", fill="gray50")
Caesalpinioid_world_map <- ggplot() + coord_fixed() + world_map +
  geom_point(data = Caesalpinioideae_gbif_occurrences_filtered, 
             aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", 
             size = 0.5) +
  theme_bw()
ggsave("Caesalpinioideae_world_map.png")  
