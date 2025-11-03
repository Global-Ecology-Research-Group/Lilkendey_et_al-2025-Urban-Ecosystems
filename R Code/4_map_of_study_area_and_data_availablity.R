# Script to make study area map

library(tidyverse)
library(sf)
library(rnaturalearth)
library(osmdata)
library(tmap)
library(biscale)
library(cowplot)
library(leaflet)
library(htmlwidgets)

# determine the 129 park study areas
data <- read_csv("Data/Angiosperm_vs_Pollinator_Diversity/species_richness_wide_fl.csv")

# read in park_cov
parks <- st_read("Data/ParkServe/parks_in_urban_areas.shp")

# filter parks to the parks that we used in the study
parks <- parks[parks$ParkID %in% data$ParkID,]

parks <- parks %>%
  group_by(ParkID) %>%
  summarise(NAME20=first(NAME20))

# read in urban areas
urban_areas <- st_read("Data/Urban_Areas/Urban_areas.shp")
urban_areas <- st_transform(urban_areas, st_crs(florida))
urban_areas_fl <- st_intersection(urban_areas, florida)
urban_areas_fl <- urban_areas_fl %>%
  dplyr::select(NAME20)

# get count of parks in urban areas
parks_urban_areas <- parks %>%
  group_by(NAME20) %>%
  summarise(count=n()) %>%
  as.data.frame()

urban_areas_fl <- left_join(urban_areas_fl, parks_urban_areas, by=c("NAME20"))

# check that everything is adding up
sum(urban_areas_fl$count, na.rm=TRUE)

# plot florida
# Get US states as an sf object
us_states <- ne_states(country = "United States of America", returnclass = "sf")

# filter for Florida
florida <- us_states[us_states$name == "Florida", ]

# transform to border CRS
florida <- st_transform(florida, st_crs(border))

# Map urban areas with number of urban parks ------------------------------------------------

# remove urban areas that do not have parks that meet our criteria
urban_areas_fl_comp <- urban_areas_fl[complete.cases(urban_areas_fl$count),]

# define color pallete
purple_palette <- c("#e7d4e8", "#6f4e87", "#3a0f42")

ggplot() +
  geom_sf(data = florida, fill = "grey95", color = "grey20", lwd = 1) +
  geom_sf(data = urban_areas_fl_comp, aes(fill = count)) +
  scale_fill_gradientn(colors = purple_palette, name = "Number of Parks",
                       guide = guide_colorbar(
                         barwidth = 4,
                         barheight = 15)) +
  theme_void() +
  theme(
    legend.title = element_text(size = 34),
    legend.text = element_text(size = 32),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )


# Map parks just in South Florida to use as an example ------------------------------------------------

# South Florida border
south_fl <- urban_areas_fl_comp[urban_areas_fl_comp$NAME20=="Miami--Fort Lauderdale, FL",]

# data for South Florida parks
parks_sf <- parks %>%
  filter(NAME20=="Miami--Fort Lauderdale, FL",)
parks_sf <- left_join(parks_sf, data, by=c("ParkID"))

# create a bi class for bivariate mapping
parks_data_bi <- bi_class(parks_sf, x = Pollinators, y = `Angiosperms all`, style = "quantile", dim = 3)
parks_data_bi$bi_class <- factor(parks_data_bi$bi_class)

# let's see the hex values
bi_pal(pal = "PurpleOr", dim = 3, preview= FALSE)

map <- ggplot() +
  geom_sf(data = south_fl, fill = "grey95", color = "grey20", lwd = 1) +
  geom_sf(data = parks_data_bi,aes(fill = bi_class, color = bi_class), size = 0.5) +
  bi_scale_fill(pal = "PurpleOr", dim = 3) +
  bi_scale_color(pal = "PurpleOr", dim = 3) +  
  theme_void() +
  theme(
    legend.position = "none", 
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

legend <- bi_legend(pal = "PurpleOr",
                    dim = 3,
                    xlab = "More Pollinators",
                    ylab = "More Angiosperms",
                    size = 20)

# combine map and legend
cowplot::plot_grid(
  map, legend,
  rel_widths = c(1, 0.3),
  nrow = 1
)

# let's see the hex values
bi_pal(pal = "PurpleOr", dim = 3, preview= FALSE)

