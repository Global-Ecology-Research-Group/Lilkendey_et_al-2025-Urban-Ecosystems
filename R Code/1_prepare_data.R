# Define study areas as parks in urban areas of Florida, filter iNaturalist data to these urban areas, and
# then use the iNaturalist API to obtain list of non-native angiosperm species

library(dplyr)
library(readr)
library(sf)
library(rnaturalearth)
library(leaflet)
library(httr)
library(jsonlite)
library(dplyr)
library(purrr)
sf_use_s2(FALSE)

# Get extract urban areas in Florida from the ParkServe data -----------------------------------

# The code in this section cannot be run because the ParkServe_Parks.shp is too large to store on the server
# However, the outputted parks_in_urban_areas.shp file is provided.

# get polygon of Florida
us_states <- ne_states(country = "united states of america", returnclass = "sf")
florida <- us_states %>% dplyr::filter(name == "Florida")

# get the polygons of urban areas that we want to extract iNaturalist data
urban_areas <- st_read("Data/Urban_Areas/Urban_areas.shp")
parks <- st_read("Data/ParkServe/ParkServe_Parks.shp")

# crop files to Florida
# Transform to Florida's CRS
urban_areas_fl <- st_transform(urban_areas, st_crs(florida))
parks_fl <- st_transform(parks, st_crs(florida))

# make geometries valid
parks_fl_crop <- st_make_valid(parks_fl)

urban_areas_crop <- st_intersection(urban_areas_fl, florida)
parks_crop <- st_intersection(parks_fl_crop, florida)
parks_in_urban <- st_intersection(parks_fl_crop, urban_areas_crop)

# save the parks in urban areas shapefile
st_write(parks_in_urban, "Data/ParkServe/parks_in_urban_areas.shp")
parks_in_urban <- st_read("Data/ParkServe/parks_in_urban_areas.shp")


## Plot urban areas in Florida ---------------------------------------------

# plot the data
plot(st_geometry(florida), col = "lightgray")
plot(st_geometry(parks_in_urban), col = "darkgreen", add = TRUE)

# use leaflet to further investigate the urban areas
leaflet() %>%
  addProviderTiles("CartoDB.Positron") %>%
  addPolygons(data = florida, fillColor = "lightgray", weight = 1, color = "black") %>%
  addPolygons(data = parks_in_urban, fillColor = "darkgreen", fillOpacity = 0.6, color = "green", weight = 0.5)


# Filter iNaturalist data to urban areas -----------------------------------------

# Raw iNaturalist data is not provided in this repository because the file is too large to store.
# This data was downloaded on June 26, 2025 using the iNaturalist data export tool.
# We do, however, provide the cleaned data file, named florida_pollinators_and_angiosperms_taxon.csv.

# get iNaturalist data - this code will not work since raw data is not provided
data_files <- list.files()
data <- bind_rows(lapply(data_files, read_csv))

# filter to only include pollinators and angiosperms, and research grade data

data_filtered <- data %>%
  filter(taxon_superfamily_name == "Apoidea" |
           taxon_family_name == "Bombyliidae" |
           taxon_subfamily_name == "Cetoniinae" |
           taxon_order_name == "Lepidoptera" |
           taxon_subfamily_name == "Lepturinae" |
           taxon_subphylum_name == "Angiospermae")

nrow(data_filtered)

data_rg <- data_filtered %>%
  filter(quality_grade=="research")

data_clean <- data_rg %>%
  mutate(group=ifelse(taxon_subphylum_name == "Angiospermae", 
                      "Angiosperms", "Pollinators")) %>%
  dplyr::select(id, created_at, observed_on, latitude, longitude, taxon_species_name, group,
                starts_with("taxon_")) %>%
  rename(scientific_name=taxon_species_name)

saveRDS(data_clean, "Data/Pollinator_and_Angiosperm_Data/florida_pollinators_and_angiosperms_taxon.RDS")

# Obtain list of non-native angiosperm species from the iNaturalist API -----------------------------------------

# set up the API call
base_url <- "https://api.inaturalist.org/v1/observations/species_counts"
per_page <- 200
page <- 1
all_pages <- list()

# set up loop to extract all data
repeat {
  message("Fetching page ", page)
  
  # Make GET request with page number
  response <- GET(base_url, query = list(
    taxon_id = 47125,
    place_id = 21,
    introduced = "true",
    quality_grade = "research",
    per_page = per_page,
    page = page
  ))
  
  stop_for_status(response)
  
  data_json <- content(response, as = "text", encoding = "UTF-8")
  data_parsed <- fromJSON(data_json, flatten = TRUE)
  
  results_df <- as_tibble(data_parsed$results)
  
  # Stop if no results
  if (nrow(results_df) == 0) {
    message("No more results, stopping.")
    break
  }
  
  # Append this page's results
  all_pages[[page]] <- results_df
  
  # Calculate total pages from total_results if available
  total_results <- data_parsed$total_results
  total_pages <- ceiling(total_results / per_page)
  
  if (page >= total_pages) {
    message("Reached last page: ", page)
    break
  }
  
  page <- page + 1
}

# Combine all pages into one tibble
combined_df <- bind_rows(all_pages)

print(combined_df)

combined <- combined_df %>% as.data.frame()

write_csv(combined, "Data/Pollinator_and_Angiosperm_Data/all_nonnative_sp.csv")
