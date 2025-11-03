This repository contains all the code and files used in the analyses reported in "Proportion of native plants is a key predictor of pollinator richness in urban greenspaces" article, which is published in Urban Ecosystems.

# Data Folder

This folder contains all data used in the analyses. The only datasets not included, due to storage limitations, are the raw iNaturalist data for the state of Florida and the full ParkServe dataset. The raw iNaturalist data was downloaded using the iNaturalist data export tool on June 26, 2025. However, we provide a filtered version of this dataset, which includes only the taxonomic groups of interest and observations classified as “Research Grade.” We additionally provide the filtered ParkServe greenspace boundaries to only include those located within urban areas in Florida. The following is a description of data present in this repository:

**Dynamic_World_Habitat/dynamic_world_habitat_fl.csv** – A CSV file listing all urban parks in our study area and the percentage of each *Dynamic World* habitat class within those greenspaces. These percentages were calculated using Google Earth Engine. Specifically, we applied a pixel-based histogram reducer to count the number of 100 square-meter land cover pixels belonging to each habitat type and then calculated the percentage of each type relative to the total number of pixels within each polygon. **Dataset citation:** Brown, C. F., Brumby, S. P., Guzder-Williams, B., *et al.* (2022). *Dynamic World: Near real-time global 10 m land use land cover mapping.* *Scientific Data, 9*(251). <https://doi.org/10.1038/s41597-022-01307-4>

**ParkServe/parks_in_urban_areas.shp** – Filtered ParkServe greenspace boundaries to only include those located within urban areas in Florida. The full ParkServe dataset is too large to store in this repository, but it can downloaded online. **Dataset citation:** Trust for Public Land (2025) ParkServe [database]. Land and People Lab. <https://tpl.org/parkserve>

**Pollinator_and_Angiosperm_Data/all_nonnative_sp.csv** – A CSV file output from `1_prepare_data.R` which contains taxonomic information on all introduced angiosperm species in Florida.

**Pollinator_and_Angiosperm_Data/florida_pollinators_and_angiosperms_taxon.RDS** – An RDS file containing all Research Grade iNaturalist observations of pollinators (Superfamily Apoidea, Family Bombyliidae, Subfamily Cetoniinae, Order Lepidoptera, and Subfamily Lepturinae) and angiosperm species in the state of Florida through June 26, 2025. This data was obtained using the iNaturalist data export tool, which was accessed on June 26, 2025.

**Pollinator_and_Angiosperm_Data/species_richness_wide_fl.csv** – A CSV file output from `2_data_analysis.R` which contains pollinator and angiosperm richness values and number of iNaturalist observations per park for the 129 parks used in this study. This file was used to show trends in richness values and observation intensity in the `4_map_of_study_area_and_data_availablity.R` script.

**Urban_Areas/Urban_areas.shp** – Shapefile containing boundaries of urban areas in Florida. This data was obtained from the U.S. Census Bureau. **Dataset Citation:** United States Census Bureau (2023) 2020 Census urban area TIGER/Line shapefiles. Retrieved September 2024. <https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html>

# Figures Folder

This folder contains all the main figures and supplemental figures presented in the article.

# R Folder

This folder contains 4 R scripts which can be used to repeat the results presented in the article and supplemental material. Below is a description of each script:

**1_prepare_data** – This script defines urban greenspaces in Florida using the `ParkServe_Parks.shp` and `Urban_areas.shp` files, reads in raw iNaturalist data for the entire states, filters the data to include only the taxonomic groups of interest and "Research Grade" observations, and generates a list of non-native angiosperm species in Florida. The non-native species list was obtained on July 17, 2025. If this script is run in the future, it may produce a slightly different list as species are added to or removed from the list of introduced plants in Florida by the iNaturalist community.

**2_data_analysis** – This script contains all the code used to analyze pollinator and angiosperm data from urban greenspaces in Florida to generate the main results presented in the article.

**3_supplemental_analysis_by_subgroup** – This script contains all the code used to repeat the main analysis reported in the article for four subgroups - Apoidea, Lepidoptera, butterflies (Family Hesperiidae, Papilionidae, Pieridae, Lycaenidae, Riodinidae, and Nymphalidae), and moths (Lepidoptera and not family Hesperiidae, Papilionidae, Pieridae, Lycaenidae, Riodinidae, and Nymphalidae).

**4_map_of_study_area_and_data_availability** – This script is used to map the `species_richness_wide_fl.csv` data to show trends in richness values and observation intensity across urban greenspaces in Florida.
