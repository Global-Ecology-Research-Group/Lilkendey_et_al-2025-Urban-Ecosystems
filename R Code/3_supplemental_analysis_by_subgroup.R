# In this script, we use iNaturalist data to assess the relationship between pollinators, 
# native and non-native angiosperms, and several environmental covariates in urban greenspaces by subgroups.

# This analysis was conducted because some taxonomic groups are more likely to reach research grade status than others
# and some have more consistency with identification than others

### Read in packages 
library(tidyverse)
library(sf)
library(mgcv)
library(piecewiseSEM)
library(multcompView)
library(visreg)
library(semPlot)
library(MASS)
library(ggrepel)
library(janitor)
library(exactextractr)
library(raster)
library(broom)
library(patchwork)
library(ggtext)
library(gratia)
set.seed(129) # for reproducibility

# Read in data and clean --------------------------------------------------

# read in relevant shapefiles - urban areas and parks in urban areas
urban_areas <- st_read("tmlilkendey/Data/Urban_Areas/Urban_areas.shp")
parks_in_urban <- st_read("Data/ParkServe/parks_in_urban_areas.shp")

# read in florida data of pollinators and angiosperms
data_clean <- readRDS("Data/Angiosperm_vs_Pollinator_Diversity/florida_pollinators_and_angiosperms_taxon.RDS")

# make the data a shapefile
data_sf <- st_as_sf(data_clean, coords = c("longitude", "latitude"), crs = 4326)

# select relevant data from the parks_in_urban shapefile
parks_in_urban_rel <- parks_in_urban %>%
  dplyr::select(ParkID, Park_Size_)

# add the parks in urban information to the iNaturalist data shapefile
data_fl_parks <- st_join(data_sf, parks_in_urban_rel)

# remove all data that not in a park
data_fl_parks_only <- data_fl_parks[complete.cases(data_fl_parks$ParkID),]

# we will use the urban areas as random variable later for modeling, so let's 
# add that information into the inaturalist data
data_fl_ua <- st_join(data_fl_parks_only, urban_areas %>% dplyr::select(NAME20))

# polish the data and filter so all observations are after 2014
parks_gbif_assigned <- data_fl_ua %>%
  mutate(
    longitude = st_coordinates(.)[, 1],
    latitude = st_coordinates(.)[, 2]
  ) %>% 
  as.data.frame() %>%
  rename(urban_area=NAME20)  %>% 
  mutate(observed_on = ymd(observed_on)) %>%
  filter(observed_on >= ymd("2014-01-01"))

# range of dates encompassed
summary(parks_gbif_assigned$observed_on)
sd(parks_gbif_assigned$observed_on)

# total number of parks in Florida urban areas
cat("Total Florida urban parks:", length(unique(parks_gbif_assigned$ParkID)))

# total number of observations in urban parks
cat("Total observations in urban parks:", nrow(parks_gbif_assigned))

# Get number of observations by group
obs_summary <- parks_gbif_assigned %>%
  as.data.frame() %>%
  filter(group %in% c("Pollinators", "Angiosperms")) %>%
  group_by(ParkID, group) %>%
  summarise(number_obs = n(), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = number_obs, values_fill = 0)

# Let's try filtering this data so that each group has at least 50 observations
included_parks <- obs_summary %>%
  filter(Angiosperms>=50 & Pollinators>=50) %>%
  rename("angiosperm_obs"=Angiosperms, "pollinator_obs"=Pollinators)
nrow(included_parks)
# this gives us 129 parks

# Filter the full data to these "included_parks"
cleaned_data <- parks_gbif_assigned %>%
  filter(ParkID %in% included_parks$ParkID,
         complete.cases(scientific_name))
nrow(cleaned_data)

# Prepare data for modeling ------------------------

# Quantify number of pollinator species per park
pollinators <- cleaned_data %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(ParkID) %>%
  summarise(number_of_species = n_distinct(scientific_name)) %>%
  mutate(group="Pollinators")

# Quantify number of native and non-native angiosperm species

# Read in the Broward county non-native angiosperm observation data from iNaturalist
non_native_angiosperm <- read.csv("Data/Pollinator_and_Angiosperm_Data/all_nonnative_sp.csv")

# Generate a list of non-native angiosperms in Broward county
non_natives <- non_native_angiosperm %>%
  dplyr::select(taxon.name) %>%
  distinct() %>% # get only unique species
  mutate(status="Non-native")

cleaned_data %>%
  as.data.frame() %>%
  dplyr::filter(group=="Angiosperms") %>%
  left_join(., non_natives, by=c("scientific_name"="taxon.name")) %>%
  replace_na(list(status="Native")) %>%
  group_by(status) %>%
  summarise(count=n(),
            percentage=n()/nrow(cleaned_data %>% dplyr::filter(group=="Angiosperms"))*100)

# Quantify number of native and non-native angiosperm species per park
angiosperms <- cleaned_data %>%
  dplyr::filter(group=="Angiosperms") %>%
  left_join(., non_natives, by=c("scientific_name"="taxon.name")) %>%
  replace_na(list(status="Native")) %>%
  group_by(ParkID, status) %>%
  summarise(number_of_species = n_distinct(scientific_name)) %>%
  mutate(group="Angiosperms")

# Quantify number of angiosperm species per park, regardless of indigeneity 
angiosperms_all <- angiosperms %>%
  group_by(ParkID) %>%
  summarize(number_of_species=sum(number_of_species)) %>%
  mutate(group="Angiosperms all")

# Create data set from the above data sets: pollinators, angiosperms, and angiosperms_all. For each park, 
# the data set includes the number of pollinator species, the number of angiosperms total, and the number of native
# versus non-native angiosperms. 
species_count <- pollinators %>%
  bind_rows(angiosperms_all) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Native")) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Non-native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Non-native"))


# Count the total number of pollinator and angiosperm observations per park
annotated_data <- cleaned_data %>%
  left_join(non_natives, by=c("scientific_name"="taxon.name")) %>%
  mutate(status = replace_na(status, "Native"))

# Now summarize by ParkID and group/status
observations_count <- annotated_data %>%
  mutate(group = as.character(group)) %>%
  group_by(ParkID) %>%
  summarise(
    pollinator_obs = sum(group == "Pollinators"),
    angiosperm_obs = sum(group == "Angiosperms"),
    native_angiosperm_obs = sum(group == "Angiosperms" & status == "Native"),
    nonnative_angiosperm_obs = sum(group == "Angiosperms" & status == "Non-native"),
    .groups = "drop"
  )

# Reformat species richness data (observation_count dataset) to wide format. This data set quantifies 
# the biodiversity (richness), and number of pollinator/angiosperm observations for the final
# set of greenspaces. Also, add the observation counts of pollinators and angiosperms to the data frame
species_richness_wide <- species_count %>%
  pivot_wider(names_from = "group", values_from = "number_of_species", 
              values_fill = list(number_of_species = 0)) %>%
  left_join(observations_count, by = "ParkID") 

# add urban area name
urban_area_names <- parks_gbif_assigned %>% 
  dplyr::select(ParkID, longitude, latitude, urban_area) %>%
  group_by(ParkID) %>%
  summarise(lon=first(longitude), lat=first(latitude), urban_area=first(urban_area))

species_richness_wide <- species_richness_wide %>%
  left_join(urban_area_names, by="ParkID")

# Analyze pollinator richness and angiosperm richness by subgroup ------------------------

# let's create subgroups
data_taxon <- cleaned_data %>%
  mutate(group_name = case_when(
    taxon_superfamily_name == "Apoidea" ~ "Apoidea",
    taxon_family_name == "Bombyliidae" ~ "Bombyliidae",
    taxon_subfamily_name == "Cetoniinae" ~ "Cetoniinae",
    taxon_order_name == "Lepidoptera" ~ "Lepidoptera",
    taxon_subfamily_name == "Lepturinae" ~ "Lepturinae",
    taxon_subphylum_name == "Angiospermae" ~ "Angiospermae",
    TRUE ~ NA_character_),
    butterfly_moth = ifelse(taxon_family_name %in% c("Hesperiidae", "Papilionidae",
                                                     "Pieridae", "Lycaenidae",
                                                     "Riodinidae", "Nymphalidae"), "butterfly", NA),
    butterfly_moth = ifelse(taxon_order_name == "Lepidoptera" &
                              !taxon_family_name %in% c("Hesperiidae", "Papilionidae",
                                                     "Pieridae", "Lycaenidae",
                                                     "Riodinidae", "Nymphalidae"), "moth", butterfly_moth))

# Quantify number of pollinator groups in the dataset
# distinct species in dataset
distinct_sp <- data_taxon %>%
                      dplyr::filter(group=="Pollinators")

# number of species per group
data_taxon %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(group_name) %>%
  summarise(number_of_species = n_distinct(scientific_name),
            prop = n_distinct(scientific_name)/n_distinct(distinct_sp$scientific_name)*100) %>%
  mutate(group="Pollinators")

# number of observations per group
data_taxon %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(group_name) %>%
  summarise(number_of_observations = n(),
            prop = n()/nrow(data_taxon %>%
                              dplyr::filter(group=="Pollinators"))*100) %>%
  mutate(group="Pollinators")

# since there are a lot of observations in the Lepidoptera group, let's also examine differences between
# butterflies versus moths
distinct_sp_lep <- data_taxon %>%
  dplyr::filter(group=="Pollinators",
                taxon_order_name == "Lepidoptera")

# number of species
data_taxon %>%
  dplyr::filter(group=="Pollinators",
                taxon_order_name == "Lepidoptera") %>%
  group_by(butterfly_moth) %>%
  summarise(number_of_species = n_distinct(scientific_name),
            prop = n_distinct(scientific_name)/n_distinct(distinct_sp$scientific_name)*100) 

# number of observations
data_taxon %>%
  dplyr::filter(group=="Pollinators",
                taxon_order_name == "Lepidoptera") %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(butterfly_moth) %>%
  summarise(number_of_observations = n(),
            prop = n()/nrow(data_taxon %>%
                              dplyr::filter(group=="Pollinators"))*100) %>%
  mutate(group="Pollinators")

# what are the most prominent families in Lepidoptera
(family_sums <- data_taxon %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(group_name, taxon_family_name) %>%
  summarise(number_of_species = n_distinct(scientific_name),
            number_of_observations = n()))

# what about subfaimiles?
(subfam_sums <- data_taxon %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(group_name, taxon_subfamily_name) %>%
  summarise(number_of_species = n_distinct(scientific_name),
            number_of_observations = n()))

# Quantify number of pollinator species per park by group
pollinators <- data_taxon %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(ParkID, group_name) %>%
  summarise(number_of_obs = n(),
            number_of_species=n_distinct(scientific_name)) %>%
  mutate(group="Pollinators")
pollinators

# do the same for leps
pollinators_lep <- data_taxon %>%
  dplyr::filter(group=="Pollinators",
                taxon_order_name == "Lepidoptera") %>%
  group_by(ParkID, butterfly_moth) %>%
  summarise(number_of_obs = n(),
            number_of_species=n_distinct(scientific_name))
pollinators_lep

# let's filter for more than 50 observations
pollinators_filtered <- pollinators %>%
  filter(number_of_obs >= 50)

pollinators_filtered %>% 
  group_by(group_name) %>%
  summarise(count=n())
# 12 parks meet the criteria for Apoidea and 113 meet the criteria for Lepidoptera

# repeat for butterflies and moths
pollinators_filtered_lep <- pollinators_lep %>%
  filter(number_of_obs >= 50)

pollinators_filtered_lep %>% 
  group_by(butterfly_moth) %>%
  summarise(count=n())
# 80 parks meet the criteria for butterflies and 28 meet the criteria for moths

# Quantify number of native and non-native angiosperm species per park
angiosperms <- data_taxon %>%
  dplyr::filter(group=="Angiosperms") %>%
  left_join(., non_natives, by=c("scientific_name"="taxon.name")) %>%
  replace_na(list(status="Native")) %>%
  group_by(ParkID, status) %>%
  summarise(number_of_species = n_distinct(scientific_name)) %>%
  mutate(group="Angiosperms")

summary(angiosperms$number_of_species)
sd(angiosperms$number_of_species)

# Quantify number of angiosperm species per park, regardless of indigeneity 
angiosperms_all <- angiosperms %>%
  group_by(ParkID) %>%
  summarize(number_of_species=sum(number_of_species)) %>%
  mutate(group="Angiosperms all")

# Create data set from the above data sets: pollinators, angiosperms, and angiosperms_all. For each park, 
# the data set includes the number of pollinator species, the number of angiosperms total, and the number of native
# versus non-native angiosperms. 
species_count <- pollinators_filtered %>%
  bind_rows(angiosperms_all) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Native")) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Non-native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Non-native")) %>%
  mutate(group = if_else(
    group == "Pollinators" & !is.na(group_name),
    paste(group, group_name, sep = "_"),
    group
  )) %>%
  dplyr::select(-number_of_obs, -group_name)
species_count

# repeat for butterflies and moths
species_count_lep <- pollinators_filtered_lep %>%
  bind_rows(angiosperms_all) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Native")) %>%
  bind_rows(angiosperms %>%
              dplyr::filter(status=="Non-native") %>%
              dplyr::select(-status) %>%
              mutate(group="Angiosperms Non-native")) %>%
  mutate(group = if_else(
    !is.na(butterfly_moth),
    paste(butterfly_moth),
    group
  )) %>%
  dplyr::select(-number_of_obs, -butterfly_moth)
species_count_lep

# Count the total number of pollinator and angiosperm observations per park
annotated_data <- data_taxon %>%
  filter(ParkID %in% pollinators_filtered$ParkID) %>%
  left_join(non_natives, by=c("scientific_name"="taxon.name")) %>%
  mutate(status = replace_na(status, "Native"))

# Now summarize by ParkID and group/status
observations_count <- annotated_data %>%
  mutate(group = as.character(group)) %>%
  group_by(ParkID) %>%
  summarise(
    pollinator_obs_lep = sum(group == "Pollinators" & group_name == "Lepidoptera"),
    pollinator_obs_ap = sum(group == "Pollinators" & group_name == "Apoidea"),
    pollinator_obs_but = sum(taxon_family_name %in% c("Hesperiidae", "Papilionidae",
                                                      "Pieridae", "Lycaenidae",
                                                      "Riodinidae", "Nymphalidae")),
    pollinator_obs_moth = sum(taxon_order_name == "Lepidoptera" &
                                !taxon_family_name %in% c("Hesperiidae", "Papilionidae",
                                                          "Pieridae", "Lycaenidae",
                                                          "Riodinidae", "Nymphalidae")),
    angiosperm_obs = sum(group == "Angiosperms"),
    native_angiosperm_obs = sum(group == "Angiosperms" & status == "Native"),
    nonnative_angiosperm_obs = sum(group == "Angiosperms" & status == "Non-native"),
    .groups = "drop"
  ) %>%
  mutate(pollinator_obs_lep=ifelse(pollinator_obs_lep<50, NA, pollinator_obs_lep),
         pollinator_obs_ap=ifelse(pollinator_obs_ap<50, NA, pollinator_obs_ap),
         pollinator_obs_but=ifelse(pollinator_obs_but<50, NA, pollinator_obs_but),
         pollinator_obs_moth=ifelse(pollinator_obs_moth<50, NA, pollinator_obs_moth))
observations_count


# Reformat species richness data (observation_count dataset) to wide format. This data set quantifies 
# the biodiversity (richness), and number of pollinator/angiosperm observations for the final
# set of greenspaces. Also, add the observation counts of pollinators and angiosperms to the data frame
species_richness_wide <- species_count %>%
  pivot_wider(names_from = "group", values_from = "number_of_species", 
              values_fill = list(number_of_species = NA)) %>%
  left_join(observations_count, by = "ParkID") 

# add urban area name
urban_area_names <- parks_gbif_assigned %>% 
  dplyr::select(ParkID, longitude, latitude, urban_area) %>%
  group_by(ParkID) %>%
  summarise(lon=first(longitude), lat=first(latitude), urban_area=first(urban_area))

species_richness_wide <- species_richness_wide %>%
  left_join(urban_area_names, by="ParkID")


# repeate for butterflies and moths
species_richness_wide_lep <- species_count_lep %>%
  pivot_wider(names_from = "group", values_from = "number_of_species", 
              values_fill = list(number_of_species = NA)) %>%
  left_join(observations_count, by = c("ParkID")) 

# add urban area name
urban_area_names <- parks_gbif_assigned %>% 
  dplyr::select(ParkID, longitude, latitude, urban_area) %>%
  group_by(ParkID) %>%
  summarise(lon=first(longitude), lat=first(latitude), urban_area=first(urban_area))

species_richness_wide_lep <- species_richness_wide_lep %>%
  left_join(urban_area_names, by="ParkID")

## Lepidoptera ------------------------

# model species richness for Lepidoptera

### All Angiosperms ------------------------

# log transform skewed values
species_richness_wide_a <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Lepidoptera)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_lep),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area)
  ) %>%
  rename(Angiosperms.all=`Angiosperms all`) 

# check for multicolinearity before modeling
vars <- species_richness_wide_a %>%
  as.data.frame() %>%
  dplyr::select(Angiosperms.all, log_angiosperm_obs, log_pollinator_obs, Pollinators_Lepidoptera)

cor(vars, use = "complete.obs")
# nothing too concerning in terms of the predictor variables

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio <- gam(Angiosperms.all ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
               method="REML",
               data = species_richness_wide_a, 
               family = nb())
summary(m_angio)
gam.check(m_angio)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_a$Angiosperms_all.fitted <- predict(m_angio, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll <- gam(Pollinators_Lepidoptera ~ Angiosperms_all.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll)
gam.check(m_poll)

plot(m_poll)

# let's compare this to the null model 
m_poll_null <- gam(Pollinators_Lepidoptera ~  log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_a, 
                   family = nb())
summary(m_poll_null)

AIC(m_poll_null)-AIC(m_poll)


# now visualize the relationship

# Predict with standard errors
pred <- predict(m_poll, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_a,
             aes(x = Angiosperms_all.fitted, y = Pollinators_Lepidoptera),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_all.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_all.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Angiosperm Species Richness",
       y = "Lepidoptera Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$Pollinators_Lepidoptera), 75))

### Proportion of Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_np <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Lepidoptera)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_lep),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_native = `Angiosperms Native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_np$prop_native)

# check for multicolinearity before modeling
vars <- species_richness_wide_np %>%
  as.data.frame() %>%
  dplyr::select(prop_native, log_angiosperm_obs, log_pollinator_obs, Pollinators_Lepidoptera)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_np <- gam(prop_native ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_np, 
                  family = nb())
summary(m_angio_np)
gam.check(m_angio_np)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_np$prop_native.fitted <- predict(m_angio_np, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_np <- gam(Pollinators_Lepidoptera ~ prop_native.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_np, 
                 family = nb())
summary(m_poll_np)
gam.check(m_poll_np)

# compare to a null model
m_poll_np_null <- gam(Pollinators_Lepidoptera ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                      method="REML",
                      data = species_richness_wide_np, 
                      family = nb())
summary(m_poll_np_null)

AIC(m_poll_np_null)-AIC(m_poll_np)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_np,
               data.frame(
                 prop_native.fitted = seq(min(prop_native.fitted, na.rm=TRUE),
                                          max(prop_native.fitted, na.rm=TRUE),
                                          length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_np, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_np,
             aes(x = prop_native.fitted, y = Pollinators_Lepidoptera),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_native.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_native.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Native Species",
       y = "Lepidoptera Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_np$Pollinators_Lepidoptera), 75))

### Proportion of Non-Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_nnp <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Lepidoptera)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_lep),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_nonnative = `Angiosperms Non-native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_nnp$prop_nonnative)

# check for multicolinearity before modeling
vars <- species_richness_wide_nnp %>%
  as.data.frame() %>%
  dplyr::select(prop_nonnative, log_angiosperm_obs, log_pollinator_obs, Pollinators_Lepidoptera)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_nnp <- gam(prop_nonnative ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_nnp, 
                   family = nb())
summary(m_angio_nnp)
gam.check(m_angio_nnp)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_nnp$prop_nonnative.fitted <- predict(m_angio_nnp, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_nnp <- gam(Pollinators_Lepidoptera ~ prop_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_nnp, 
                  family = nb())
summary(m_poll_nnp)
gam.check(m_poll_nnp)

# compare to a null model
m_poll_nnp_null <- gam(Pollinators_Lepidoptera ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                       method="REML",
                       data = species_richness_wide_nnp, 
                       family = nb())
summary(m_poll_nnp_null)

AIC(m_poll_nnp_null)-AIC(m_poll_nnp)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_nnp,
               data.frame(
                 prop_nonnative.fitted = seq(min(prop_nonnative.fitted, na.rm=TRUE),
                                             max(prop_nonnative.fitted, na.rm=TRUE),
                                             length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_nnp, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nnp,
             aes(x = prop_nonnative.fitted, y = Pollinators_Lepidoptera),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_nonnative.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_nonnative.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Non-Native Species",
       y = "Lepidoptera Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$Pollinators_Lepidoptera), 75))

## Butterflies ------------------------

### All Angiosperms ------------------------

# log transform skewed values
species_richness_wide_a <- species_richness_wide_lep %>%
  filter(complete.cases(butterfly)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_but),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area)
  ) %>%
  rename(Angiosperms.all=`Angiosperms all`) 

# check for multicolinearity before modeling
vars <- species_richness_wide_a %>%
  as.data.frame() %>%
  dplyr::select(Angiosperms.all, log_angiosperm_obs, log_pollinator_obs, butterfly)

cor(vars, use = "complete.obs")
# nothing too concerning in terms of the predictor variables

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio <- gam(Angiosperms.all ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
               method="REML",
               data = species_richness_wide_a, 
               family = nb())
summary(m_angio)
gam.check(m_angio)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_a$Angiosperms_all.fitted <- predict(m_angio, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll <- gam(butterfly ~ Angiosperms_all.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll)
gam.check(m_poll)

# let's compare this to the null model 
m_poll_null <- gam(butterfly ~  log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_a, 
                   family = nb())
summary(m_poll_null)

AIC(m_poll_null)-AIC(m_poll)


# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_a,
               data.frame(
                 Angiosperms_all.fitted = seq(min(Angiosperms_all.fitted, na.rm=TRUE),
                                              max(Angiosperms_all.fitted, na.rm=TRUE),
                                              length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_a,
             aes(x = Angiosperms_all.fitted, y = butterfly),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_all.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_all.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Angiosperm Species Richness",
       y = "Butterfly Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$butterfly), 75))

### Proportion of Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_np <- species_richness_wide_lep %>%
  filter(complete.cases(butterfly)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_but),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_native = `Angiosperms Native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_np$prop_native)

# check for multicolinearity before modeling
vars <- species_richness_wide_np %>%
  as.data.frame() %>%
  dplyr::select(prop_native, log_angiosperm_obs, log_pollinator_obs, butterfly)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_np <- gam(prop_native ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_np, 
                  family = nb())
summary(m_angio_np)
gam.check(m_angio_np)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_np$prop_native.fitted <- predict(m_angio_np, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_np <- gam(butterfly ~ prop_native.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_np, 
                 family = nb())
summary(m_poll_np)
gam.check(m_poll_np)

# compare to a null model
m_poll_np_null <- gam(butterfly ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                      method="REML",
                      data = species_richness_wide_np, 
                      family = nb())
summary(m_poll_np_null)

AIC(m_poll_np_null)-AIC(m_poll_np)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_np,
               data.frame(
                 prop_native.fitted = seq(min(prop_native.fitted, na.rm=TRUE),
                                          max(prop_native.fitted, na.rm=TRUE),
                                          length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_np, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_np,
             aes(x = prop_native.fitted, y = butterfly),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_native.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_native.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Native Species",
       y = "Butterfly Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_np$butterfly), 75))

### Proportion of Non-Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_nnp <- species_richness_wide_lep %>%
  filter(complete.cases(butterfly)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_but),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_nonnative = `Angiosperms Non-native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_nnp$prop_nonnative)

# check for multicolinearity before modeling
vars <- species_richness_wide_nnp %>%
  as.data.frame() %>%
  dplyr::select(prop_nonnative, log_angiosperm_obs, log_pollinator_obs, butterfly)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_nnp <- gam(prop_nonnative ~ log_angiosperm_obs + s(lat, lon, k = 50, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_nnp, 
                   family = nb())
summary(m_angio_nnp)
gam.check(m_angio_nnp)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_nnp$prop_nonnative.fitted <- predict(m_angio_nnp, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_nnp <- gam(butterfly ~ prop_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_nnp, 
                  family = nb())
summary(m_poll_nnp)
gam.check(m_poll_nnp)

# compare to a null model
m_poll_nnp_null <- gam(butterfly ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                       method="REML",
                       data = species_richness_wide_nnp, 
                       family = nb())
summary(m_poll_nnp_null)

AIC(m_poll_nnp_null)-AIC(m_poll_nnp)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_nnp,
               data.frame(
                 prop_nonnative.fitted = seq(min(prop_nonnative.fitted, na.rm=TRUE),
                                             max(prop_nonnative.fitted, na.rm=TRUE),
                                             length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_nnp, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nnp,
             aes(x = prop_nonnative.fitted, y = butterfly),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_nonnative.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_nonnative.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Non-Native Species",
       y = "Butterfly Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$butterfly), 75))

## Moths ------------------------

### All Angiosperms ------------------------

# log transform skewed values
species_richness_wide_a <- species_richness_wide_lep %>%
  filter(complete.cases(moth)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_moth),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area)
  ) %>%
  rename(Angiosperms.all=`Angiosperms all`) 

hist(species_richness_wide_a$log_pollinator_obs)

# check for multicolinearity before modeling
vars <- species_richness_wide_a %>%
  as.data.frame() %>%
  dplyr::select(Angiosperms.all, log_angiosperm_obs, log_pollinator_obs, moth)

cor(vars, use = "complete.obs")
# nothing too concerning in terms of the predictor variables

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio <- gam(Angiosperms.all ~ log_angiosperm_obs + s(lat, lon, k = 27, bs="tp"), 
               method="REML",
               data = species_richness_wide_a, 
               family = nb())
summary(m_angio)
gam.check(m_angio)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_a$Angiosperms_all.fitted <- predict(m_angio, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll <- gam(moth ~ Angiosperms_all.fitted + log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll)
gam.check(m_poll)


# let's compare this to the null model 
m_poll_null <- gam(moth ~  log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_a, 
                   family = nb())
summary(m_poll_null)

AIC(m_poll_null)-AIC(m_poll)


# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_a,
               data.frame(
                 Angiosperms_all.fitted = seq(min(Angiosperms_all.fitted, na.rm=TRUE),
                                              max(Angiosperms_all.fitted, na.rm=TRUE),
                                              length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_a,
             aes(x = Angiosperms_all.fitted, y = moth),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_all.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_all.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Angiosperm Species Richness",
       y = "Moth Species Richness") +
  theme_classic()  +
  coord_cartesian(ylim = c(min(species_richness_wide_a$moth), 75))

### Proportion of Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_np <- species_richness_wide_lep %>%
  filter(complete.cases(moth)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_moth),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_native = `Angiosperms Native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_np$prop_native)

# check for multicolinearity before modeling
vars <- species_richness_wide_np %>%
  as.data.frame() %>%
  dplyr::select(prop_native, log_angiosperm_obs, log_pollinator_obs, moth)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_np <- gam(prop_native ~ log_angiosperm_obs + s(lat, lon, k = 27, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_np, 
                  family = nb())
summary(m_angio_np)
gam.check(m_angio_np)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_np$prop_native.fitted <- predict(m_angio_np, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_np <- gam(moth ~ prop_native.fitted + log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_np, 
                 family = nb())
summary(m_poll_np)
gam.check(m_poll_np)

# compare to a null model
m_poll_np_null <- gam(moth ~ log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
                      method="REML",
                      data = species_richness_wide_np, 
                      family = nb())
summary(m_poll_np_null)

AIC(m_poll_np_null)-AIC(m_poll_np)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_np,
               data.frame(
                 prop_native.fitted = seq(min(prop_native.fitted, na.rm=TRUE),
                                          max(prop_native.fitted, na.rm=TRUE),
                                          length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_np, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_np,
             aes(x = prop_native.fitted, y = moth),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_native.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_native.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Native Species",
       y = "Moth Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_np$moth), 75))

### Proportion of Non-Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_nnp <- species_richness_wide_lep %>%
  filter(complete.cases(moth)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_moth),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_nonnative = `Angiosperms Non-native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_nnp$prop_nonnative)

# check for multicolinearity before modeling
vars <- species_richness_wide_nnp %>%
  as.data.frame() %>%
  dplyr::select(prop_nonnative, log_angiosperm_obs, log_pollinator_obs, moth)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_nnp <- gam(prop_nonnative ~ log_angiosperm_obs + s(lat, lon, k = 27, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_nnp, 
                   family = nb())
summary(m_angio_nnp)
gam.check(m_angio_nnp)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_nnp$prop_nonnative.fitted <- predict(m_angio_nnp, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_nnp <- gam(moth ~ prop_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_nnp, 
                  family = nb())
summary(m_poll_nnp)
gam.check(m_poll_nnp)

# compare to a null model
m_poll_nnp_null <- gam(moth ~ log_pollinator_obs + s(lat, lon, k=27, bs="tp"), 
                       method="REML",
                       data = species_richness_wide_nnp, 
                       family = nb())
summary(m_poll_nnp_null)

AIC(m_poll_nnp_null)-AIC(m_poll_nnp)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_nnp,
               data.frame(
                 prop_nonnative.fitted = seq(min(prop_nonnative.fitted, na.rm=TRUE),
                                             max(prop_nonnative.fitted, na.rm=TRUE),
                                             length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_nnp, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nnp,
             aes(x = prop_nonnative.fitted, y = moth),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_nonnative.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_nonnative.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Non-Native Species",
       y = "Moth Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$moth), 100))

## Apoidea ------------------------

### All Angiosperms ------------------------

# log transform skewed values
species_richness_wide_a <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Apoidea)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_ap),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area)
  ) %>%
  rename(Angiosperms.all=`Angiosperms all`) 

# check for multicolinearity before modeling
vars <- species_richness_wide_a %>%
  as.data.frame() %>%
  dplyr::select(Angiosperms.all, log_angiosperm_obs, log_pollinator_obs, Pollinators_Apoidea)

cor(vars, use = "complete.obs")
# nothing too concerning in terms of the predictor variables

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio <- gam(Angiosperms.all ~ log_angiosperm_obs + s(lat, lon, k = 12, bs="tp"), 
               method="REML",
               data = species_richness_wide_a, 
               family = nb())
summary(m_angio)
gam.check(m_angio)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_a$Angiosperms_all.fitted <- predict(m_angio, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll <- gam(Pollinators_Apoidea ~ Angiosperms_all.fitted + log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll)
gam.check(m_poll)

# let's compare this to the null model 
m_poll_null <- gam(Pollinators_Apoidea ~  log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_a, 
                   family = nb())
summary(m_poll_null)

AIC(m_poll_null)-AIC(m_poll)


# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_a,
               data.frame(
                 Angiosperms_all.fitted = seq(min(Angiosperms_all.fitted, na.rm=TRUE),
                                              max(Angiosperms_all.fitted, na.rm=TRUE),
                                              length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_a,
             aes(x = Angiosperms_all.fitted, y = Pollinators_Apoidea),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_all.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_all.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Angiosperm Species Richness",
       y = "Apoidea Species Richness") +
  theme_classic()

### Proportion of Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_np <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Apoidea)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_ap),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_native = `Angiosperms Native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_np$prop_native)

# check for multicolinearity before modeling
vars <- species_richness_wide_np %>%
  as.data.frame() %>%
  dplyr::select(prop_native, log_angiosperm_obs, log_pollinator_obs, Pollinators_Apoidea)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_np <- gam(prop_native ~ log_angiosperm_obs + s(lat, lon, k = 12, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_np, 
                  family = nb())
summary(m_angio_np)
gam.check(m_angio_np)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_np$prop_native.fitted <- predict(m_angio_np, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_np <- gam(Pollinators_Apoidea ~ prop_native.fitted + log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_np, 
                 family = nb())
summary(m_poll_np)
gam.check(m_poll_np)

# compare to a null model
m_poll_np_null <- gam(Pollinators_Lepidoptera ~ log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
                      method="REML",
                      data = species_richness_wide_np, 
                      family = nb())
summary(m_poll_np_null)

AIC(m_poll_np_null)-AIC(m_poll_np)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_np,
               data.frame(
                 prop_native.fitted = seq(min(prop_native.fitted, na.rm=TRUE),
                                          max(prop_native.fitted, na.rm=TRUE),
                                          length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_np, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_np,
             aes(x = prop_native.fitted, y = Pollinators_Apoidea),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_native.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_native.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Native Species",
       y = "Apoidea Species Richness") +
  theme_classic() 

### Proportion of Non-Native Angiosperms ------------------------

# log transform skewed values
species_richness_wide_nnp <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Apoidea)) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs_ap),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_nonnative = `Angiosperms Non-native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_nnp$prop_nonnative)

# check for multicolinearity before modeling
vars <- species_richness_wide_nnp %>%
  as.data.frame() %>%
  dplyr::select(prop_nonnative, log_angiosperm_obs, log_pollinator_obs, Pollinators_Apoidea)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_nnp <- gam(prop_nonnative ~ log_angiosperm_obs + s(lat, lon, k = 12, bs="tp"), 
                   method="REML",
                   data = species_richness_wide_nnp, 
                   family = nb())
summary(m_angio_nnp)
gam.check(m_angio_nnp)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_nnp$prop_nonnative.fitted <- predict(m_angio_nnp, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_nnp <- gam(Pollinators_Apoidea ~ prop_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
                  method="REML",
                  data = species_richness_wide_nnp, 
                  family = nb())
summary(m_poll_nnp)
gam.check(m_poll_nnp)

# compare to a null model
m_poll_nnp_null <- gam(Pollinators_Lepidoptera ~ log_pollinator_obs + s(lat, lon, k=12, bs="tp"), 
                       method="REML",
                       data = species_richness_wide_nnp, 
                       family = nb())
summary(m_poll_nnp_null)

AIC(m_poll_nnp_null)-AIC(m_poll_nnp)

# now visualize the relationship

# Create a new dataset varying Angiosperms_all.fitted
newdat <- with(species_richness_wide_nnp,
               data.frame(
                 prop_nonnative.fitted = seq(min(prop_nonnative.fitted, na.rm=TRUE),
                                             max(prop_nonnative.fitted, na.rm=TRUE),
                                             length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_nnp, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nnp,
             aes(x = prop_nonnative.fitted, y = Pollinators_Apoidea),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = prop_nonnative.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = prop_nonnative.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Proportion of Non-Native Species",
       y = "Apoidea Species Richness") +
  theme_classic()

# Model pollinator richness with landcover attributes ------------------------

## Lepidoptera ------------------------

### Adding landcover values ----------------------------------------------------

# Read in data from dynamic world
landcover <- read_csv("Data/Dynamic_World_Habitat/dynamic_world_habitat_fl.csv")

# add this data to our data frame
lc_data <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Lepidoptera)) %>%
  left_join(landcover, by = "ParkID") %>%
  group_by(ParkID) %>%
  summarise(Shape_Area=sum(Shape_Area), trees=mean(trees), snow_and_ice=mean(snow_and_ice),
            water=mean(water), shrub_and_scrub=mean(shrub_and_scrub), built=mean(built),
            crops=mean(crops), grass=mean(grass), flooded_vegetation=mean(flooded_vegetation),
            bare=mean(bare), pollinator_obs=sum(pollinator_obs_lep), Pollinators=sum(Pollinators_Lepidoptera),
            lon=first(lon), lat=first(lat))

# let's examine the landcover data to figure out which should be included in the model
summary(lc_data[,c("trees", "snow_and_ice", "water", "shrub_and_scrub", "built", 
                   "bare", "crops", "grass", "flooded_vegetation")])
# the only variables with a median > 0 are trees, water, built, and grass so we will use these variables

# the features that are most common in parks are built, grass, trees, and water so we will model with those
# let's examine their distribution
hist(lc_data$built)
hist(lc_data$grass)
hist(lc_data$trees)
hist(lc_data$water)

### Model the data ------------------------

# log transform number of observations for both pollinators and angiosperms as well as park size
richness_covariates <- lc_data %>%
  mutate(log_obs_pol=log(pollinator_obs),
         log_parksize=log(Shape_Area/10000)) %>%
  dplyr::select(Pollinators, built, grass, trees, water, log_parksize, log_obs_pol, lat, lon)


# let's start with a loaded model
model_all <- gam(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_all)
gam.check(model_all)
# gam check looks good

# now let's try different combinations of the model
model_null <- gam(Pollinators ~ log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_null)
gam.check(model_null)

model_1 <- gam(Pollinators ~ built + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_1)
gam.check(model_1)

model_2 <- gam(Pollinators ~ grass + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_2)
gam.check(model_2)

model_3 <- gam(Pollinators ~ trees + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_3)
gam.check(model_3)

model_4 <- gam(Pollinators ~ water + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_4)
gam.check(model_4)

model_5 <- gam(Pollinators ~ log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_5)
gam.check(model_5)

model_6 <- gam(Pollinators ~ built + trees + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_6)
gam.check(model_6)

AIC(model_all, model_null, model_1, model_2, model_3, model_4, model_5, model_6)
# The model with impervious cover explains the most

# top model compared to null
AIC(model_null)-AIC(model_all)

# create a function for predicting the data for plotting purposes
get_prediction_data <- function(formula, data, predictor_name, label) {
  library(mgcv)
  library(dplyr)
  
  model <- gam(formula, method="REML",
               data = data,
               family = nb())
  
  # Create new data for prediction
  newdata <- data.frame(
    built = mean(data$built),
    grass = mean(data$grass),
    trees = mean(data$trees),
    water = mean(data$water),
    log_parksize = mean(data$log_parksize),
    log_obs_pol = mean(data$log_obs_pol),
    lat = mean(data$lat),
    lon = mean(data$lon)
  )
  
  # vary the predictor of interest
  newdata <- newdata[rep(1, 100), ] # repeat 100 rows
  newdata[[predictor_name]] <- seq(min(data[[predictor_name]]), max(data[[predictor_name]]), length.out = 100)
  
  # Predict on link scale with SE
  pred <- predict(model, newdata, type = "link", se.fit = TRUE)
  
  newdata <- newdata %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      ci_lower = exp(fit - 1.96 * se),
      ci_upper = exp(fit + 1.96 * se),
      Pollinators = exp(fit),
      model = label,
      x = .data[[predictor_name]]
    )
  
  return(newdata %>% dplyr::select(x, Pollinators, ci_lower, ci_upper, model))
}

# Plot the data
pred_built <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "built", "Impervious Cover")
pred_grass <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "grass", "Grass Cover")
pred_trees <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "trees", "Tree Cover")
pred_water <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "water", "Water Cover")
pred_parksize <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "log_parksize", "Log10(Greenspace Size)")

combined_preds <- bind_rows(pred_built, pred_grass, pred_trees, pred_water, pred_parksize)


# adjust the pollinator richness based on number of observations

# Get predicted pollinator richness based on sampling effort only
richness_covariates$pred_null <- predict(model_null, type = "response")

# Reshape richness_covariates to long format for predictor variables
observed_long <- richness_covariates %>%
  dplyr::select(built, grass, trees, water, log_parksize, pred_null) %>%
  pivot_longer(
    cols = c(built, grass, trees, water, log_parksize),
    names_to = "model",
    values_to = "x"
  ) %>%
  mutate(model = case_when(
    model == "built" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.002 SE, <i>P</i> = 0.003</span>",
    model == "grass" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.003 SE, <i>P</i> = 0.561</span>",
    model == "trees" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.002 SE, <i>P</i> = 0.002</span>",
    model == "water" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0009 ± 0.003 SE, <i>P</i> = 0.729</span>",
    model == "log_parksize" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.035 ± 0.032 SE, <i>P</i> = 0.023</span>"
  ))

combined_preds <- combined_preds %>%
  mutate(model = case_when(
    model == "Impervious Cover" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.002 SE, <i>P</i> = 0.003</span>",
    model == "Grass Cover" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.003 SE, <i>P</i> = 0.561</span>",
    model == "Tree Cover" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.002 SE, <i>P</i> = 0.002</span>",
    model == "Water Cover" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0009 ± 0.003 SE, <i>P</i> = 0.729</span>",
    model == "Log10(Greenspace Size)" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.035 ± 0.032 SE, <i>P</i> = 0.023</span>"
  ))

# reorder plots so greenspace area is last
facet_levels <- c(
  "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.002 SE, <i>P</i> = 0.003</span>",
  "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.003 SE, <i>P</i> = 0.561</span>",
  "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.002 SE, <i>P</i> = 0.002</span>",
  "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0009 ± 0.003 SE, <i>P</i> = 0.729</span>",
  "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.035 ± 0.032 SE, <i>P</i> = 0.023</span>"
)

observed_long <- observed_long %>%
  mutate(model = factor(model, levels = facet_levels))

combined_preds <- combined_preds %>%
  mutate(model = factor(model, levels = facet_levels))

# Now plot using this observed_long for points (x = predictor values, y = adjusted richness residuals)
ggplot() +
  geom_ribbon(data = combined_preds, aes(x = x, ymin = ci_lower, ymax = ci_upper), fill = "lightblue") +
  geom_line(data = combined_preds, aes(x = x, y = Pollinators), color = "blue", size = 1) +
  geom_point(data = observed_long, aes(x = x, y = pred_null), color = "black", alpha = 0.6) +
  facet_wrap(~ model, scales = "free_x", ncol = 2) +
  labs(x = "Predictor Variable",
       y = "Lepidoptera Species Richness") +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("Figures/Figure_S7.jpeg", height=8, width=8, units="in")

# Create a parameter estimate plot from the full model
# Tidy the model for plotting
model_df <- broom::tidy(model_all, conf.int = TRUE, parametric=TRUE)

# Filter out the intercept if you don't want to include it
model_df <- model_df %>%
  filter(term != "(Intercept)",
         term != "log_obs_pol") %>%
  mutate(term = case_match(term,
                           "built" ~ "Impervious Cover (%)*",
                           "grass" ~ "Grass Cover (%)",
                           "trees" ~ "Tree Cover (%)**",
                           "water" ~ "Water Cover (%)",
                           "log_parksize" ~ "Log10(Greenspace Size)*"))

# Create the parameter estimate plot
ggplot(model_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Coefficient Estimate (± 95% CI)",
    y = ""
  ) +
  theme_classic(base_size = 14)

## Butterflies ------------------------

### Adding landcover values ----------------------------------------------------

# Read in data from dynamic world
landcover <- read_csv("Data/Dynamic_World_Habitat/dynamic_world_habitat_fl.csv")

# add this data to our data frame
lc_data <- species_richness_wide_lep %>%
  filter(complete.cases(butterfly)) %>%
  left_join(landcover, by = "ParkID") %>%
  group_by(ParkID) %>%
  summarise(Shape_Area=sum(Shape_Area), trees=mean(trees), snow_and_ice=mean(snow_and_ice),
            water=mean(water), shrub_and_scrub=mean(shrub_and_scrub), built=mean(built),
            crops=mean(crops), grass=mean(grass), flooded_vegetation=mean(flooded_vegetation),
            bare=mean(bare), pollinator_obs=sum(pollinator_obs_but), Pollinators=sum(butterfly),
            lon=first(lon), lat=first(lat))

# let's examine the landcover data to figure out which should be included in the model
summary(lc_data[,c("trees", "snow_and_ice", "water", "shrub_and_scrub", "built", 
                   "bare", "crops", "grass", "flooded_vegetation")])
# the only variables with a median > 0 are trees, water, built, and grass so we will use these variables

# the features that are most common in parks are built, grass, trees, and water so we will model with those
# let's examine their distribution
hist(lc_data$built)
hist(lc_data$grass)
hist(lc_data$trees)
hist(lc_data$water)

### Model the data ------------------------

# log transform number of observations for both pollinators and angiosperms as well as park size
richness_covariates <- lc_data %>%
  mutate(log_obs_pol=log(pollinator_obs),
         log_parksize=log(Shape_Area/10000)) %>%
  dplyr::select(Pollinators, built, grass, trees, water, log_parksize, log_obs_pol, lat, lon)


# let's start with a loaded model
model_all <- gam(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_all)
gam.check(model_all)
# gam check looks good

# now let's try different combinations of the model
model_null <- gam(Pollinators ~ log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
                  method="REML",
                  data = richness_covariates,
                  family=nb())
summary(model_null)
gam.check(model_null)

model_1 <- gam(Pollinators ~ built + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_1)
gam.check(model_1)

model_2 <- gam(Pollinators ~ grass + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_2)
gam.check(model_2)

model_3 <- gam(Pollinators ~ trees + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_3)
gam.check(model_3)

model_4 <- gam(Pollinators ~ water + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_4)
gam.check(model_4)

model_5 <- gam(Pollinators ~ log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_5)
gam.check(model_5)

model_6 <- gam(Pollinators ~ built + trees + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_6)
gam.check(model_6)

AIC(model_all, model_null, model_1, model_2, model_3, model_4, model_5, model_6)
# The model with impervious cover explains the most

# top model compared to null
AIC(model_null)-AIC(model_all)

# create a function for predicting the data for plotting purposes
get_prediction_data <- function(formula, data, predictor_name, label) {
  library(mgcv)
  library(dplyr)
  
  model <- gam(formula, method="REML",
               data = data,
               family = nb())
  
  # Create new data for prediction
  newdata <- data.frame(
    built = mean(data$built),
    grass = mean(data$grass),
    trees = mean(data$trees),
    water = mean(data$water),
    log_parksize = mean(data$log_parksize),
    log_obs_pol = mean(data$log_obs_pol),
    lat = mean(data$lat),
    lon = mean(data$lon)
  )
  
  # vary the predictor of interest
  newdata <- newdata[rep(1, 100), ] # repeat 100 rows
  newdata[[predictor_name]] <- seq(min(data[[predictor_name]]), max(data[[predictor_name]]), length.out = 100)
  
  # Predict on link scale with SE
  pred <- predict(model, newdata, type = "link", se.fit = TRUE)
  
  newdata <- newdata %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      ci_lower = exp(fit - 1.96 * se),
      ci_upper = exp(fit + 1.96 * se),
      Pollinators = exp(fit),
      model = label,
      x = .data[[predictor_name]]
    )
  
  return(newdata %>% dplyr::select(x, Pollinators, ci_lower, ci_upper, model))
}

# Plot the data
pred_built <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "built", "Impervious Cover")
pred_grass <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "grass", "Grass Cover")
pred_trees <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "trees", "Tree Cover")
pred_water <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "water", "Water Cover")
pred_parksize <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=50, bs="tp"), richness_covariates, "log_parksize", "Log10(Greenspace Size)")

combined_preds <- bind_rows(pred_built, pred_grass, pred_trees, pred_water, pred_parksize)


# adjust the pollinator richness based on number of observations

# Get predicted pollinator richness based on sampling effort only
richness_covariates$pred_null <- predict(model_null, type = "response")

# Reshape richness_covariates to long format for predictor variables
observed_long <- richness_covariates %>%
  dplyr::select(built, grass, trees, water, log_parksize, pred_null) %>%
  pivot_longer(
    cols = c(built, grass, trees, water, log_parksize),
    names_to = "model",
    values_to = "x"
  ) %>%
  mutate(model = case_when(
    model == "built" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.396</span>",
    model == "grass" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0001 ± 0.002 SE, <i>P</i> = 0.940</span>",
    model == "trees" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0007 ± 0.001 SE, <i>P</i> = 0.575</span>",
    model == "water" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.635</span>",
    model == "log_parksize" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.012 ± 0.015 SE, <i>P</i> = 0.393</span>"
  ))

combined_preds <- combined_preds %>%
  mutate(model = case_when(
    model == "Impervious Cover" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.396</span>",
    model == "Grass Cover" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0001 ± 0.002 SE, <i>P</i> = 0.940</span>",
    model == "Tree Cover" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0007 ± 0.001 SE, <i>P</i> = 0.575</span>",
    model == "Water Cover" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.635</span>",
    model == "Log10(Greenspace Size)" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.012 ± 0.015 SE, <i>P</i> = 0.393</span>"
  ))

# reorder plots so greenspace area is last
facet_levels <- c(
  "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.396</span>",
  "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0001 ± 0.002 SE, <i>P</i> = 0.940</span>",
  "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.0007 ± 0.001 SE, <i>P</i> = 0.575</span>",
  "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.002 SE, <i>P</i> = 0.635</span>",
  "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.012 ± 0.015 SE, <i>P</i> = 0.393</span>"
)

observed_long <- observed_long %>%
  mutate(model = factor(model, levels = facet_levels))

combined_preds <- combined_preds %>%
  mutate(model = factor(model, levels = facet_levels))

# Now plot using this observed_long for points (x = predictor values, y = adjusted richness residuals)
ggplot() +
  geom_ribbon(data = combined_preds, aes(x = x, ymin = ci_lower, ymax = ci_upper), fill = "lightblue") +
  geom_line(data = combined_preds, aes(x = x, y = Pollinators), color = "blue", size = 1) +
  geom_point(data = observed_long, aes(x = x, y = pred_null), color = "black", alpha = 0.6) +
  facet_wrap(~ model, scales = "free_x", ncol = 2) +
  labs(x = "Predictor Variable",
       y = "Butterfly Species Richness") +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("Figures/Figure_S8.jpeg", height=8, width=8, units="in")

# Create a parameter estimate plot from the full model
# Tidy the model for plotting
model_df <- broom::tidy(model_all, conf.int = TRUE, parametric=TRUE)

# Filter out the intercept if you don't want to include it
model_df <- model_df %>%
  filter(term != "(Intercept)",
         term != "log_obs_pol") %>%
  mutate(term = case_match(term,
                           "built" ~ "Impervious Cover (%)*",
                           "grass" ~ "Grass Cover (%)",
                           "trees" ~ "Tree Cover (%)**",
                           "water" ~ "Water Cover (%)",
                           "log_parksize" ~ "Log10(Greenspace Size)*"))

# Create the parameter estimate plot
ggplot(model_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Coefficient Estimate (± 95% CI)",
    y = ""
  ) +
  theme_classic(base_size = 14)

## Moths ------------------------

### Adding landcover values ----------------------------------------------------

# Read in data from dynamic world
landcover <- read_csv("Data/Dynamic_World_Habitat/dynamic_world_habitat_fl.csv")

# add this data to our data frame
lc_data <- species_richness_wide_lep %>%
  filter(complete.cases(moth)) %>%
  left_join(landcover, by = "ParkID") %>%
  group_by(ParkID) %>%
  summarise(Shape_Area=sum(Shape_Area), trees=mean(trees), snow_and_ice=mean(snow_and_ice),
            water=mean(water), shrub_and_scrub=mean(shrub_and_scrub), built=mean(built),
            crops=mean(crops), grass=mean(grass), flooded_vegetation=mean(flooded_vegetation),
            bare=mean(bare), pollinator_obs=sum(pollinator_obs_moth), Pollinators=sum(moth),
            lon=first(lon), lat=first(lat))

# let's examine the landcover data to figure out which should be included in the model
summary(lc_data[,c("trees", "snow_and_ice", "water", "shrub_and_scrub", "built", 
                   "bare", "crops", "grass", "flooded_vegetation")])
# the only variables with a median > 0 are trees, water, built, and grass so we will use these variables

# the features that are most common in parks are built, grass, trees, and water so we will model with those
# let's examine their distribution
hist(lc_data$built)
hist(lc_data$grass)
hist(lc_data$trees)
hist(lc_data$water)

### Model the data ------------------------

# log transform number of observations for both pollinators and angiosperms as well as park size
richness_covariates <- lc_data %>%
  mutate(log_obs_pol=log(pollinator_obs),
         log_parksize=log(Shape_Area/10000)) %>%
  dplyr::select(Pollinators, built, grass, trees, water, log_parksize, log_obs_pol, lat, lon)


# let's start with a loaded model
model_all <- gam(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_all)
gam.check(model_all)
# gam check looks good

# now let's try different combinations of the model
model_null <- gam(Pollinators ~ log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
                  method="REML",
                  data = richness_covariates,
                  family=nb())
summary(model_null)
gam.check(model_null)

model_1 <- gam(Pollinators ~ built + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_1)
gam.check(model_1)

model_2 <- gam(Pollinators ~ grass + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_2)
gam.check(model_2)

model_3 <- gam(Pollinators ~ trees + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_3)
gam.check(model_3)

model_4 <- gam(Pollinators ~ water + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_4)
gam.check(model_4)

model_5 <- gam(Pollinators ~ log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_5)
gam.check(model_5)

model_6 <- gam(Pollinators ~ built + trees + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_6)
gam.check(model_6)

AIC(model_all, model_null, model_1, model_2, model_3, model_4, model_5, model_6)
# The model with impervious cover explains the most

# top model compared to null
AIC(model_null)-AIC(model_all)

# create a function for predicting the data for plotting purposes
get_prediction_data <- function(formula, data, predictor_name, label) {
  library(mgcv)
  library(dplyr)
  
  model <- gam(formula, method="REML",
               data = data,
               family = nb())
  
  # Create new data for prediction
  newdata <- data.frame(
    built = mean(data$built),
    grass = mean(data$grass),
    trees = mean(data$trees),
    water = mean(data$water),
    log_parksize = mean(data$log_parksize),
    log_obs_pol = mean(data$log_obs_pol),
    lat = mean(data$lat),
    lon = mean(data$lon)
  )
  
  # vary the predictor of interest
  newdata <- newdata[rep(1, 100), ] # repeat 100 rows
  newdata[[predictor_name]] <- seq(min(data[[predictor_name]]), max(data[[predictor_name]]), length.out = 100)
  
  # Predict on link scale with SE
  pred <- predict(model, newdata, type = "link", se.fit = TRUE)
  
  newdata <- newdata %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      ci_lower = exp(fit - 1.96 * se),
      ci_upper = exp(fit + 1.96 * se),
      Pollinators = exp(fit),
      model = label,
      x = .data[[predictor_name]]
    )
  
  return(newdata %>% dplyr::select(x, Pollinators, ci_lower, ci_upper, model))
}

# Plot the data
pred_built <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), richness_covariates, "built", "Impervious Cover")
pred_grass <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), richness_covariates, "grass", "Grass Cover")
pred_trees <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), richness_covariates, "trees", "Tree Cover")
pred_water <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), richness_covariates, "water", "Water Cover")
pred_parksize <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=27, bs="tp"), richness_covariates, "log_parksize", "Log10(Greenspace Size)")

combined_preds <- bind_rows(pred_built, pred_grass, pred_trees, pred_water, pred_parksize)


# adjust the pollinator richness based on number of observations

# Get predicted pollinator richness based on sampling effort only
richness_covariates$pred_null <- predict(model_null, type = "response")

# Reshape richness_covariates to long format for predictor variables
observed_long <- richness_covariates %>%
  dplyr::select(built, grass, trees, water, log_parksize, pred_null) %>%
  pivot_longer(
    cols = c(built, grass, trees, water, log_parksize),
    names_to = "model",
    values_to = "x"
  ) %>%
  mutate(model = case_when(
    model == "built" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.589 SE, <i>P</i> = 0.691</span>",
    model == "grass" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.011 SE, <i>P</i> = 0.574</span>",
    model == "trees" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.005 SE, <i>P</i> = 0.777</span>",
    model == "water" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.007 SE, <i>P</i> = 0.575</span>",
    model == "log_parksize" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.058 ± 0.031 SE, <i>P</i> = 0.0592</span>"
  ))

combined_preds <- combined_preds %>%
  mutate(model = case_when(
    model == "Impervious Cover" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.589 SE, <i>P</i> = 0.691</span>",
    model == "Grass Cover" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.011 SE, <i>P</i> = 0.574</span>",
    model == "Tree Cover" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.005 SE, <i>P</i> = 0.777</span>",
    model == "Water Cover" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.007 SE, <i>P</i> = 0.575</span>",
    model == "Log10(Greenspace Size)" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.058 ± 0.031 SE, <i>P</i> = 0.0592</span>"
  ))

# reorder plots so greenspace area is last
facet_levels <- c(
  "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.589 SE, <i>P</i> = 0.691</span>",
  "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.011 SE, <i>P</i> = 0.574</span>",
  "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.001 ± 0.005 SE, <i>P</i> = 0.777</span>",
  "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.007 SE, <i>P</i> = 0.575</span>",
  "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.058 ± 0.031 SE, <i>P</i> = 0.0592</span>"
)

observed_long <- observed_long %>%
  mutate(model = factor(model, levels = facet_levels))

combined_preds <- combined_preds %>%
  mutate(model = factor(model, levels = facet_levels))

# Now plot using this observed_long for points (x = predictor values, y = adjusted richness residuals)
ggplot() +
  geom_ribbon(data = combined_preds, aes(x = x, ymin = ci_lower, ymax = ci_upper), fill = "lightblue") +
  geom_line(data = combined_preds, aes(x = x, y = Pollinators), color = "blue", size = 1) +
  geom_point(data = observed_long, aes(x = x, y = pred_null), color = "black", alpha = 0.6) +
  facet_wrap(~ model, scales = "free_x", ncol = 2) +
  labs(x = "Predictor Variable",
       y = "Moth Species Richness") +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("Figures/Figure_S9.jpeg", height=8, width=8, units="in")

# Create a parameter estimate plot from the full model
# Tidy the model for plotting
model_df <- broom::tidy(model_all, conf.int = TRUE, parametric=TRUE)

# Filter out the intercept if you don't want to include it
model_df <- model_df %>%
  filter(term != "(Intercept)",
         term != "log_obs_pol") %>%
  mutate(term = case_match(term,
                           "built" ~ "Impervious Cover (%)*",
                           "grass" ~ "Grass Cover (%)",
                           "trees" ~ "Tree Cover (%)**",
                           "water" ~ "Water Cover (%)",
                           "log_parksize" ~ "Log10(Greenspace Size)*"))

# Create the parameter estimate plot
ggplot(model_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Coefficient Estimate (± 95% CI)",
    y = ""
  ) +
  theme_classic(base_size = 14)

## Apoidea ------------------------

### Adding landcover values ----------------------------------------------------

# Read in data from dynamic world
landcover <- read_csv("Data/Dynamic_World_Habitat/dynamic_world_habitat_fl.csv")

# add this data to our data frame
lc_data <- species_richness_wide %>%
  filter(complete.cases(Pollinators_Apoidea)) %>%
  left_join(landcover, by = "ParkID") %>%
  group_by(ParkID) %>%
  summarise(Shape_Area=sum(Shape_Area), trees=mean(trees), snow_and_ice=mean(snow_and_ice),
            water=mean(water), shrub_and_scrub=mean(shrub_and_scrub), built=mean(built),
            crops=mean(crops), grass=mean(grass), flooded_vegetation=mean(flooded_vegetation),
            bare=mean(bare), pollinator_obs=sum(pollinator_obs_ap), Pollinators=sum(Pollinators_Apoidea),
            lon=first(lon), lat=first(lat))

# let's examine the landcover data to figure out which should be included in the model
summary(lc_data[,c("trees", "snow_and_ice", "water", "shrub_and_scrub", "built", 
                   "bare", "crops", "grass", "flooded_vegetation")])
# the only variables with a median > 0 are trees, water, built, and grass so we will use these variables

# the features that are most common in parks are built, grass, trees, and water so we will model with those
# let's examine their distribution
hist(lc_data$built)
hist(lc_data$grass)
hist(lc_data$trees)
hist(lc_data$water)

### Model the data ------------------------

# log transform number of observations for both pollinators and angiosperms as well as park size
richness_covariates <- lc_data %>%
  mutate(log_obs_pol=log(pollinator_obs),
         log_parksize=log(Shape_Area/10000)) %>%
  dplyr::select(Pollinators, built, grass, trees, water, log_parksize, log_obs_pol, lat, lon)

# let's start with a loaded model
model_all <- gam(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
                 method="REML",
                 data = richness_covariates,
                 family=nb())
summary(model_all)
gam.check(model_all)
# gam check looks good

# now let's try different combinations of the model
model_null <- gam(Pollinators ~ log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
                  method="REML",
                  data = richness_covariates,
                  family=nb())
summary(model_null)
gam.check(model_null)

model_1 <- gam(Pollinators ~ built + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_1)
gam.check(model_1)

model_2 <- gam(Pollinators ~ grass + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_2)
gam.check(model_2)

model_3 <- gam(Pollinators ~ trees + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_3)
gam.check(model_3)

model_4 <- gam(Pollinators ~ water + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_4)
gam.check(model_4)

model_5 <- gam(Pollinators ~ log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_5)
gam.check(model_5)

model_6 <- gam(Pollinators ~ built + trees + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), 
               method="REML",
               data = richness_covariates,
               family=nb())
summary(model_6)
gam.check(model_6)

AIC(model_all, model_null, model_1, model_2, model_3, model_4, model_5, model_6)
# The model with impervious cover explains the most

# top model compared to null
AIC(model_null)-AIC(model_all)

# create a function for predicting the data for plotting purposes
get_prediction_data <- function(formula, data, predictor_name, label) {
  library(mgcv)
  library(dplyr)
  
  model <- gam(formula, method="REML",
               data = data,
               family = nb())
  
  # Create new data for prediction
  newdata <- data.frame(
    built = mean(data$built),
    grass = mean(data$grass),
    trees = mean(data$trees),
    water = mean(data$water),
    log_parksize = mean(data$log_parksize),
    log_obs_pol = mean(data$log_obs_pol),
    lat = mean(data$lat),
    lon = mean(data$lon)
  )
  
  # vary the predictor of interest
  newdata <- newdata[rep(1, 100), ] # repeat 100 rows
  newdata[[predictor_name]] <- seq(min(data[[predictor_name]]), max(data[[predictor_name]]), length.out = 100)
  
  # Predict on link scale with SE
  pred <- predict(model, newdata, type = "link", se.fit = TRUE)
  
  newdata <- newdata %>%
    mutate(
      fit = pred$fit,
      se = pred$se.fit,
      ci_lower = exp(fit - 1.96 * se),
      ci_upper = exp(fit + 1.96 * se),
      Pollinators = exp(fit),
      model = label,
      x = .data[[predictor_name]]
    )
  
  return(newdata %>% dplyr::select(x, Pollinators, ci_lower, ci_upper, model))
}

# Plot the data
pred_built <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), richness_covariates, "built", "Impervious Cover")
pred_grass <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), richness_covariates, "grass", "Grass Cover")
pred_trees <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), richness_covariates, "trees", "Tree Cover")
pred_water <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), richness_covariates, "water", "Water Cover")
pred_parksize <- get_prediction_data(Pollinators ~ built + grass + trees + water + log_parksize + log_obs_pol  + s(lat, lon, k=12, bs="tp"), richness_covariates, "log_parksize", "Log10(Greenspace Size)")

combined_preds <- bind_rows(pred_built, pred_grass, pred_trees, pred_water, pred_parksize)


# adjust the pollinator richness based on number of observations

# Get predicted pollinator richness based on sampling effort only
richness_covariates$pred_null <- predict(model_null, type = "response")

# Reshape richness_covariates to long format for predictor variables
observed_long <- richness_covariates %>%
  dplyr::select(built, grass, trees, water, log_parksize, pred_null) %>%
  pivot_longer(
    cols = c(built, grass, trees, water, log_parksize),
    names_to = "model",
    values_to = "x"
  ) %>%
  mutate(model = case_when(
    model == "built" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.181 SE, <i>P</i> = 0.746</span>",
    model == "grass" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.031 ± 0.031 SE, <i>P</i> = 0.321</span>",
    model == "trees" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.017 ± 0.015 SE, <i>P</i> = 0.266</span>",
    model == "water" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.039 SE, <i>P</i> = 0.878</span>",
    model == "log_parksize" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.010 ± 0.129 SE, <i>P</i> = 0.939</span>"
  ))

combined_preds <- combined_preds %>%
  mutate(model = case_when(
    model == "Impervious Cover" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.181 SE, <i>P</i> = 0.746</span>",
    model == "Grass Cover" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.031 ± 0.031 SE, <i>P</i> = 0.321</span>",
    model == "Tree Cover" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.017 ± 0.015 SE, <i>P</i> = 0.266</span>",
    model == "Water Cover" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.039 SE, <i>P</i> = 0.878</span>",
    model == "Log10(Greenspace Size)" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.010 ± 0.129 SE, <i>P</i> = 0.939</span>"
  ))

# reorder plots so greenspace area is last
facet_levels <- c(
  "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.181 SE, <i>P</i> = 0.746</span>",
  "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.031 ± 0.031 SE, <i>P</i> = 0.321</span>",
  "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.017 ± 0.015 SE, <i>P</i> = 0.266</span>",
  "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.006 ± 0.039 SE, <i>P</i> = 0.878</span>",
  "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.010 ± 0.129 SE, <i>P</i> = 0.939</span>"
)

observed_long <- observed_long %>%
  mutate(model = factor(model, levels = facet_levels))

combined_preds <- combined_preds %>%
  mutate(model = factor(model, levels = facet_levels))

# Now plot using this observed_long for points (x = predictor values, y = adjusted richness residuals)
ggplot() +
  geom_ribbon(data = combined_preds, aes(x = x, ymin = ci_lower, ymax = ci_upper), fill = "lightblue") +
  geom_line(data = combined_preds, aes(x = x, y = Pollinators), color = "blue", size = 1) +
  geom_point(data = observed_long, aes(x = x, y = pred_null), color = "black", alpha = 0.6) +
  facet_wrap(~ model, scales = "free_x", ncol = 2) +
  labs(x = "Predictor Variable",
       y = "Apoidea Species Richness") +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("Figures/Figure_S10.jpeg", height=8, width=8, units="in")

# Create a parameter estimate plot from the full model
# Tidy the model for plotting
model_df <- broom::tidy(model_all, conf.int = TRUE, parametric=TRUE)

# Filter out the intercept if you don't want to include it
model_df <- model_df %>%
  filter(term != "(Intercept)",
         term != "log_obs_pol") %>%
  mutate(term = case_match(term,
                           "built" ~ "Impervious Cover (%)*",
                           "grass" ~ "Grass Cover (%)",
                           "trees" ~ "Tree Cover (%)**",
                           "water" ~ "Water Cover (%)",
                           "log_parksize" ~ "Log10(Greenspace Size)*"))

# Create the parameter estimate plot
ggplot(model_df, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "Coefficient Estimate (± 95% CI)",
    y = ""
  ) +
  theme_classic(base_size = 14)