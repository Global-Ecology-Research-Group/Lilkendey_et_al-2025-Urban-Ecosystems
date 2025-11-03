# In this script, we use iNaturalist data to assess the relationship between pollinators, 
# native and non-native angiosperms, and several environmental covariates in urban greenspaces.

# Read in packages 
library(tidyverse)
library(sf)
library(mgcv)
library(MuMIn)
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

# In this section, we will prepare the data for analysis

# Read in relevant shapefiles - urban areas and parks in urban areas
urban_areas <- st_read("Data/Urban_Areas/Urban_areas.shp")
parks_in_urban <- st_read("Data/ParkServe/parks_in_urban_areas.shp")

# read in florida data of pollinators and angiosperms
data_clean <- readRDS("Data/Pollinator_and_Angiosperm_Data/florida_pollinators_and_angiosperms_taxon.RDS")

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

# Filter this data so that each group has at least 50 observations
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

# Quantify species richness and other measures for modeling ------------------------

# Quantify the number and proportion of observations that are Angiosperms versus Pollinators
cat("Number of pollinator observations:", nrow(cleaned_data[cleaned_data$group=="Pollinators",]))
cat("Percentage of observations that are pollinators:", nrow(cleaned_data[cleaned_data$group=="Pollinators",])/nrow(cleaned_data)*100)

cat("Number of angiosperm observations:", nrow(cleaned_data[cleaned_data$group=="Angiosperms",]))
cat("Percentage of observations that are angiosperms:", nrow(cleaned_data[cleaned_data$group=="Angiosperms",])/nrow(cleaned_data)*100)

# Quantify the total number of Pollinators species present in this study
cat("Total number of pollinator species in this study:", length(unique(cleaned_data[cleaned_data$group=="Pollinators",]$scientific_name)))

# Quantify the total number of angiosperm species present in this study
cat("Total number of angiosperm species in this study:", length(unique(cleaned_data[cleaned_data$group=="Angiosperms",]$scientific_name)))

# Quantify number of pollinator species per park
pollinators <- cleaned_data %>%
  dplyr::filter(group=="Pollinators") %>%
  group_by(ParkID) %>%
  summarise(number_of_species = n_distinct(scientific_name)) %>%
  mutate(group="Pollinators")

summary(pollinators$number_of_species)
sd(pollinators$number_of_species)

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

# save data
write_csv(species_richness_wide, "Data/Pollinator_and_Angiosperm_Data/species_richness_wide_fl.csv")

### All Angiosperms ------------------------

# The data is now ready to model
# Start with modeling the species richness for all angiosperms

# log transform skewed values
species_richness_wide_a <- species_richness_wide %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area)
  ) %>%
  rename(Angiosperms.all=`Angiosperms all`) 

# check for multicolinearity before modeling
vars <- species_richness_wide_a %>%
  dplyr::select(Angiosperms.all, log_angiosperm_obs, log_pollinator_obs, Pollinators)

cor(vars, use = "complete.obs")
# nothing too concerning

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
m_poll <- gam(Pollinators ~ Angiosperms_all.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll)
gam.check(m_poll)

# let's compare this to the null model 
m_poll_null <- gam(Pollinators ~  log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_a, 
              family = nb())
summary(m_poll_null)

AIC(m_poll_null)-AIC(m_poll)

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
             aes(x = Angiosperms_all.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$Pollinators), 75))

ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_a,
             aes(x = Angiosperms_all.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic() 

# what does adding 10 angiosperm species do?
# Calculate differences in predicted pollinator richness
newdat$delta_poll <- c(NA, diff(newdat$fit))
newdat$delta_angio <- c(NA, diff(newdat$Angiosperms_all.fitted))

# Approximate slope: pollinators per 1 angiosperm
newdat$slope <- newdat$delta_poll / newdat$delta_angio

# Average slope across the range
mean_slope <- mean(newdat$slope, na.rm = TRUE)

# Scale to 10 angiosperm species
increase_per_50 <- mean_slope * 50
increase_per_50

### Native Angiosperms ------------------------

# model richness for native angiosperm species

# log transform skewed values
species_richness_wide_n <- species_richness_wide %>%
  filter(native_angiosperm_obs>50) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs),
    log_native_angiosperm_obs = log(native_angiosperm_obs)
  ) %>%
  rename(Angiosperms.native=`Angiosperms Native`)

# check for multicolinearity before modeling
vars <- species_richness_wide_n %>%
  dplyr::select(Angiosperms.native, log_native_angiosperm_obs, log_pollinator_obs, Pollinators)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_n <- gam(Angiosperms.native ~ log_native_angiosperm_obs + s(lat, lon, k=50, bs="tp"), 
               method="REML",
               data = species_richness_wide_n, 
               family = nb())
summary(m_angio_n)
gam.check(m_angio_n)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_n$Angiosperms_native.fitted <- predict(m_angio_n, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_n <- gam(Pollinators ~ Angiosperms_native.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_n, 
              family = nb())
summary(m_poll_n)
gam.check(m_poll_n)

# compare to a null model
m_poll_n_null <- gam(Pollinators ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                method="REML",
                data = species_richness_wide_n, 
                family = nb())
summary(m_poll_n_null)

AIC(m_poll_n_null)-AIC(m_poll_n)

# now visualize the relationship

# Create a new dataset varying Angiosperms_native.fitted
newdat <- with(species_richness_wide_n,
               data.frame(
                 Angiosperms_native.fitted = seq(min(Angiosperms_native.fitted, na.rm=TRUE),
                                              max(Angiosperms_native.fitted, na.rm=TRUE),
                                              length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_n, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_n,
             aes(x = Angiosperms_native.fitted, y = Pollinators),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_native.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_native.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Native Angiosperm Species Richness",
       y = "Pollinator Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_n$Pollinators), 75))

### Non-Native Angiosperms ------------------------

# model the data for non-native species richness

# log transform skewed values
species_richness_wide_nn <- species_richness_wide %>%
  filter(nonnative_angiosperm_obs>50) %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs),
    log_nonnative_angiosperm_obs = log(nonnative_angiosperm_obs)
  ) %>%
  rename(Angiosperms.nonnative=`Angiosperms Non-native`)

# check for multicolinearity before modeling
vars <- species_richness_wide_nn %>%
  dplyr::select(Angiosperms.nonnative, log_nonnative_angiosperm_obs, log_pollinator_obs, Pollinators)

cor(vars, use = "complete.obs")
# nothing too concerning

# use GAM's to model the data
# Model 1: Angiosperm richness ~ sampling effort + spatial location
m_angio_nn <- gam(Angiosperms.nonnative ~ log_nonnative_angiosperm_obs + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_nn, 
                 family = nb())
summary(m_angio_nn)
gam.check(m_angio_nn)

# get the residuals from this model
# Get fitted angiosperm richness
species_richness_wide_nn$Angiosperms_nonnative.fitted <- predict(m_angio_nn, type = "response")

# Model 2: Pollinator richness ~ angiosperm richness + sampling effort + spatial location
m_poll_nn <- gam(Pollinators ~ Angiosperms_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                method="REML",
                data = species_richness_wide_nn, 
                family = nb())
summary(m_poll_nn)
gam.check(m_poll_nn)

# compare to null model
m_poll_nn_null <- gam(Pollinators ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_nn, 
                 family = nb())
summary(m_poll_nn_null)

AIC(m_poll_nn_null)-AIC(m_poll_nn)

model.sel(m_poll_nn_null, m_poll_nn)

# now visualize the relationship

# Create a new dataset varying Angiosperms_native.fitted
newdat <- with(species_richness_wide_nn,
               data.frame(
                 Angiosperms_nonnative.fitted = seq(min(Angiosperms_nonnative.fitted, na.rm=TRUE),
                                                 max(Angiosperms_nonnative.fitted, na.rm=TRUE),
                                                 length.out = 100),
                 log_pollinator_obs = mean(log_pollinator_obs, na.rm=TRUE),
                 lat = mean(lat, na.rm=TRUE),
                 lon = mean(lon, na.rm=TRUE)
               ))

# Predict with standard errors
pred <- predict(m_poll_nn, newdata = newdat, type = "link", se.fit = TRUE)

# Back-transform to response scale
newdat$fit <- exp(pred$fit)   # because NB uses log link
newdat$upper <- exp(pred$fit + 2 * pred$se.fit)
newdat$lower <- exp(pred$fit - 2 * pred$se.fit)

# Plot
ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nn,
             aes(x = Angiosperms_nonnative.fitted, y = Pollinators),
             alpha = 0.7, color = "black") +
  
  # Partial effect line
  geom_line(data = newdat,
            aes(x = Angiosperms_nonnative.fitted, y = fit),
            color = "blue", size = 1) +
  
  # Confidence ribbon
  geom_ribbon(data = newdat,
              aes(x = Angiosperms_nonnative.fitted, ymin = lower, ymax = upper),
              alpha = 0.2, fill = "lightblue") +
  
  labs(x = "Non-Native Angiosperm Species Richness",
       y = "Pollinator Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_nn$Pollinators), 75))

### Proportion of Native Angiosperms ------------------------

# model the data for proportion of native angiosperm species

# log transform skewed values
species_richness_wide_np <- species_richness_wide %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_native = `Angiosperms Native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_np$prop_native)

# check for multicolinearity before modeling
vars <- species_richness_wide_np %>%
  dplyr::select(prop_native, log_angiosperm_obs, log_pollinator_obs, Pollinators)

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
m_poll_np <- gam(Pollinators ~ prop_native.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
              method="REML",
              data = species_richness_wide_np, 
              family = nb())
summary(m_poll_np)
gam.check(m_poll_np)

# compare to a null model
m_poll_np_null <- gam(Pollinators ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
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
             aes(x = prop_native.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_np$Pollinators), 75))

ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_np,
             aes(x = prop_native.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic()

### Proportion of Non-Native Angiosperms ------------------------

# model the data for proportion of non-native angiosperm species

# log transform skewed values
species_richness_wide_nnp <- species_richness_wide %>%
  mutate(
    log_pollinator_obs = log(pollinator_obs),
    log_angiosperm_obs = log(angiosperm_obs),
    urban_area = factor(urban_area),
    prop_nonnative = `Angiosperms Non-native`/`Angiosperms all`
  ) 

# examine the distribution of proportion of native species
hist(species_richness_wide_nnp$prop_nonnative)

# check for multicolinearity before modeling
vars <- species_richness_wide_nnp %>%
  dplyr::select(prop_nonnative, log_angiosperm_obs, log_pollinator_obs, Pollinators)

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
m_poll_nnp <- gam(Pollinators ~ prop_nonnative.fitted + log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
                 method="REML",
                 data = species_richness_wide_nnp, 
                 family = nb())
summary(m_poll_nnp)
gam.check(m_poll_nnp)

# compare to a null model
m_poll_nnp_null <- gam(Pollinators ~ log_pollinator_obs + s(lat, lon, k=50, bs="tp"), 
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
             aes(x = prop_nonnative.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic() +
  coord_cartesian(ylim = c(min(species_richness_wide_a$Pollinators), 75))

ggplot() +
  # Raw data scatter
  geom_point(data = species_richness_wide_nnp,
             aes(x = prop_nonnative.fitted, y = Pollinators),
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
       y = "Pollinator Species Richness") +
  theme_classic() 

### Create a supplementary figure based on number of observations ------------------------

# create a supplmentry figure to show the relationship between species richness and number of observations

# Plot angiosperm observations and richness
(p1 <- ggplot(species_richness_wide_a, aes(x = log_angiosperm_obs, y = Angiosperms.all)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), 
                method.args = list(family = "nb"), 
                color = "blue", fill = "lightblue") +
    labs(
      x = "Log(Number of Angiosperm Observations)",
      y = "Angiosperm Species Richness"
    ) +
    theme_classic(base_size = 14))

# Plot pollinator observations and richness
(p2 <- ggplot(species_richness_wide_a, aes(x = log_pollinator_obs, y = Pollinators)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp"), 
                method.args = list(family = "nb"), 
                color = "blue", fill = "lightblue") +
    labs(
      x = "Log(Number of Pollinator Observations)",
      y = "Pollinator Species Richness",
    ) +
    theme_classic(base_size = 14))

# Combine and label
combined_plot <- (p1 / p2) +
  plot_annotation(tag_levels = 'a') &
  theme(
    plot.tag = element_text(size = 16, face = "bold")
  )

# Show combined plot
combined_plot

ggsave("Figures/Figure_S1.jpeg", height=7, width=7, units="in")

# Model pollinator richness with landcover attributes ------------------------

## Adding landcover values ----------------------------------------------------

# Read in data from dynamic world
landcover <- read_csv("Data/Dynamic_World_Habitat/dynamic_world_habitat_fl.csv")

# add this data to our data frame
lc_data <- species_richness_wide %>%
  filter(complete.cases(Pollinators)) %>%
  left_join(landcover, by = "ParkID") %>%
  group_by(ParkID) %>%
  summarise(Shape_Area=sum(Shape_Area), trees=mean(trees), snow_and_ice=mean(snow_and_ice),
            water=mean(water), shrub_and_scrub=mean(shrub_and_scrub), built=mean(built),
            crops=mean(crops), grass=mean(grass), flooded_vegetation=mean(flooded_vegetation),
            bare=mean(bare), pollinator_obs=sum(pollinator_obs), Pollinators=sum(Pollinators),
            lon=first(lon), lat=first(lat))

# let's examine the landcover data to figure out which should be included in the model
summary(lc_data[,c("trees", "snow_and_ice", "water", "shrub_and_scrub", "built", 
                   "bare", "crops", "grass", "flooded_vegetation")])
# the only variables with a median > 0 are trees, water, built, and grass so we will use these variables
vars <- c("trees", "snow_and_ice", "water", "shrub_and_scrub", 
          "built", "bare", "crops", "grass", "flooded_vegetation")

summary_df <- data.frame(
  variables = vars,
  median = sapply(lc_data[, vars], median, na.rm = TRUE),
  mean = sapply(lc_data[, vars], mean, na.rm = TRUE),
  sd   = sapply(lc_data[, vars], sd,   na.rm = TRUE),
  num_of_parks_w_lc = sapply(lc_data[, vars], function(x) sum(x > 0, na.rm = TRUE))
)

summary_df

# the features that are most common in parks are built, grass, trees, and water so we will model with those
# let's examine their distribution
hist(lc_data$built)
hist(lc_data$grass)
hist(lc_data$trees)
hist(lc_data$water)

## Model the data ------------------------

# Plots Pollinator richness, All Angiosperm richness, Native Angiosperm 
# richness, and Non-native Angiosperm richness against park size, impervious cover percent, 
# non-tree vegetation cover percent, water cover percent, number of observations,
# and population density

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
    model == "built" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.002 SE, <i>P</i> = 0.010</span>",
    model == "grass" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.002 SE, <i>P</i> = 0.391</span>",
    model == "trees" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.001 SE, <i>P</i> = 0.001</span>",
    model == "water" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.00005 ± 0.002 SE, <i>P</i> = 0.985</span>",
    model == "log_parksize" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.028 ± 0.013 SE, <i>P</i> = 0.037</span>"
  ))

combined_preds <- combined_preds %>%
  mutate(model = case_when(
    model == "Impervious Cover" ~ "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.002 SE, <i>P</i> = 0.010</span>",
    model == "Grass Cover" ~ "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.002 SE, <i>P</i> = 0.391</span>",
    model == "Tree Cover" ~ "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.001 SE, <i>P</i> = 0.001</span>",
    model == "Water Cover" ~ "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.00005 ± 0.002 SE, <i>P</i> = 0.985</span>",
    model == "Log10(Greenspace Size)" ~ "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.028 ± 0.013 SE, <i>P</i> = 0.037</span>"
  ))

# reorder plots so greenspace area is last
facet_levels <- c(
  "Impervious Surface Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.004 ± 0.002 SE, <i>P</i> = 0.010</span>",
  "Grass Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.002 ± 0.002 SE, <i>P</i> = 0.391</span>",
  "Tree Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.005 ± 0.001 SE, <i>P</i> = 0.001</span>",
  "Water Cover<br><span style='font-size:12pt; font-weight:normal;'>Estimate = -0.00005 ± 0.002 SE, <i>P</i> = 0.985</span>",
  "Log10(Greenspace Size)<br><span style='font-size:12pt; font-weight:normal;'>Estimate = 0.028 ± 0.013 SE, <i>P</i> = 0.037</span>"
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
       y = "Pollinator Species Richness") +
  theme_classic(base_size = 14) +
  theme(
    strip.text = element_markdown(size = 16, face = "bold"),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

ggsave("Figures/Figure_3.jpeg", height=8, width=8, units="in")

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
  theme_classic(base_size = 17)
