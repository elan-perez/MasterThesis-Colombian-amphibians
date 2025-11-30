###############################################################################
# MASTER'S THESIS SCRIPT: Colombian Amphibian Habitat & Conservation Analysis
# Author: Elán Pérez García
# Date: 30/11/2025
#
# Purpose:
# Full reproducible pipeline that:
# 1. Imports IUCN Colombian amphibian range shapefiles
# 2. Queries IUCN Red List for habitat & threat information
# 3. Cleans and formats habitat codes (NewCode)
# 4. Builds habitat-filtered species rasters
# 5. Add species not in IUCN through coordinates 
# 6. Produces richness rasters (for all Colombian species & threatened)
# 7. Assesses representation in protected areas (WDPA) and unprotected species
# 8. Compute Corrected Weighted Endemism
# 
#
# Notes:
# - Requires IUCN Red List API key in environment variable IUCN_REDLIST_KEY
# - Designed to run sequentially; sections are numbered for clarity
# - Sources of input data:
#   * IUCN amphibian range shapefiles: downloaded from IUCN Red List
#   * Habitat raster: IUCN Habitat Classification v2, composite level 2
#   * Colombia boundaries: rnaturalearth (Natural Earth) & GADM shapefiles
#   * Protected areas: WDPA May 2025
###############################################################################


# -------------------------------
# 0. Setup: working dir + packages
# -------------------------------
wd <- "~/Desktop/Thesis" # <-- # Change to your project folder
setwd(wd)

required_pkgs <- c("terra", "sf", "rnaturalearth", "rredlist", "dplyr", "exactextractr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if(length(missing_pkgs)) stop("Missing packages: ", paste(missing_pkgs, collapse=", "))

# Load libraries (terra preferred for raster/vector processing)
library(terra)
library(sf)
library(rnaturalearth)
library(rredlist)
library(dplyr)
library(exactextractr)

# Create directories used by the pipeline
dir.create(file.path(wd, "tmp"), showWarnings = FALSE)


# -------------------------------
# 1. Species range preparation
# -------------------------------
cat("\n--- 1. Species range preparation ---\n")

# Load IUCN amphibian range shapefiles
amp1 <- vect("AMPHIBIANS_PART1.shp")
amp2 <- vect("AMPHIBIANS_PART2.shp")

# Colombia boundary (Natural Earth)
colombia <- ne_countries(country = "Colombia", returnclass = "sf") |> vect()

# Reproject the range layers to Colombia CRS and intersect to crop
amp1_proj <- project(amp1, colombia)
amp2_proj <- project(amp2, colombia)
amp1_col <- intersect(amp1_proj, colombia)
amp2_col <- intersect(amp2_proj, colombia)
IUCN <- rbind(amp1_col, amp2_col) # IUCN species
# Get Amphibian Species of the World (ASW) endemic or natural resident species of Colombia
ASW <- read.csv("ASW_colombia_species.csv")
# Combine ASW and IUCN and match Colombian resident or endemic species found in both
common_species <- intersect(IUCN$species, ASW$species)
amphibians_native <- IUCN[IUCN$species %in% common_species, ]
# For the species found in both
writeVector(amphibians_native, "amphibians_colombia.shp", overwrite = TRUE)
cat("Saved amphibians_colombia.shp with", length(unique(amphibians_native$sci_name)), "unique species\n")


# -------------------------------
# 2. Habitat & threat extraction
# -------------------------------
cat("\n--- 2. Habitat & threat extraction ---\n")

# Read habitat code legend (IUCN mapping)
habitat.code <- read.csv("IUCN_mapping_legend.csv", fileEncoding = "latin1", stringsAsFactors = FALSE)
codes <- strsplit(habitat.code[,1], " ")
out <- do.call(rbind, lapply(seq_len(nrow(habitat.code)), function(i) c(code = codes[[i]][1], NewCode = habitat.code$NewCode[i])))
out <- as.data.frame(out, stringsAsFactors = FALSE)

# Species list 
sp_raw <- read.csv("amphibians_colombia.csv")
species_list <- unique(sp_raw$sci_name)

# Helper: query IUCN Red List (habitats and threat info)
get_rl_data <- function(sci_name){
  parts <- strsplit(sci_name, " ")[[1]]
  if (length(parts) != 2) return(NULL)
  genus <- parts[1]; species <- parts[2]
  tryCatch(rl_species_latest(genus = genus, species = species, key = Sys.getenv("IUCN_REDLIST_KEY")), error = function(e) NULL)
}

# Fetch habitats for each species
hab_list <- list(); miss_hab <- character()
for (s in species_list){
  cat("Query habitat:", s, "\n")
  dat <- get_rl_data(s)
  if (!is.null(dat) && !is.null(dat$habitats) && nrow(dat$habitats) > 0) {
    tmp <- dat$habitats
    tmp$species <- s
    hab_list[[length(hab_list) + 1]] <- tmp
  } else {
    miss_hab <- c(miss_hab, s)
  }
}

habitats_df <- bind_rows(hab_list)
if ("description" %in% names(habitats_df)) habitats_df$description <- sapply(habitats_df$description, function(x) if (is.list(x) && !is.null(x$en)) x$en else as.character(x))
habitats_merged <- merge(habitats_df, out, by = "code", all.x = TRUE)

# Save intermediate files
write.csv(habitats_df, "Species_Habitats_raw.csv", row.names = FALSE)
write.csv(habitats_merged, "Species_Habitats_withCodes.csv", row.names = FALSE)
write.csv(miss_hab, "Species_missing_Habitats.csv", row.names = FALSE)

# Fetch threat (red list category) for each species
threat_list <- list(); miss_threat <- character()
for (s in species_list){
  cat("Query threat:", s, "\n")
  dat <- get_rl_data(s)
  if (!is.null(dat) && !is.null(dat$red_list_category)){
    tmp <- dat$red_list_category
    tmp$species <- s
    threat_list[[length(threat_list) + 1]] <- tmp
  } else {
    miss_threat <- c(miss_threat, s)
  }
}

th <- bind_rows(threat_list)
# Clean threat structure to simple table
th_clean <- data.frame(
  version = sapply(th$version, `[`, 1),
  description = sapply(th$description, function(x) { if (is.list(x[[1]]) && !is.null(x[[1]]$en)) x[[1]]$en else as.character(x[[1]]) }),
  code = sapply(th$code, `[`, 1),
  species = th$species,
  stringsAsFactors = FALSE
)
th_clean$threatened <- th_clean$code %in% c("CR", "EN", "VU")

write.csv(th_clean, "Species_ThreatStatus.csv", row.names = FALSE)
write.csv(miss_threat, "Species_missing_ThreatStatus.csv", row.names = FALSE)


# -------------------------------
# 3. Fix NewCode formatting and merge
# -------------------------------
cat("\n--- 3. Fix NewCode formatting ---\n")

# Load merged habitats file
hab <- read.csv("Species_Habitats_withCodes.csv", stringsAsFactors = FALSE)
unique_codes <- unique(hab$code)

new_codes <- sapply(unique_codes, function(x){
  x <- as.character(x)
  x_clean <- gsub("_", "0", x)
  if (!grepl("^[0-9]+$", x_clean)) return(NA)
  as.numeric(sprintf("%03d", as.numeric(x_clean)))
})

lookup <- data.frame(code = unique_codes, NewCode = new_codes, stringsAsFactors = FALSE)
hab_clean <- merge(hab, lookup, by = "code", all.x = TRUE)
# tidy column names
if ("NewCode.x" %in% names(hab_clean)) hab_clean$NewCode.x <- NULL
if ("NewCode.y" %in% names(hab_clean)) names(hab_clean)[names(hab_clean) == "NewCode.y"] <- "NewCode"
write.csv(hab_clean, "Species_Habitats_with_NewCode.csv", row.names = FALSE)


# -------------------------------
# 4. Build habitat-filtered species rasters
# -------------------------------
cat("\n--- 4. Build habitat-filtered species rasters ---\n")

# Load vectors + habitat raster (cropped version will be created)
sp_polygons <- vect("amphibians_colombia.shp")

# Habitat classification raster: IUCN Habitat Classification v2, composite level 2
habitat_full <- rast("iucn_habitatclassification_composite_lvl2_ver004.tif")

# Colombia boundary from GADM (level 1) for masking habitat raster
colombia_vect <- vect("gadm41_COL_1.shp")

# Crop and mask habitat raster to Colombia 
habitat_colombia <- tryCatch({
  crop(habitat_full, colombia_vect) %>% mask(colombia_vect)
}, error = function(e){
  # if pre-cropped file exists, try to load
  if (file.exists("habitat_colombia.tif")) {
    rast("habitat_colombia.tif")
  } else stop("Could not crop habitat raster and no pre-cropped file found.")
})
writeRaster(habitat_colombia, "habitat_colombia.tif", overwrite = TRUE)

# Load hab_data with NewCode
hab_data <- read.csv("Species_Habitats_with_NewCode.csv", stringsAsFactors = FALSE)
species_list2 <- unique(hab_data$species)

# tmp dir path
tmp_dir <- file.path(wd, "tmp")

create_species_raster <- function(sp_name){
  cat("Processing species:", sp_name, "\n")
  sp_hab <- subset(hab_data, species == sp_name)
  valid_codes <- as.numeric(unique(sp_hab$NewCode))
  if (length(valid_codes) == 0 || all(is.na(valid_codes))){
    cat(" -> No valid habitat codes for", sp_name, "\n")
    return(NULL)
  }
  rcl <- cbind(valid_codes, 1)
  sp_r <- classify(habitat_colombia, rcl = rcl, others = NA)
  sp_shape <- sp_polygons[sp_polygons$sci_name == sp_name, ]
  if (nrow(sp_shape) == 0) { cat(" -> No polygon for", sp_name, "\n"); return(NULL) }
  sp_shape <- intersect(sp_shape, colombia_vect)
  if (nrow(sp_shape) == 0) { cat(" -> No overlap Colombia for", sp_name, "\n"); return(NULL) }
  sp_masked <- mask(sp_r, sp_shape)
  sp_masked[is.na(sp_masked)] <- 0
  outfile <- file.path(tmp_dir, paste0(sp_name, ".tif"))
  writeRaster(sp_masked, outfile, overwrite = TRUE)
  return(outfile)
}
for (s in species_list2) create_species_raster(s)


# -------------------------------
# 5. Species added through coordinates
# -------------------------------
cat("\n--- 6. Species added through coordinates ---\n")

# Some species present in Colombia (based on ASW/literature) were not included
# in the IUCN polygon data. These species were added manually using presence
# coordinates and converted into raster layers to be included in the analyses.

# Add them to the folder with other species rasters
setwd (tmp.dir)

# Example: Atelopus calima (single known locality with 20 km buffer)
coords <- data.frame(lon = -76.4167, lat = 3.8667)
# Convert coordinate to spatial object
point_sf <- st_as_sf(coords, coords = c("lon", "lat"), crs = 4326)
# Project to metric CRS (for buffering)
point_proj <- st_transform(point_sf, 3116)
# Create 20 km buffer
buffer_proj <- st_buffer(point_proj, dist = 20000)
# Reproject back to WGS84
buffer_wgs <- st_transform(buffer_proj, 4326)
# Load Colombia habitat raster as template
template_raster <- rast("habitat_colombia.tif")
# Rasterise buffer to match template
buffer_rast <- rasterize(vect(buffer_wgs), template_raster, field = 1, background = 0)
# Save raster to new species folder
writeRaster(buffer_rast,
            file.path(new_dir, "Atelopus_calima.tif"),
            datatype = "FLT4S",
            overwrite = TRUE)


# -------------------------------
# 6. Build richness rasters (all & threatened)
# -------------------------------
cat("\n--- 5. Build richness rasters ---\n")

build_richness <- function(species_vector, output){
  cat("Building richness ->", output, "\n")
  base <- setValues(habitat_colombia, 0)
  for (sp in species_vector){
    f <- file.path(tmp_dir, paste0(sp, ".tif"))
    if (!file.exists(f)) { cat("Missing raster for", sp, "\n"); next }
    r <- rast(f)
    base <- base + r
  }
  writeRaster(base, output, overwrite = TRUE)
  cat("Saved", output, "\n")
}

all_species <- unique(read.csv("Species.inc.Threat.rest.csv", stringsAsFactors = FALSE)$species)
build_richness(all_species, "Richness.habitats.neu.tif")

threat_info <- read.csv("Species.inc.Threat.rest.csv", stringsAsFactors = FALSE)
threatened_species <- subset(threat_info, threatened == TRUE)$species
build_richness(threatened_species, "Richness.threat.habitats.neu.tif")

# -------------------------------
# 7. WDPA representation analysis & unprotected species
# -------------------------------
cat("\n--- 7. WDPA representation analysis ---\n")

# WDPA polygons (protected areas)
# Source: World Database on Protected Areas (May 2025 version)
wdpa <- vect("WDPA_WDOECM_May2025_Public_COL_shp-polygons.shp")
setwd(file.path(wd, "tmp"))
species_files <- list.files(getwd(), pattern = "\\.tif$", full.names = TRUE)
setwd(wd)

# matrix to store presence/absence per polygon
results <- matrix(0, nrow = length(wdpa), ncol = length(species_files))
colnames(results) <- basename(species_files)
for (i in seq_along(species_files)){
  cat("Extracting WDPA for file", i, "of", length(species_files), "\n")
  r <- rast(species_files[i])
  e <- extract(r, wdpa, fun = "max", touches = TRUE, na.rm = TRUE)
  results[, i] <- e[,2]
}

N.sp <- rowSums(results > 0, na.rm = TRUE)
wdpa$N_species <- N.sp
results_df <- as.data.frame(results)
results_df$WDPAID <- wdpa$WDPAID
results_df$N_species <- N.sp
write.csv(results_df, "summary.wdpa.neu.csv", row.names = FALSE)
writeVector(wdpa, "WDPA_with_species_richness.gpkg", overwrite = TRUE)

# Unprotected species (absent from all WDPA)
res_clean <- results_df[, !(names(results_df) %in% c("WDPAID", "N_species"))]
unprotected_counts <- colSums(res_clean, na.rm = TRUE)
unprotected <- names(unprotected_counts[unprotected_counts == 0])
# convert to clean species names
unprotected_names <- gsub("\\.tif$", "", unprotected)
unprotected_names <- gsub("\\.", " ", unprotected_names)
write.csv(unprotected_names, "unprotected.sp.neu.csv", row.names = FALSE)

# create richness raster for unprotected species
mask <- setValues(habitat_colombia, 0)
for (f in unprotected){
  fp <- file.path(tmp_dir, f)
  if (file.exists(fp)) mask <- mask + rast(fp)
}
writeRaster(mask, "Richness.unprotected.sp.habitats.neu.tif", overwrite = TRUE)

# Threatened species counts inside WDPA
thred <- read.csv("Species_ThreatStatus.csv", stringsAsFactors = FALSE)
thred <- thred[, c("species", "threatened")]
names(thred) <- c("Species", "Threatened")
res <- read.csv("summary.wdpa.neu.csv", stringsAsFactors = FALSE)
n.sp <- res$N_species
rownames(res) <- res$WDPAID
res <- res[, !(names(res) %in% c("WDPAID", "N_species"))]
res_t <- as.data.frame(t(res))
res_spnames <- rownames(res_t)
res_spnames <- gsub("\\.tif$", "", res_spnames)
res_spnames <- gsub("\\.", " ", res_spnames)
res_t$res.sp <- res_spnames

out4 <- read.csv("Species.inc.Habitats.rest.csv", stringsAsFactors = FALSE)
out4$Species <- gsub("\\.tif$", "", out4$species)
out4$Species <- gsub("\\.", " ", out4$Species)
out4 <- merge(out4, thred, by = "Species")
out4 <- subset(out4, Threatened == TRUE)
sp_threat <- unique(out4$Species)
res_threat <- res_t[res_t$res.sp %in% sp_threat, ]
res_threat <- as.data.frame(t(res_threat[,-ncol(res_threat)]))
res_threat$N_species.threat <- rowSums(res_threat)
res_threat$N_species <- n.sp

wdpa2 <- vect("WDPA_WDOECM_May2025_Public_COL_shp-polygons.shp")
res_threat$WDPAID <- wdpa2$WDPAID
write.csv(res_threat, "summary.wdpa.threat.csv", row.names = FALSE)
wdpa_sp <- merge(wdpa2, res_threat, by = "WDPAID")
writeVector(wdpa_sp, "WDPA.inc.threat.sp.numbers.gpkg", overwrite = TRUE)


# -------------------------------
# 8. Corrected Weighted Endemism (CWE)
# -------------------------------
cat("\n--- 8. Compute CWE ---\n")

# Species area & weighted raster
out <- data.frame()
for(i in 1:length(species_list)){
  sp <- species_list[i]
  matches <- sp_polygons[sp_polygons$sci_name==sp, ]
  if(nrow(matches)>0){ area <- sum(expanse(matches, unit="km")); out <- rbind(out, data.frame(species=sp, area_km2=area)) }
}
write.csv(out, "tmp.csv", row.names=FALSE)

mask <- colombia; mask <- setValues(mask,0); writeRaster(mask,"mask.tif",overwrite=TRUE)
mask <- rast("mask.tif")
mask.2 <- mask

for(i in 1:nrow(out)){
  species_name <- out$species[i]; area <- out$area_km2[i]; raster_path <- paste0(wd,"/tmp/",species_name,".tif")
  if(file.exists(raster_path)){
    r.i <- rast(raster_path); r.i[r.i==1] <- 1/area; mask.2 <- mask.2 + r.i
  } else warning(paste("Missing raster for:",species_name))
}
writeRaster(mask.2,"we.2.tif",overwrite=TRUE)
removeTmpFiles(h=0.1)
mask.3 <- rast("we.2.tif")
mask.4 <- rast("Richness.habitats.neu.tif")
mask.3[is.na(mask.4)] <- NA
mask.4[is.na(mask.3)] <- NA
cwe <- mask.3/mask.4
writeRaster(cwe,"corrected_weighted_endemism.tif",overwrite=TRUE)

cat("\n=== MASTER SCRIPT COMPLETED ===\n")

