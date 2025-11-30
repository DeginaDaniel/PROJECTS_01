# ---- Modélisation MaxEnt (MaxNet) pour le Spillover de Mpox ---- #

# Charger les bibliothèques nécessaires
library(raster)   # Manipulation de rasters spatiaux
library(maxnet)   # Modélisation MaxEnt
library(tidyr)    # Manipulation des données (séparation de colonnes)
library(dplyr)    # Manipulation et filtrage de données

# ----------------------------
# Définition du répertoire contenant les inputs
road <- "D:/SpillOver Mpox/files"

# ----------------------------
# Lecture des données d'occurrences de l'espèce
presences_raw <- read.csv2(file.path(road, "mpox_species.csv"), stringsAsFactors = FALSE)

# Vérification et séparation des colonnes si nécessaire
if(ncol(presences_raw) == 1){
  # Cas où toutes les informations sont dans une seule colonne séparée par ";"
  presences <- presences_raw %>%
    separate(col = 1, into = c("genus","species","lat","lon"), sep = ";", remove = TRUE)
} else {
  presences <- presences_raw
}

# Maintien des coordonnées en numérique
presences$lat <- as.numeric(presences$lat)
presences$lon <- as.numeric(presences$lon)

# ----------------------------
# Lecture des rasters background points (environnement data)
r1 <- raster(file.path(road, "DEM_RDC.tif"))                     # Modèle numérique de terrain
r2 <- raster(file.path(road, "EVI_DRC_2022.tif"))                # Indice de végétation EVI
r3 <- raster(file.path(road, "MODIS_MCD12Q1_LandCover_2022_DRC_500m.tif"))  # Couverture terrestre
r4 <- raster(file.path(road, "NDVI_MeanNDVI_DRC_2022.tif"))      # NDVI moyen
r5 <- raster(file.path(road, "precipitation_mean_2010-2021_2.tif"))
r6 <- raster(file.path(road, "temperature_mean_2010-2021_2.tif"))

# ----------------------------
# Nettoyage des données d'occurrences
presences <- presences %>% filter(!is.na(lat) & !is.na(lon) & !is.na(species))

# ----------------------------
# Préparation du RasterStack pour la modélisation
# Utilisation de r2 comme référence spatiale pour aligner raster r1
ref <- r2

# Reprojection et resampling de tous les rasters pour qu'ils aient la même étendue et résolution
r1_fix <- projectRaster(r1, ref, method = "bilinear")
r3_fix <- projectRaster(r3, ref, method = "ngb") # facteur/catégoriel → nearest neighbor (method appropriée pour interpolation des données catégorielles car conserve les classes)
r4_fix <- projectRaster(r4, ref, method = "bilinear")
r5_fix <- projectRaster(r5, ref, method = "bilinear")
r6_fix <- projectRaster(r6, ref, method = "bilinear")

# Création du RasterStack aligné
myraster <- stack(r1_fix, r2, r3_fix, r4_fix, r5_fix, r6_fix)

# Conversion de la variable catégorielle LandCover en facteur
myraster[[3]] <- as.factor(myraster[[3]])

# ----------------------------
# Définition du répertoire des outputs
parent_dir <- dirname(road)
outdir <- file.path(parent_dir, "MaxEnt_results")
dir.create(outdir, showWarnings = FALSE)

# ----------------------------
# Boucle de modélisation MaxEnt pour chaque espèce
species_list <- unique(presences$species)

for(sp in species_list){
  
  sp_data <- presences[presences$species == sp, ]
  
  # Vérifier le nombre minimal de points d'occurrence pour modéliser
  if(nrow(sp_data) < 8){
    cat("Attention :", sp, "a moins de 8 occurrences. Modèle non réalisé.\n")
    next
  }
  
  cat("Traitement de l'espèce :", sp, "2\n")
  
  # Extraction des coordonnées
  coords <- sp_data[, c("lon", "lat")]
  
  # Extraction des background points pour les présences
  pres_env <- raster::extract(myraster, coords)
  pres_labels <- rep(1, nrow(pres_env))  # Étiquettes de présence = 1
  
  # ----------------------------
  # Génération de pseudo-absences aléatoires
  n_bg <- 10000
  bg_cells <- sample(ncell(myraster), n_bg)
  bg_coords <- xyFromCell(myraster, bg_cells)
  bg_env <- raster::extract(myraster, bg_coords)
  bg_labels <- rep(0, nrow(bg_env))  # Étiquettes d'absence = 0
  
  # Combinaison des données de présence et de pseudo-absence
  train_data <- rbind(pres_env, bg_env)
  train_labels <- c(pres_labels, bg_labels)
  train_data <- as.data.frame(train_data)
  
  # Nettoyage : suppression des lignes avec valeurs manquantes
  valid <- complete.cases(train_data)
  train_data <- train_data[valid, ]
  train_labels <- train_labels[valid]
  
  # Suppression des colonnes constantes (sans variation)
  train_data <- train_data[, sapply(train_data, function(x) length(unique(x)) > 1)]
  
  # Vérification qu'il reste des variables informatives
  if(ncol(train_data) == 0){
    cat("Aucune variable informative pour", sp, ". Modèle non réalisé.\n")
    next
  }
  
  # ----------------------------
  # Ajustement du modèle MaxEnt
  mx_model <- maxnet(p = train_labels, data = train_data,
                     f = maxnet.formula(train_labels, train_data))
  
  # Prédiction sur l'ensemble du RasterStack
  pred_raster <- raster::predict(object = myraster, model = mx_model, type = "cloglog", progress = "text")
  
  # Sauvegarde du raster de prédiction (output)
  out_file <- file.path(outdir, paste0(sp, "_pred2.tif"))
  writeRaster(pred_raster, filename = out_file, format = "GTiff", overwrite = TRUE)
  
  cat("Terminé pour :", sp, "2\n")
}