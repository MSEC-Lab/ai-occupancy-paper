########### Expert vs AI psi posterior projections Example

################# projecting psi from models in Guatemala per-species outputs:

# Load Packages:

suppressPackageStartupMessages({
  library(camtrapR)
  library(rjags)
  library(coda)
  library(terra)
  suppressWarnings(try(library(nimble), silent = TRUE))
})

# Set Directory:

base_dir   <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters"

# Folder with your GT *_fit.rds and *_mod.rds
models_dir <- file.path(base_dir, "_GT_mod_and_fit_files")

# Scaled 300 m WGS84 covariate stack for GT
stack_fp   <- file.path(
  base_dir,
  "_GT_covariates_WGS84_300m_TRIM_SCALED",
  "GT_covariates_WGS84_300m_STACK_TRIM_SCALED.tif"
)

# Where to write tiled projections + mosaics
out_root   <- file.path(base_dir, "_GT_projections_TILES_FINAL")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

# =================== SETTINGS ===================
draws    <- 100            # MCMC draws per tile
nx       <- 3              # tiles horizontally
ny       <- 3              # tiles vertically
datatype <- "FLT4S"
gdal_opt <- c("COMPRESS=DEFLATE", "PREDICTOR=3", "ZLEVEL=6")

# Superset of GT covariates prepared (models will only use what they need)
req_covs_superset <- c("d2hs","d2hs_squared","elevation","slope","d2w","canopy","evi","precipitation")

terraOptions(memfrac = 0.8)

# =================== LOAD STACK ===================
stopifnot(file.exists(stack_fp))
stk_all <- rast(stack_fp)

# Keep only the covariates we actually have (intersection)
req_covs <- intersect(req_covs_superset, names(stk_all))
if (!length(req_covs)) {
  stop("Scaled stack has none of the expected covariates.\nHave: ",
       paste(names(stk_all), collapse=", "),
       "\nExpected any of: ",
       paste(req_covs_superset, collapse=", "))
}
stk <- stk_all[[req_covs]]
message("Using covariates in stack: ", paste(req_covs, collapse=", "))

# =================== HELPERS ===================

# 1) Tile extents
mk_tile_extents <- function(x, nx = 2, ny = 2) {
  e  <- ext(x)
  xs <- seq(e$xmin, e$xmax, length.out = nx + 1)
  ys <- seq(e$ymin, e$ymax, length.out = ny + 1)
  exts <- vector("list", nx * ny)
  k <- 1
  for (j in 1:ny) for (i in 1:nx) {
    exts[[k]] <- ext(xs[i], xs[i+1], ys[j], ys[j+1])
    k <- k + 1
  }
  exts
}

# 2) Robust field extractor that works across list/S3/S4/env WITHOUT isS4()
get_field <- function(obj, name) {
  # list/data.frame
  if (is.list(obj) && !is.null(obj[[name]])) return(obj[[name]])
  # environment
  if (is.environment(obj) && exists(name, obj, inherits = FALSE)) return(get(name, obj))
  # S4-ish: safely probe slotNames(); if it errors, it's not S4
  sn_try <- try(methods::slotNames(obj), silent = TRUE)
  if (!inherits(sn_try, "try-error")) {
    if (length(sn_try) && name %in% sn_try) {
      return(methods::slot(obj, name))
    }
  }
  # Try generic [[ ]] (some S3)
  out <- try(obj[[name]], silent = TRUE)
  if (!inherits(out, "try-error")) return(out)
  NULL
}

# 3) Try to detect which covariates the model expects; otherwise return full stack
subset_stack_for_model <- function(stk, mod) {
  candidate_fields <- c(
    "covariates","covariate_names","x_names","X_names",
    "predictor_names","predictors","features",
    # common nested places
    "data","dat","model","settings","params","X","x","covX"
  )

  found <- character(0)

  for (nm in candidate_fields) {
    v <- get_field(mod, nm)
    if (is.null(v)) next
    if (is.matrix(v) || is.data.frame(v)) {
      cn <- colnames(v); if (length(cn)) found <- c(found, cn)
    } else if (is.character(v)) {
      found <- c(found, v)
    } else if (is.list(v)) {
      for (w in v) {
        if (is.matrix(w) || is.data.frame(w)) {
          cn <- colnames(w); if (length(cn)) found <- c(found, cn)
        } else if (is.character(w)) {
          found <- c(found, w)
        }
      }
    }
  }

  found <- unique(as.character(found))
  keep  <- intersect(names(stk), found)

  if (length(keep)) {
    message("  • Model covariates detected: ", paste(keep, collapse=", "))
    stk[[keep]]
  } else {
    message("  • Could not detect model covariate names; using full prepared stack: ",
            paste(names(stk), collapse=", "))
    stk
  }
}

# 4) Predict wrapper
safe_predict_tile <- function(mod, fit, tile, draws) {
  preds <- try(predict(object=mod, mcmc.list=fit, x=tile,
                       type="psi", batch=TRUE, draws=draws), silent=TRUE)
  if (inherits(preds, "try-error")) {
    stop("predict() failed. Tile covariates: ", paste(names(tile), collapse=", "),
         "\nUnderlying error: ", conditionMessage(attr(preds, "condition")))
  }
  preds
}

# 5) Tile writer (per model)
write_tiles_for_model <- function(fit_fp, mod_fp, stk, tile_exts, out_dir,
                                  draws = 100, datatype = "FLT4S", gdal_opt = NULL) {
  message("\n=== ", basename(fit_fp), " ===")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  tile_dir <- file.path(out_dir, "_tiles"); dir.create(tile_dir, showWarnings = FALSE)

  fit <- readRDS(fit_fp)
  mod <- readRDS(mod_fp)

  # reduce stack if we can detect the model's covariates
  stk_use <- subset_stack_for_model(stk, mod)

  mean_tiles <- character(0)
  sd_tiles   <- character(0)

  for (i in seq_along(tile_exts)) {
    te   <- tile_exts[[i]]
    tile <- crop(stk_use, te)
    if (ncell(tile) == 0) next

    nas_first <- global(is.na(tile[[1]]), "sum", na.rm = TRUE)[1,1]
    if (!is.na(nas_first) && nas_first == ncell(tile)) next

    preds <- safe_predict_tile(mod, fit, tile, draws)

    pm <- if (inherits(preds$mean, "SpatRaster")) preds$mean else rast(preds$mean)
    ps <- if (inherits(preds$sd,   "SpatRaster")) preds$sd   else rast(preds$sd)

    if (!identical(names(pm), names(ps))) {
      stop("Within-tile name mismatch (mean vs sd) at tile ", i)
    }
    if (any(!nzchar(names(pm)))) {
      stop("Empty band names found in tile ", i, "; prediction should carry species names.")
    }

    mf <- file.path(tile_dir, sprintf("psi_mean_tile_%03d.tif", i))
    sf <- file.path(tile_dir, sprintf("psi_sd_tile_%03d.tif",   i))
    writeRaster(pm, mf, overwrite=TRUE, NAflag=-9999,
                wopt=list(gdal=gdal_opt, datatype=datatype))
    writeRaster(ps, sf, overwrite=TRUE, NAflag=-9999,
                wopt=list(gdal=gdal_opt, datatype=datatype))

    mean_tiles <- c(mean_tiles, mf)
    sd_tiles   <- c(sd_tiles,   sf)

    message(sprintf("  tile %d/%d written", i, length(tile_exts)))
    rm(preds, pm, ps, tile); gc(FALSE)
  }

  if (!length(mean_tiles)) stop("No tiles were written (all tiles empty?).")

  list(mean_tiles = mean_tiles, sd_tiles = sd_tiles)
}

# 6) Mosaic & outputs 
finalize_mosaic_from_tiles <- function(mean_tiles, sd_tiles, out_dir,
                                       datatype = "FLT4S",
                                       gdal_opt = c("COMPRESS=DEFLATE","PREDICTOR=3","ZLEVEL=6")) {
  stopifnot(length(mean_tiles) == length(sd_tiles), length(mean_tiles) > 0)

  first_mean <- rast(mean_tiles[1])
  sp_names   <- names(first_mean)
  if (!length(sp_names) || any(!nzchar(sp_names))) {
    stop("Tile bands have missing names; cannot determine species.")
  }

  check_names <- function(files, ref_names, label) {
    for (f in files) {
      nm <- names(rast(f))
      if (!identical(nm, ref_names)) {
        stop("Band-name mismatch across ", label, " tiles.\n",
             "Reference names: ", paste(ref_names, collapse=", "), "\n",
             "Mismatched file: ", f, "\n",
             "Its names:       ", paste(nm, collapse=", "), "\n",
             "Fix tile generation before mosaicking.")
      }
    }
  }
  check_names(mean_tiles, sp_names, "mean")
  check_names(sd_tiles,   sp_names, "sd")

  sb_dir <- file.path(out_dir, "_tmp_singlebands")
  dir.create(sb_dir, showWarnings = FALSE)

  mean_species_fps <- character(length(sp_names))
  sd_species_fps   <- character(length(sp_names))

  for (k in seq_along(sp_names)) {
    nm <- sp_names[k]

    mean_sb_files <- character(0)
    sd_sb_files   <- character(0)
    for (t in seq_along(mean_tiles)) {
      r_m <- rast(mean_tiles[t])[[nm]]
      f_m <- file.path(sb_dir, sprintf("mean_%s_tile%03d.tif", nm, t))
      writeRaster(r_m, f_m, overwrite=TRUE, NAflag=-9999,
                  wopt=list(gdal=gdal_opt, datatype=datatype))
      mean_sb_files <- c(mean_sb_files, f_m)

      r_s <- rast(sd_tiles[t])[[nm]]
      f_s <- file.path(sb_dir, sprintf("sd_%s_tile%03d.tif", nm, t))
      writeRaster(r_s, f_s, overwrite=TRUE, NAflag=-9999,
                  wopt=list(gdal=gdal_opt, datatype=datatype))
      sd_sb_files <- c(sd_sb_files, f_s)
    }

    vrm <- if (length(mean_sb_files) == 1) rast(mean_sb_files) else vrt(mean_sb_files)
    vrs <- if (length(sd_sb_files)   == 1) rast(sd_sb_files)   else vrt(sd_sb_files)

    mean_species_fps[k] <- file.path(out_dir, paste0("psi_mean_", nm, ".tif"))
    sd_species_fps[k]   <- file.path(out_dir, paste0("psi_sd_",   nm, ".tif"))

    writeRaster(vrm, mean_species_fps[k], overwrite=TRUE, NAflag=-9999,
                wopt=list(gdal=gdal_opt, datatype=datatype))
    writeRaster(vrs, sd_species_fps[k],   overwrite=TRUE, NAflag=-9999,
                wopt=list(gdal=gdal_opt, datatype=datatype))
  }

  r_mean <- rast(mean_species_fps); names(r_mean) <- sp_names
  r_sd   <- rast(sd_species_fps);   names(r_sd)   <- sp_names

  mean_fp <- file.path(out_dir, "psi_mean.tif")
  sd_fp   <- file.path(out_dir, "psi_sd.tif")
  writeRaster(r_mean, mean_fp, overwrite=TRUE, NAflag=-9999,
              wopt=list(gdal=gdal_opt, datatype=datatype))
  writeRaster(r_sd,   sd_fp,   overwrite=TRUE, NAflag=-9999,
              wopt=list(gdal=gdal_opt, datatype=datatype))

  psi_comm_mean <- app(r_mean, mean, na.rm=TRUE); names(psi_comm_mean) <- "psi_community_mean"
  writeRaster(psi_comm_mean, file.path(out_dir, "psi_mean_COMMUNITY.tif"),
              overwrite=TRUE, NAflag=-9999, wopt=list(gdal=gdal_opt, datatype=datatype))

  pdf(file.path(out_dir, "psi_MEAN_per_species.pdf"), width=11, height=8.5)
  for (nm in names(r_mean)) plot(r_mean[[nm]], main=nm, col=hcl.colors(100), zlim=c(0,1))
  dev.off()

  pdf(file.path(out_dir, "psi_SD_per_species.pdf"), width=11, height=8.5)
  for (nm in names(r_sd)) plot(r_sd[[nm]], main=nm, col=hcl.colors(100), zlim=c(0,0.5))
  dev.off()

  pdf(file.path(out_dir, "psi_COMMUNITY_MEAN.pdf"), width=11, height=8.5)
  plot(psi_comm_mean, main="Community mean psi", col=hcl.colors(100), zlim=c(0,1))
  dev.off()

  invisible(list(mean_fp = mean_fp, sd_fp = sd_fp,
                 mean_species = mean_species_fps, sd_species = sd_species_fps))
}

# =================== RUN ALL MODEL PAIRS (GT: *_fit.rds ↔ *_mod.rds) ===================
tile_exts <- mk_tile_extents(stk, nx = nx, ny = ny)

fit_files <- list.files(models_dir, pattern = "_fit\\.rds$", full.names = TRUE, ignore.case = TRUE)
stopifnot(length(fit_files) > 0)

for (fit_fp in fit_files) {
  mod_fp <- file.path(models_dir, sub("_fit\\.rds$", "_mod.rds", basename(fit_fp), ignore.case = TRUE))
  if (!file.exists(mod_fp)) { warning("Missing mod for: ", basename(fit_fp)); next }

  out_dir <- file.path(out_root, sub("\\.rds$", "", basename(fit_fp), ignore.case = TRUE))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  tiles <- write_tiles_for_model(fit_fp, mod_fp, stk, tile_exts, out_dir,
                                 draws = draws, datatype = datatype, gdal_opt = gdal_opt)

  finalize_mosaic_from_tiles(tiles$mean_tiles, tiles$sd_tiles, out_dir,
                             datatype = datatype, gdal_opt = gdal_opt)

  # Optionally remove tiles to save space:
  # unlink(file.path(out_dir, "_tiles"), recursive = TRUE, force = TRUE)
}










################ clipping to study extent:



suppressPackageStartupMessages({
  library(terra)
  library(sf)
})

### Inputs 
proj_root <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL"
areas_shp <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/areas/Areas.shp/Areas.shp"
suffix    <- "_CLIPPED_TO_AREAS"  # outputs go under each model folder in this subdir
datatype  <- "FLT4S"              # float32 for psi
touches   <- FALSE                # TRUE to keep rim cells; FALSE = cell center inside polygon


# GeoTIFF options
gdal_opt_float <- c("COMPRESS=DEFLATE","PREDICTOR=3","ZLEVEL=6")
gdal_opt_mask  <- c("COMPRESS=DEFLATE","ZLEVEL=6")

stopifnot(dir.exists(proj_root), file.exists(areas_shp))

# Find model folders that actually have psi outputs 

all_mean <- list.files(proj_root, pattern="^psi_mean\\.tif$", full.names=TRUE,
                       recursive=TRUE, ignore.case=TRUE)
all_mean <- all_mean[!grepl("_CLIPPED_|_MASK", all_mean, ignore.case=TRUE)]
stopifnot(length(all_mean) > 0)
model_dirs <- unique(dirname(all_mean))
cat("Found", length(model_dirs), "model folder(s) with psi_mean.tif\n")

# Reference ψ grid (CRS/extent/resolution)
r_ref <- rast(all_mean[1])
cat("\n[ψ REF GRID]\n  Extent   : ",
    paste(signif(as.vector(ext(r_ref)), 6), collapse=", "),
    "\n  Resolution: ",
    paste(signif(res(r_ref), 8), collapse=" x "), "\n\n", sep="")

## helpers 
extent_overlaps <- function(E1, E2){
  e1 <- c(xmin(E1), xmax(E1), ymin(E1), ymax(E1))
  e2 <- c(xmin(E2), xmax(E2), ymin(E2), ymax(E2))
  if (any(!is.finite(c(e1, e2)))) return(FALSE)
  (e1[1] < e2[2]) && (e1[2] > e2[1]) && (e1[3] < e2[4]) && (e1[4] > e2[3])
}

## read Areas.shp with sf, standardize to WGS84

areas_sf <- suppressWarnings(st_read(areas_shp, quiet = TRUE))
if (is.na(st_crs(areas_sf))) {
  warning("Areas.shp has empty CRS; assigning WGS84 (EPSG:4326) as a fallback. Verify!")
  st_crs(areas_sf) <- 4326
}

areas_ll  <- st_transform(areas_sf, 4326)

# Convert to terra and project directly to the raster grid CRS
areas_ll_v  <- vect(areas_ll)
areas_ref_v <- project(areas_ll_v, r_ref)  # terra::project with raster target ensures exact CRS match

cat("[AREAS]\n  CRS match to ψ grid: ", compareCRS(areas_ref_v, r_ref), "\n", sep="")
cat("  Areas extent (ψ CRS): ",
    paste(signif(as.vector(ext(areas_ref_v)), 6), collapse=", "),
    "\n", sep="")

if (!extent_overlaps(ext(areas_ref_v), ext(r_ref))) {
  stop("Areas.shp polygon(s) and reference raster do not overlap after reprojection.\n",
       "Check that Areas.shp is the intended Guatemala study extent.")
}

# rasterize a mask on ψ grid 
mask_ref <- rast(r_ref[[1]]); mask_ref[] <- NA_real_
mask_ref <- rasterize(areas_ref_v, mask_ref, field=1, background=NA, touches=touches)
mask_ref <- trim(mask_ref)

# Write the mask (byte)
mask_dir <- file.path(proj_root, "_AREAS_MASK"); dir.create(mask_dir, showWarnings=FALSE, recursive=TRUE)
mask_fp  <- file.path(mask_dir, "AREAS_MASK_aligned_to_reference.tif")
writeRaster(mask_ref, mask_fp, overwrite=TRUE, NAflag=-9999,
            wopt=list(datatype="INT1U", gdal=gdal_opt_mask))
cat("Mask written:\n  ", mask_fp, "\n\n", sep="")

# clip
clip_model <- function(model_dir, mask_r,
                       datatype="FLT4S",
                       gdal_opt_float=c("COMPRESS=DEFLATE","PREDICTOR=3","ZLEVEL=6")) {

  mean_fp <- file.path(model_dir, "psi_mean.tif")
  sd_fp   <- file.path(model_dir, "psi_sd.tif")
  if (!file.exists(mean_fp) || !file.exists(sd_fp)) {
    message("Skipping (missing psi_mean/psi_sd): ", model_dir)
    return(invisible(FALSE))
  }

  r_mean <- rast(mean_fp)
  r_sd   <- rast(sd_fp)

  # Align to mask grid if needed
  if (!compareGeom(r_mean, mask_r, stopOnError = FALSE)) r_mean <- resample(r_mean, mask_r, method="bilinear")
  if (!compareGeom(r_sd,   mask_r, stopOnError = FALSE)) r_sd   <- resample(r_sd,   mask_r, method="bilinear")

  # Crop to mask extent and apply mask
  r_mean_c <- mask(crop(r_mean, ext(mask_r), snap="out"), mask_r)
  r_sd_c   <- mask(crop(r_sd,   ext(mask_r), snap="out"), mask_r)

  out_dir <- file.path(model_dir, suffix)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # Write multi-band clipped stacks
  writeRaster(r_mean_c, file.path(out_dir, "psi_mean_CLIPPED.tif"),
              overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype=datatype, gdal=gdal_opt_float))
  writeRaster(r_sd_c,   file.path(out_dir, "psi_sd_CLIPPED.tif"),
              overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype=datatype, gdal=gdal_opt_float))

  # Per-species single-band files
  for (nm in names(r_mean_c)) {
    writeRaster(r_mean_c[[nm]], file.path(out_dir, paste0("psi_mean_", nm, "_CLIPPED.tif")),
                overwrite=TRUE, NAflag=-9999,
                wopt=list(datatype=datatype, gdal=gdal_opt_float))
    writeRaster(r_sd_c[[nm]],   file.path(out_dir, paste0("psi_sd_",   nm, "_CLIPPED.tif")),
                overwrite=TRUE, NAflag=-9999,
                wopt=list(datatype=datatype, gdal=gdal_opt_float))
  }

  # Community mean + PDFs
  psi_comm_mean <- app(r_mean_c, mean, na.rm=TRUE); names(psi_comm_mean) <- "psi_community_mean"
  writeRaster(psi_comm_mean, file.path(out_dir, "psi_mean_COMMUNITY_CLIPPED.tif"),
              overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype=datatype, gdal=gdal_opt_float))

  pdf(file.path(out_dir, "psi_MEAN_per_species_CLIPPED.pdf"), width=11, height=8.5)
  for (nm in names(r_mean_c)) plot(r_mean_c[[nm]], main=paste0(nm," (psi mean)"),
                                   col=hcl.colors(100), zlim=c(0,1))
  dev.off()

  pdf(file.path(out_dir, "psi_SD_per_species_CLIPPED.pdf"), width=11, height=8.5)
  for (nm in names(r_sd_c)) plot(r_sd_c[[nm]], main=paste0(nm," (psi sd)"),
                                 col=hcl.colors(100), zlim=c(0,0.5))
  dev.off()

  pdf(file.path(out_dir, "psi_COMMUNITY_MEAN_CLIPPED.pdf"), width=11, height=8.5)
  plot(psi_comm_mean, main="Community mean psi", col=hcl.colors(100), zlim=c(0,1))
  dev.off()

  invisible(TRUE)
}

# run for every model 
for (md in model_dirs) {
  message("\nClipping: ", basename(md))
  ok <- try(clip_model(md, mask_ref, datatype=datatype, gdal_opt_float=gdal_opt_float), silent=TRUE)
  if (inherits(ok, "try-error")) {
    warning("Clip failed for ", basename(md), ": ", conditionMessage(ok))
  }
}














#################################################################### comparing AI vs expert reviewed psi projections

################# finally compare between expert and AI projections:


#
suppressPackageStartupMessages(library(terra))

# Inputs:
expert_example <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/GT_expertReview_ltd_DryOnly_d2hs_fit/_CLIPPED_TO_AREAS/psi_mean_Cuniculus paca_CLIPPED.tif"
ai_example     <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/GT_CleanedPairs50_ltd_DryOnly_d2hs_fit/_CLIPPED_TO_AREAS/psi_mean_Cuniculus paca_CLIPPED.tif"

stopifnot(file.exists(expert_example), file.exists(ai_example))

# Folders that contain all species-specific, clipped rasters
expert_dir <- dirname(expert_example)
ai_dir     <- dirname(ai_example)

# Put output at the same level as _HEX_MASK 
proj_root  <- normalizePath(file.path(expert_dir, "..", ".."), winslash = "/", mustWork = FALSE)
out_dir    <- file.path(proj_root, "_EXPERT_vs_AI")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# find species present in BOTH expert and ai 
get_species <- function(dir){
  fs <- list.files(dir, pattern = "^psi_mean_.*_CLIPPED\\.tif$", full.names = FALSE)
  # extract species name between "psi_mean_" and "_CLIPPED.tif"
  sub("^psi_mean_(.*)_CLIPPED\\.tif$", "\\1", fs)
}

sp_expert <- get_species(expert_dir)
sp_ai     <- get_species(ai_dir)
shared    <- intersect(sp_expert, sp_ai)

if (!length(shared)) stop("No shared species between Expert and AI folders.")

cat("Shared species (n=", length(shared), "):\n  ", paste(shared, collapse=", "), "\n", sep="")

# helpers
align_to <- function(src, tgt){
  # ensures src matches tgt grid (CRS/res/extent)
  if (!compareCRS(src, tgt)) src <- project(src, tgt)
  if (!compareGeom(src, tgt, stopOnError = FALSE)) src <- resample(src, tgt, method = "bilinear")
  src
}

# Common color ramps
pal_prob <- hcl.colors(100, "YlGnBu")
pal_div  <- hcl.colors(100, "Blue-Red 2")

# process each shared species
pdf_fp <- file.path(out_dir, "Expert_vs_AI_per_species.pdf")
pdf(pdf_fp, width = 8.5, height = 11)  # one species per page, 3 panels stacked

for (sp in shared) {
  # input files
  fp_exp <- file.path(expert_dir, paste0("psi_mean_", sp, "_CLIPPED.tif"))
  fp_ai  <- file.path(ai_dir,     paste0("psi_mean_", sp, "_CLIPPED.tif"))

  if (!file.exists(fp_exp) || !file.exists(fp_ai)) next

  r_exp <- rast(fp_exp)
  r_ai  <- rast(fp_ai)

  # Align AI to Expert grid
  r_ai_a <- align_to(r_ai, r_exp)

  # Outputs (copies + difference)
  out_exp <- file.path(out_dir, paste0("psi_mean_", sp, "_Expert.tif"))
  out_ai  <- file.path(out_dir, paste0("psi_mean_", sp, "_AI.tif"))
  out_dif <- file.path(out_dir, paste0("psi_mean_", sp, "_Difference.tif"))  # Expert - AI

  writeRaster(r_exp, out_exp, overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=DEFLATE")))
  writeRaster(r_ai_a, out_ai, overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=DEFLATE")))

  r_diff <- r_exp - r_ai_a
  writeRaster(r_diff, out_dif, overwrite=TRUE, NAflag=-9999,
              wopt=list(datatype="FLT4S", gdal=c("COMPRESS=DEFLATE")))

  # ---- page plot (Expert top; Difference middle; AI bottom) ----
  par(mfrow = c(3,1), mar = c(3,4,3,6), xpd = NA)

  plot(r_exp, col = pal_prob, zlim = c(0,1),
       main = paste0(sp, " — Expert (psi mean)"))
  box()

  # symmetric diff scale around 0 (cap to [-1, 1])
  zmax <- 1
  plot(clamp(r_diff, lower=-zmax, upper=zmax), col = pal_div, zlim = c(-zmax, zmax),
       main = paste0(sp, " — Difference (Expert − AI)"))
  box()

  plot(r_ai_a, col = pal_prob, zlim = c(0,1),
       main = paste0(sp, " — AI (psi mean)"))
  box()
}

dev.off()









############### comparison graphics improvement:

# --- Expert vs AI psi comparison (robust axes; fixed diff scale; consistent palettes) ---
suppressPackageStartupMessages(library(terra))

# EDIT: folder that already contains these per-species rasters:
#   psi_mean_<SPECIES>_Expert.tif
#   psi_mean_<SPECIES>_AI.tif
#   psi_mean_<SPECIES>_Difference.tif
cmp_dir <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/_EXPERT_vs_AI"
out_pdf <- file.path(cmp_dir, "Expert_vs_AI_comparison_V1.pdf")

stopifnot(dir.exists(cmp_dir))

# --- helpers ---
list_by_suffix <- function(sfx) {
  fs <- list.files(cmp_dir, pattern = paste0("^psi_mean_.*_", sfx, "\\.tif$"),
                   full.names = TRUE)
  if (!length(fs)) return(data.frame(species=character(), file=character()))
  spp <- sub(paste0("^psi_mean_(.*)_", sfx, "\\.tif$"), "\\1",
             basename(fs), perl = TRUE)
  data.frame(species = spp, file = fs, stringsAsFactors = FALSE)
}

# Axis ticks that are robust to ymax < ymin (and tiny extents)
draw_axes_int <- function(r){
  e <- ext(r)
  x0 <- min(e$xmin, e$xmax); x1 <- max(e$xmin, e$xmax)
  y0 <- min(e$ymin, e$ymax); y1 <- max(e$ymin, e$ymax)
  xt <- if ((x1 - x0) >= 1) seq(ceiling(x0), floor(x1), by = 1) else pretty(c(x0, x1))
  yt <- if ((y1 - y0) >= 1) seq(ceiling(y0), floor(y1), by = 1) else pretty(c(y0, y1))
  axis(1, at = xt, labels = xt, las = 1)
  axis(2, at = yt, labels = yt, las = 1)
  box()
}

# --- gather files ---
tab_exp <- list_by_suffix("Expert")
tab_ai  <- list_by_suffix("AI")
tab_dif <- list_by_suffix("Difference")

common <- Reduce(intersect, list(tab_exp$species, tab_ai$species, tab_dif$species))
stopifnot(length(common) > 0)

# order pages by species name
tab_exp <- tab_exp[match(sort(common), tab_exp$species), ]
tab_ai  <- tab_ai [match(tab_exp$species, tab_ai$species), ]
tab_dif <- tab_dif[match(tab_exp$species, tab_dif$species), ]

# --- palettes & fixed scales ---
pal_psi   <- hcl.colors(100, palette = "viridis")     # blue -> yellow
pal_diff  <- hcl.colors(100, palette = "Blue-Red 3")  # diverging
zlim_psi  <- c(0, 1)
zlim_diff <- c(-0.5, 0.5)  # FIXED shared range across ALL species, as requested

# --- make PDF ---
pdf(out_pdf, width = 10.5, height = 13)  # tall page (3 stacked panels)
op <- par(mfrow = c(3,1), mar = c(5,5,3,6), xpd = NA, cex.axis = 0.9, cex.main = 1.0)

for (i in seq_len(nrow(tab_exp))) {
  sp <- tab_exp$species[i]
  rE <- rast(tab_exp$file[i])
  rA <- rast(tab_ai$file[i])
  rD <- rast(tab_dif$file[i])

  # align (rarely needed, but safe)
  if (!compareGeom(rA, rE, stopOnError = FALSE)) rA <- resample(rA, rE, method = "bilinear")
  if (!compareGeom(rD, rE, stopOnError = FALSE)) rD <- resample(rD, rE, method = "bilinear")

  # 1) Expert (fixed 0–1)
  plot(rE, main = paste0(sp, " — Expert (psi)"),
       col = pal_psi, zlim = zlim_psi, axes = FALSE)
  draw_axes_int(rE)

  # 2) Difference (Expert − AI) with strictly fixed shared range [-0.5, +0.5]
  plot(clamp(rD, lower = zlim_diff[1], upper = zlim_diff[2]),
       main = paste0(sp, " — Difference (Expert − AI)"),
       col = pal_diff, zlim = zlim_diff, axes = FALSE)
  draw_axes_int(rD)

  # 3) AI (fixed 0–1)
  plot(rA, main = paste0(sp, " — AI (psi)"),
       col = pal_psi, zlim = zlim_psi, axes = FALSE)
  draw_axes_int(rA)
}

par(op); dev.off()
cat("Wrote comparison PDF:\n  ", out_pdf, "\n", sep = "")









################ differences pdf


# --- Difference-only (Expert − AI) PDF with fixed shared scale (-0.5 .. +0.5)
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(colorspace)   # for diverging_hcl palette
  library(scales)       # for oob = squish
  library(grid)         # for unit()
})

# Folder that already contains: psi_mean_<SPECIES>_Difference.tif
cmp_dir <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/_EXPERT_vs_AI"
out_pdf <- file.path(cmp_dir, "Expert_vs_AI_Difference_ONLY_fixed_0.5.pdf")

stopifnot(dir.exists(cmp_dir))

# List difference rasters
dif_fs  <- list.files(cmp_dir, pattern="^psi_mean_.*_Difference\\.tif$", full.names=TRUE)
stopifnot(length(dif_fs) > 0)
species <- sub("^psi_mean_(.*)_Difference\\.tif$", "\\1", basename(dif_fs), perl=TRUE)
ord     <- order(species); dif_fs <- dif_fs[ord]; species <- species[ord]

# Fixed legend/scale across ALL species
zlim <- c(-0.5, 0.5)
pal  <- diverging_hcl(11, palette = "Blue-Red 3") # symmetric diverging

# Safe integer ticks (handles flipped extents; falls back if span < 1)
int_ticks <- function(a, b) {
  lo <- min(a, b); hi <- max(a, b)
  if (is.finite(lo) && is.finite(hi) && (hi - lo) >= 1) {
    seq(ceiling(lo), floor(hi), by = 1)
  } else {
    pretty(c(lo, hi), n = 5)
  }
}

# Make a plot for one raster path
make_plot <- function(fp, title_text) {
  r <- rast(fp)
  # Clamp so the legend is strictly [-0.5, 0.5]
  r <- clamp(r, lower = zlim[1], upper = zlim[2], values = TRUE)

  # Tile sizing to avoid "uneven interval" warning
  rs <- res(r); dx <- abs(rs[1]); dy <- abs(rs[2])

  # Raster -> df
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "diff"

  # Whole-degree ticks (robust to flipped Y)
  ex     <- ext(r)
  xticks <- int_ticks(ex$xmin, ex$xmax)
  yticks <- int_ticks(ex$ymin, ex$ymax)

  ggplot(df, aes(x = x, y = y, fill = diff)) +
    geom_tile(width = dx, height = dy) +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = pal,
      limits  = zlim,
      oob     = squish,
      breaks  = seq(-0.5, 0.5, by = 0.1),
      labels  = number_format(accuracy = 0.1),
      name    = "Expert − AI"
    ) +
    scale_x_continuous(breaks = xticks) +
    scale_y_continuous(breaks = yticks) +
    guides(fill = guide_colorbar(barheight = unit(5, "cm"), barwidth = unit(0.6, "cm"))) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid  = element_blank(),
      axis.title  = element_blank(),
      plot.title  = element_text(face = "bold"),
      legend.title= element_text(),
      legend.text = element_text(),
      panel.border= element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
}

# Multi-page PDF (one species per page)
if (capabilities("cairo")) {
  grDevices::cairo_pdf(out_pdf, width = 11, height = 8.5)
} else {
  pdf(out_pdf, width = 11, height = 8.5)
}
for (i in seq_along(dif_fs)) {
  print(make_plot(dif_fs[i], paste0(species[i], " — Difference (Expert − AI)")))
}
dev.off()

cat("Wrote difference-only PDF with fixed scale [-0.5, 0.5]:\n  ", out_pdf, "\n", sep = "")






############## updated full code:


# --- Expert vs AI psi comparison (fixed palettes, axes, and scales; ggplot tiles) ---

suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(colorspace)   # for diverging_hcl
  library(scales)       # for oob = squish, label formats
  library(gridExtra)    # grid.arrange
  library(grid)
})

# ============ EDIT ============
cmp_dir <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/_EXPERT_vs_AI"
out_pdf <- file.path(cmp_dir, "Expert_vs_AI_comparison.pdf")
# ==============================

stopifnot(dir.exists(cmp_dir))

# --- list helpers ---
list_by_suffix <- function(sfx) {
  fs <- list.files(cmp_dir, pattern = paste0("^psi_mean_.*_", sfx, "\\.tif$"),
                   full.names = TRUE)
  if (!length(fs)) return(data.frame(species=character(), file=character()))
  spp <- sub(paste0("^psi_mean_(.*)_", sfx, "\\.tif$"), "\\1",
             basename(fs), perl = TRUE)
  data.frame(species = spp, file = fs, stringsAsFactors = FALSE)
}

tab_exp <- list_by_suffix("Expert")
tab_ai  <- list_by_suffix("AI")
tab_dif <- list_by_suffix("Difference")

# only species present in all three
common <- Reduce(intersect, list(tab_exp$species, tab_ai$species, tab_dif$species))
stopifnot(length(common) > 0)

tab_exp <- tab_exp[match(sort(common), tab_exp$species), ]
tab_ai  <- tab_ai [match(tab_exp$species, tab_ai$species), ]
tab_dif <- tab_dif[match(tab_exp$species, tab_dif$species), ]

# --- fixed scales & palettes ---
zlim_psi  <- c(0, 1)
zlim_diff <- c(-0.5, 0.5)   # shared across all species

pal_psi   <- hcl.colors(100, palette = "viridis")          # blue -> yellow
pal_diff  <- diverging_hcl(11, palette = "Blue-Red 3")     # symmetric diverging

# robust whole-degree tick helper (handles flipped extents; falls back if span < 1°)
int_ticks <- function(a, b) {
  lo <- min(a, b); hi <- max(a, b)
  if (is.finite(lo) && is.finite(hi) && (hi - lo) >= 1) {
    seq(ceiling(lo), floor(hi), by = 1)
  } else {
    pretty(c(lo, hi), n = 5)
  }
}

# --- helper to make a ggplot raster panel with whole-degree axes ---
make_panel <- function(fp, title_text, palette, zlim, is_diff = FALSE) {
  r <- rast(fp)

  # Clamp only for Difference so legend stays fixed and continuous
  if (is_diff) r <- clamp(r, lower = zlim[1], upper = zlim[2], values = TRUE)

  rs <- res(r); dx <- abs(rs[1]); dy <- abs(rs[2])   # avoid "uneven intervals" & negatives
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "val"

  ex <- ext(r)
  xticks <- int_ticks(ex$xmin, ex$xmax)
  yticks <- int_ticks(ex$ymin, ex$ymax)

  # legend breaks
  if (identical(zlim, zlim_psi)) {
    brks <- seq(0, 1, by = 0.1)
    labf <- number_format(accuracy = 0.1)
    leg_name <- "ψ"
  } else {
    brks <- seq(-0.5, 0.5, by = 0.1)
    labf <- number_format(accuracy = 0.1)
    leg_name <- "Expert − AI"
  }

  ggplot(df, aes(x = x, y = y, fill = val)) +
    geom_tile(width = dx, height = dy) +
    coord_fixed(expand = FALSE) +
    scale_fill_gradientn(
      colours = palette,
      limits  = zlim,
      oob     = squish,
      breaks  = brks,
      labels  = labf,
      name    = leg_name
    ) +
    scale_x_continuous(breaks = xticks) +
    scale_y_continuous(breaks = yticks) +
    guides(fill = guide_colorbar(barheight = unit(5, "cm"), barwidth = unit(0.6, "cm"))) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(),
      legend.text  = element_text(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
    )
}

# --- write a multi-page PDF (one species per page; Expert / Difference / AI) ---
if (capabilities("cairo")) grDevices::cairo_pdf(out_pdf, width = 10.5, height = 13) else pdf(out_pdf, width = 10.5, height = 13)

for (i in seq_len(nrow(tab_exp))) {
  sp   <- tab_exp$species[i]
  # align geometries via terra to guarantee identical grids
  rE <- rast(tab_exp$file[i])
  rA <- rast(tab_ai$file[i])
  rD <- rast(tab_dif$file[i])
  if (!compareGeom(rA, rE, stopOnError = FALSE)) rA <- resample(rA, rE, method = "bilinear")
  if (!compareGeom(rD, rE, stopOnError = FALSE)) rD <- resample(rD, rE, method = "bilinear")
  tfE <- tempfile(fileext = ".tif"); writeRaster(rE, tfE, overwrite = TRUE)
  tfA <- tempfile(fileext = ".tif"); writeRaster(rA, tfA, overwrite = TRUE)
  tfD <- tempfile(fileext = ".tif"); writeRaster(rD, tfD, overwrite = TRUE)

  pE <- make_panel(tfE, paste0(sp, " — Expert (ψ)"), pal_psi,  zlim_psi,  is_diff = FALSE)
  pD <- make_panel(tfD, paste0(sp, " — Difference"), pal_diff, zlim_diff, is_diff = TRUE)
  pA <- make_panel(tfA, paste0(sp, " — AI (ψ)"),    pal_psi,  zlim_psi,  is_diff = FALSE)

  grid.arrange(pE, pD, pA, ncol = 1, heights = c(1, 1, 1))
}

dev.off()
cat("Wrote comparison PDF:\n  ", out_pdf, "\n", sep = "")















##################### single png images for projections


# ==== High-res Expert & AI PNGs via ggplot — 2% buffer + centered legend ====
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(scales)
  library(colorspace)
})

# -------- EDIT (paths) --------
cmp_dir <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/_EXPERT_vs_AI"
out_dir <- file.path(cmp_dir, "_PNG_export_ggplot_03")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ---------------------------------------

list_by_suffix <- function(sfx){
  fs <- list.files(cmp_dir, pattern = paste0("^psi_mean_.*_", sfx, "\\.tif$"),
                   full.names = TRUE)
  if (!length(fs)) return(data.frame(species=character(), file=character()))
  spp <- sub(paste0("^psi_mean_(.*)_", sfx, "\\.tif$"), "\\1", basename(fs), perl=TRUE)
  data.frame(species=spp, file=fs, stringsAsFactors=FALSE)
}
deg_ticks <- function(lo, hi){
  lo2 <- min(lo, hi); hi2 <- max(lo, hi); span <- hi2 - lo2
  step <- if (is.finite(span) && span >= 2) 1 else 0.5
  seq(ceiling(lo2/step)*step, floor(hi2/step)*step, by=step)
}
align_to <- function(src, tgt){
  if (!compareCRS(src, tgt)) src <- project(src, tgt)
  if (!compareGeom(src, tgt, stopOnError = FALSE)) src <- resample(src, tgt, method="bilinear")
  src
}
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

# fixed palette/scale
pal_psi <- hcl.colors(256, "viridis")
zlim    <- c(0, 1)
brks    <- seq(0, 1, by = 0.25)
labf    <- number_format(accuracy = 0.25, trim = TRUE)

make_png_01 <- function(fp, title_text, out_png){
  r  <- rast(fp)[[1]]
  rs <- res(r); dx <- abs(rs[1]); dy <- abs(rs[2])
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "psi"

  ex     <- ext(r)
  xticks <- deg_ticks(ex$xmin, ex$xmax)
  yticks <- deg_ticks(ex$ymin, ex$ymax)

  p <- ggplot(df, aes(x = x, y = y, fill = psi)) +
    geom_raster(width = dx, height = dy, interpolate = FALSE) +
    # allow expansion & add ~2% buffer on both axes
    coord_fixed(clip = "off") +
    scale_x_continuous(breaks = xticks, expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = yticks, expand = expansion(mult = 0.02)) +
    scale_fill_gradientn(
      colours = pal_psi,
      limits  = zlim,
      oob     = squish,
      breaks  = brks,
      labels  = labf,
      name    = "ψ"
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      axis.title   = element_blank(),
      plot.title   = element_text(face = "bold"),
      legend.title = element_text(),
      legend.text  = element_text(),
      plot.margin  = margin(24, 56, 24, 28),     # a touch of breathing room
      # place legend slightly above center on the right edge
      legend.position      = c(1.02, 0.53),
      legend.justification = c("left","center")
    ) +
    guides(
      fill = guide_colorbar(
        title.position = "top",
        ticks = TRUE,
        ticks.colour = "black",
        barheight = unit(7, "cm"),
        barwidth  = unit(0.55, "cm"),
        direction = "vertical"
      )
    )

  ggsave(out_png, p, width = 7.5, height = 5.7, dpi = 600, units = "in", device = "png")
}

# gather & export
tab_exp <- list_by_suffix("Expert")
tab_ai  <- list_by_suffix("AI")
common  <- intersect(tab_exp$species, tab_ai$species)
stopifnot(length(common) > 0)

tab_exp <- tab_exp[match(sort(common), tab_exp$species), ]
tab_ai  <- tab_ai [match(tab_exp$species, tab_ai$species), ]

for (i in seq_len(nrow(tab_exp))) {
  sp   <- tab_exp$species[i]
  fexp <- tab_exp$file[i]
  fai  <- tab_ai$file[i]

  rE <- rast(fexp)
  rA <- align_to(rast(fai), rE)
  tfA <- tempfile(fileext = ".tif"); writeRaster(rA, tfA, overwrite = TRUE)

  out_exp <- file.path(out_dir, paste0("psi_mean_", safe_name(sp), "_Expert.png"))
  out_ai  <- file.path(out_dir, paste0("psi_mean_", safe_name(sp), "_AI.png"))

  make_png_01(fexp, paste0(sp, " — Expert (ψ)"), out_exp)
  make_png_01(tfA,  paste0(sp, " — AI (ψ)"),     out_ai)
}

cat("PNGs written to:\n  ", out_dir, "\n", sep = "")






############# single difference PNG files




# ==== Difference (Expert − AI) PNG export with fixed [-0.5, 0.5] scale ====
suppressPackageStartupMessages({
  library(terra)
  library(ggplot2)
  library(scales)
  library(colorspace)
  library(grid)  # for unit(), margin()
})

# -------- EDIT (paths) --------
cmp_dir <- "C:/Users/Home/Documents/1_AI_MSOMs/Covariate_SpatialLayers_Rasters/_GT_projections_TILES_FINAL/_EXPERT_vs_AI"
out_dir <- file.path(cmp_dir, "_PNG_export_Difference")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
# ---------------------------------------

# list per-species difference rasters created earlier:
dif_tab <- {
  fs  <- list.files(cmp_dir, pattern = "^psi_mean_.*_Difference\\.tif$", full.names = TRUE)
  stopifnot(length(fs) > 0)
  spp <- sub("^psi_mean_(.*)_Difference\\.tif$", "\\1", basename(fs), perl = TRUE)
  data.frame(species = spp, file = fs, stringsAsFactors = FALSE)
}

# helpers
deg_ticks <- function(lo, hi) {
  lo2 <- min(lo, hi); hi2 <- max(lo, hi); span <- hi2 - lo2
  step <- if (is.finite(span) && span >= 2) 1 else 0.5
  seq(ceiling(lo2/step)*step, floor(hi2/step)*step, by = step)
}
safe_name <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

# fixed scale + palette
zlim  <- c(-0.5, 0.5)
brks  <- seq(-0.5, 0.5, by = 0.25)
labf  <- number_format(accuracy = 0.25, trim = TRUE)
pal_d <- diverging_hcl(201, palette = "Blue-Red 3")  # smooth version of your PDF palette

make_png_diff <- function(fp, title_text, out_png) {
  r  <- rast(fp)[[1]]
  r  <- clamp(r, lower = zlim[1], upper = zlim[2], values = TRUE)

  rs <- res(r); dx <- abs(rs[1]); dy <- abs(rs[2])
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "diff"

  ex     <- ext(r)
  xticks <- deg_ticks(ex$xmin, ex$xmax)
  yticks <- deg_ticks(ex$ymin, ex$ymax)

  p <- ggplot(df, aes(x = x, y = y, fill = diff)) +
    geom_raster(width = dx, height = dy, interpolate = FALSE) +
    # 2% breathing room around the map, keep aspect ratio and border box
    coord_fixed(clip = "off") +
    scale_x_continuous(breaks = xticks, expand = expansion(mult = 0.02)) +
    scale_y_continuous(breaks = yticks, expand = expansion(mult = 0.02)) +
    scale_fill_gradientn(
      colours = pal_d,
      limits  = zlim,
      oob     = squish,
      breaks  = brks,
      labels  = labf,
      name    = NULL                 # <-- drop legend title
    ) +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
      axis.title   = element_blank(),
      plot.title   = element_text(face = "bold"),
      legend.title = element_blank(),                           # <-- ensure no title
      legend.text  = element_text(),
      # increase right margin so legend labels never clip
      plot.margin  = margin(24, 76, 24, 28),                    # <-- was 56; now wider
      legend.position      = c(1.02, 0.53),                     # just outside, centered
      legend.justification = c("left","center")
    ) +
    guides(
      fill = guide_colorbar(
        ticks = TRUE, ticks.colour = "black",
        barheight = unit(7, "cm"),
        barwidth  = unit(0.55, "cm"),
        direction = "vertical"
      )
    )

  ggsave(out_png, p, width = 7.5, height = 5.7, dpi = 600, units = "in", device = "png")
}

# export one PNG per species
ord <- order(dif_tab$species)
dif_tab <- dif_tab[ord, ]

for (i in seq_len(nrow(dif_tab))) {
  sp <- dif_tab$species[i]
  fp <- dif_tab$file[i]
  out_png <- file.path(out_dir, paste0("psi_mean_", safe_name(sp), "_Difference.png"))
  make_png_diff(fp, paste0(sp, " — Difference (Expert − AI)"), out_png)
}

cat("Difference PNGs written to:\n  ", out_dir, "\n", sep = "")


