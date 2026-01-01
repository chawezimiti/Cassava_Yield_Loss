#Run the Africa code and remove data  rm(data)
#Then Run the nigeria soil code and remove data 
#This gets you the model functions


#================== PREDICTIONS ================================================
library(readr)
library(dplyr)
library(BLA)
library(aplpack)
library(MASS)
library(ggplot2)
library(terra) # mapping
library(sf) # mapping
source("Extra_functions.R")

# 1) Read in data_pred ==============================================================

data_pred <- read_csv("all_merged_nigeria_pred_grid.csv")
names(data_pred)

n <- dim(data_pred)[1]
#Growing_year <- rep(c(2020,2021,2022), each=n)
Growing_year <- rep(c(2020), each=n)
data_pred <- rbind(data_pred, data_pred, data_pred)
data_pred <- cbind(data_pred,Growing_year )



# 2) data_pred manipulation and Feature engineering ================================================================================

length(which(is.na(data_pred$Growing_year)))# check for NA in the experiment year. To be used to aggregate climatic variables at site.
unique(data_pred$Growing_year)

data_pred <- data_pred %>% mutate(Growing_year = as.integer(Growing_year)) #converts the experimental year column to integer
years <- sort(unique(na.omit(data_pred$Growing_year)))

# (check that planting_month is the same)
data_pred$planting_month <- data_pred$SowMonth_Alice  # change to the appropriate column that contains sowing dates
data_pred$planting_month <- as.character(ceiling(data_pred$planting_month))

length(which(is.na(data_pred$planting_month)))# check for NA in the Planting month. PM is used as start point for aggregation of climatic variables
unique(data_pred$planting_month) # check if the values for PM are sensible between 1 and 12

#data_pred$planting_month <- ifelse(data_pred$planting_month==0,1,data_pred$planting_month)
#unique(data_pred$planting_month) # check if the values for PM are sensible between 1 and 12

## a) Cumulative Precipitation =================================================

# The Cumulative rainfall is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), veg stage (2-4 months),
# bulk (5-8 months) and mat (9-12 months).

rain_cols <- grep("^(?:pmm|ppm)_\\d{4}_\\d{1,2}$",
                  names(data_pred),
                  value = TRUE,
                  ignore.case = TRUE)

if (length(rain_cols) == 0L) {
  data_pred$pmm_1_2  <- NA_real_
  data_pred$pmm_3_4  <- NA_real_
  data_pred$pmm_5_8  <- NA_real_
  data_pred$pmm_9_12 <- NA_real_
} else {
  ## 2) Parse year + month from column names (e.g. pmm_1990_3)
  rx <- "^(?:pmm|ppm)_(\\d{4})_(\\d{1,2})$"
  mm_yy <- do.call(
    rbind,
    regmatches(rain_cols, regexec(rx, rain_cols, ignore.case = TRUE))
  )
  # mm_yy columns: full match, year, month
  yr_int  <- as.integer(mm_yy[, 2])
  mon_int <- as.integer(mm_yy[, 3])
  
  # Pretend each is the 1st of the month: YYYY-MM-01
  rain_dates <- as.Date(sprintf("%04d-%02d-01", yr_int, mon_int))
  
  # Order both dates and column names chronologically
  ord <- order(rain_dates)
  rain_dates <- rain_dates[ord]
  rain_cols  <- rain_cols[ord]
  
  ## 3) Convert rainfall values into a numeric matrix
  pmm_mat <- as.matrix(data_pred[rain_cols])
  storage.mode(pmm_mat) <- "double"
  
  ## 4) Precompute 12-month windows by (Growing_year, planting_month)
  data_pred$Growing_year <- as.integer(
    ifelse(nchar(data_pred$Growing_year) > 4,
           substr(data_pred$Growing_year, 1, 4),
           data_pred$Growing_year)
  )
  
  years  <- sort(unique(data_pred$Growing_year))
  months <- sprintf("%02d", 1:12)
  
  idx_by_yrmon <- vector("list", length = length(years) * 12L)
  names(idx_by_yrmon) <- as.vector(outer(years, months, paste, sep = "_"))
  
  for (yr in years) {
    for (m in months) {
      start <- as.Date(sprintf("%04d-%s-01", yr, m))
      end   <- as.Date(sprintf("%04d-%s-01", yr + 1L, m))  # exclusive
      idx   <- which(rain_dates >= start & rain_dates < end)
      idx_by_yrmon[[paste0(yr, "_", m)]] <- idx
    }
  }
  
  ## 5) Compute rainfall for 1–2, 3–4, 5–8, 9–12 months after planting (check that planting_month is the same)
  keys <- paste0(
    data_pred$Growing_year, "_",
    sprintf("%02d", as.integer(data_pred$planting_month))
  )
  idx_list <- idx_by_yrmon[keys]
  
  res_mat <- t(mapply(
    function(i, idx) {
      # If no data at all for the 12-month window
      if (length(idx) == 0L) return(rep(NA_real_, 4))
      
      dates_i <- rain_dates[idx]
      vals_i  <- pmm_mat[i, idx]
      
      if (all(is.na(vals_i))) return(rep(NA_real_, 4))
      
      # Planting date for this row (month 1)
      yr <- data_pred$Growing_year[i]
      pm <- as.integer(data_pred$planting_month[i])
      plant_date <- as.Date(sprintf("%04d-%02d-01", yr, pm))
      
      # Month offset from planting month (0 = planting month)
      off_year  <- as.integer(format(dates_i, "%Y")) - yr
      off_month <- as.integer(format(dates_i, "%m")) - pm
      month_off <- off_year * 12 + off_month  # ideally 0..11
      
      sum_range <- function(lo, hi) {
        sel <- which(month_off >= lo & month_off <= hi)
        if (length(sel) == 0L) return(NA_real_)
        v <- vals_i[sel]
        if (all(is.na(v))) return(NA_real_)
        sum(v, na.rm = TRUE)
      }
      
      p1_2  <- sum_range(0, 0)   # months 1
      p3_4  <- sum_range(1, 3)   # months 2–4
      p5_8  <- sum_range(4, 8)   # months 5–8
      p9_12 <- sum_range(9, 11)  # months 9–12
      
      c(p1_2, p3_4, p5_8, p9_12)
    },
    seq_len(nrow(data_pred)), idx_list,
    SIMPLIFY = TRUE
  ))
  
  data_pred$pmm_est  <- res_mat[, 1]
  data_pred$pmm_veg  <- res_mat[, 2]
  data_pred$pmm_bulk  <- res_mat[, 3]
  data_pred$pmm_mat <- res_mat[, 4]
}

# First exploration of plots of yield~rainfall by regions of Africa (East, West, South, central)


## b) Temperature Indices =============================================================

# Here we create Growing degree days, Heat stress degree days and chilling degree days.
# The GDD, HSDD and CDD is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), veg stage (2-4 months),
# bulk (5-8 months) and mat (9-12 months).
localIntPM=as.integer(data_pred$planting_month)

data_pred$planting_date <- as.Date(
  sprintf("%04d-%02d-01", data_pred$Growing_year, localIntPM)
)

## Identify daily Tmin/Tmax columns and common dates ---------------------------

tmin_cols  <- grep("^tmin_\\d{4}_\\d{2}_\\d{2}$", names(data_pred), value = TRUE)
tmax_cols  <- grep("^tmax_\\d{4}_\\d{2}_\\d{2}$", names(data_pred), value = TRUE)

tmin_dates <- as.Date(sub("^tmin_", "", tmin_cols), format = "%Y_%m_%d")
tmax_dates <- as.Date(sub("^tmax_", "", tmax_cols), format = "%Y_%m_%d")

common_dates <- sort(as.Date(intersect(tmin_dates, tmax_dates)))

if (length(common_dates) == 0L) {
  # no climatic overlap → everything NA
  data_pred$GDD_est          <- NA_real_
  data_pred$GDD_veg     <- NA_real_
  data_pred$GDD_bulk        <- NA_real_
  data_pred$GDD_mat     <- NA_real_
  
  data_pred$HSDD_est         <- NA_real_
  data_pred$HSDD_veg    <- NA_real_
  data_pred$HSDD_bulk       <- NA_real_
  data_pred$HSDD_mat    <- NA_real_
  
  data_pred$CDD_est          <- NA_real_
  data_pred$CDD_veg     <- NA_real_
  data_pred$CDD_bulk        <- NA_real_
  data_pred$CDD_mat     <- NA_real_
  
} else {
  
  ## Align daily temperature matrices ----------------------------------
  
  tmin_common_cols <- paste0("tmin_", format(common_dates, "%Y_%m_%d"))
  tmax_common_cols <- paste0("tmax_", format(common_dates, "%Y_%m_%d"))
  
  tmin_mat <- as.matrix(data_pred[tmin_common_cols]); storage.mode(tmin_mat) <- "double"
  tmax_mat <- as.matrix(data_pred[tmax_common_cols]); storage.mode(tmax_mat) <- "double"
  tmean_mat <- (tmin_mat + tmax_mat) / 2
  
  ## Thresholds + stages (cassava) -------------------------------------
  
  Tbase_gdd   <- 10  # GDD base
  Tcrit_heat  <- 32  # heat-stress threshold 
  Tbase_chill <- 18 # cold-stress threshold 
  
  # stages in days after planting (DAP); upper bound is exclusive
  stages <- list(
    est       = c(0,   30),   # 0–60 DAP
    veg  = c(30,  120),  # 60–120 DAP
    bulk     = c(120, 240),  # 120–240 DAP
    mat  = c(240, 365)   # 240–365 DAP
  )
  
  
  ## Row-wise calculation of GDD, HSDD, CDD by stage -------------------
  
  res <- mapply(
    function(i) {
      plant_date <- data_pred$planting_date[i]
      
      if (is.na(plant_date)) {
        return(rep(NA_real_, 12L))
      }
      
      out <- numeric(12L)
      names(out) <- c(
        "GDD_est", "GDD_veg", "GDD_bulk", "GDD_mat",
        "HSDD_est","HSDD_veg","HSDD_bulk","HSDD_mat",
        "CDD_est", "CDD_veg", "CDD_bulk", "CDD_mat"
      )
      
      k <- 1L
      
      for (st in names(stages)) {
        offset <- stages[[st]]
        start  <- plant_date + offset[1]
        end    <- plant_date + offset[2]
        
        idx <- which(common_dates >= start & common_dates < end)
        
        if (length(idx) == 0L) {
          out[k + 0L] <- NA_real_ # GDD
          out[k + 4L] <- NA_real_ # HSDD
          out[k + 8L] <- NA_real_ # CDD
        } else {
          tmean <- tmean_mat[i, idx]
          tmax  <- tmax_mat[i, idx]
          
          # 1) GDD
          gdd  <- sum(pmax(0, tmean - Tbase_gdd),   na.rm = TRUE)
          
          # 2) HSDD (using Tmax; swap to tmean if you prefer)
          hsdd <- sum(pmax(0, tmax  - Tcrit_heat),  na.rm = TRUE)
          
          # 3) CDD (cold stress)
          cdd  <- sum(pmax(0, Tbase_chill - tmean), na.rm = TRUE)
          
          out[k + 0L] <- gdd
          out[k + 4L] <- hsdd
          out[k + 8L] <- cdd
        }
        
        k <- k + 1L
      }
      
      out
    },
    seq_len(nrow(data_pred))
  )
  
  res <- t(res)
  
  ## Attach back to data_pred -----------------------------------------------
  
  data_pred$GDD_est   <- res[, "GDD_est"]
  data_pred$GDD_veg   <- res[, "GDD_veg"]
  data_pred$GDD_bulk  <- res[, "GDD_bulk"]
  data_pred$GDD_mat   <- res[, "GDD_mat"]
  data_pred$GDD_total <- rowSums(data_pred[, c("GDD_est","GDD_veg","GDD_bulk","GDD_mat")],na.rm = TRUE)
  
  data_pred$HSDD_est   <- res[, "HSDD_est"]
  data_pred$HSDD_veg   <- res[, "HSDD_veg"]
  data_pred$HSDD_bulk  <- res[, "HSDD_bulk"]
  data_pred$HSDD_mat   <- res[, "HSDD_mat"]
  data_pred$HSDD_total <- rowSums(data_pred[, c("HSDD_est","HSDD_veg","HSDD_bulk","HSDD_mat")], na.rm = TRUE)
  
  data_pred$CDD_est    <- res[, "CDD_est"]
  data_pred$CDD_veg    <- res[, "CDD_veg"]
  data_pred$CDD_bulk   <- res[, "CDD_bulk"]
  data_pred$CDD_mat    <- res[, "CDD_mat"]
  data_pred$CDD_total  <- rowSums(data_pred[, c("CDD_est","CDD_veg","CDD_bulk","CDD_mat")], na.rm = TRUE)
}

# First exploration of plots of yield~rainfall by regions of Africa (East, West, South, central)


## c) standard precipitation-evaporation index===========================================

# Determine the average spei index over the growing period of 12 months
# The spei is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), veg stage (2-4 months),
# bulk (5-8 months) and mat (9-12 months).

spei_cols  <- grep("^spei_\\d{4}.\\d{2}.\\d{1,2}$", names(data_pred), value = TRUE)## Identify monthly SPEI columns with dates like spei_2007-10-2

if (length(spei_cols) == 0L) {
  data_pred$spei_est  <- NA_real_
  data_pred$spei_veg  <- NA_real_
  data_pred$spei_bulk  <- NA_real_
  data_pred$spei_mat <- NA_real_
  
} else {
  
  ## Convert to Date
  spei_dates <- as.Date(sub("^spei_", "", spei_cols), format = "%Y.%m.%d")
  
  ## Order by date and align columns
  ord <- order(spei_dates)
  spei_dates <- spei_dates[ord]
  spei_cols  <- spei_cols[ord]
  
  ## SPEI matrix
  spei_mat <- as.matrix(data_pred[spei_cols])
  storage.mode(spei_mat) <- "double"
  
  ## Clean experiment year (handle "2013-2016")
  data_pred$Growing_year <- ifelse(
    nchar(data_pred$Growing_year) > 4,
    substr(data_pred$Growing_year, 1, 4),
    data_pred$Growing_year
  )
  data_pred$Growing_year <- as.integer(data_pred$Growing_year)
  
  years  <- sort(unique(data_pred$Growing_year))
  months <- sprintf("%02d", 1:12)
  
  ## Build full year-month index
  idx_by_yrmon <- vector("list", length = length(years) * 12L)
  names(idx_by_yrmon) <- as.vector(outer(years, months, paste, sep = "_"))
  
  for (yr in years) {
    for (m in months) {
      start <- as.Date(sprintf("%04d-%s-01", yr, m))
      end   <- as.Date(sprintf("%04d-%s-01", yr + 1L, m))  # exclusive
      idx   <- which(spei_dates >= start & spei_dates < end)
      idx_by_yrmon[[paste0(yr, "_", m)]] <- idx
    }
  }
  
  ## Helper to compute period-averages
  get_period_mean <- function(i, year, p_start, p_end) {
    this_year <- as.integer(year)
    key_months <- sprintf("%02d", p_start:p_end)
    
    idx_vec <- unlist(idx_by_yrmon[paste0(this_year, "_", key_months)])
    if (length(idx_vec) == 0L) return(NA_real_)
    vals <- spei_mat[i, idx_vec]
    if (all(is.na(vals))) return(NA_real_)
    mean(vals, na.rm = TRUE)
  }
  
  ## Apply to each row
  data_pred$spei_est <- mapply(function(i, yr)
    get_period_mean(i, yr, 0, 1),
    seq_len(nrow(data_pred)), data_pred$Growing_year)
  
  data_pred$spei_veg <- mapply(function(i, yr)
    get_period_mean(i, yr, 2, 4),
    seq_len(nrow(data_pred)), data_pred$Growing_year)
  
  data_pred$spei_bulk <- mapply(function(i, yr)
    get_period_mean(i, yr, 5, 8),
    seq_len(nrow(data_pred)), data_pred$Growing_year)
  
  data_pred$spei_mat <- mapply(function(i, yr)
    get_period_mean(i, yr, 9, 12),
    seq_len(nrow(data_pred)), data_pred$Growing_year)
}




### Prediction==================================================================

#a) Latitude prediction

radiation <- predictBL(model_lat, data_pred$lat)

Africa <- st_read("ng.shp") %>% st_transform(4326)

# Convert your data (lon, lat, yield) into sf points
pts <- st_as_sf(data_pred, coords = c("long", "lat"), crs = 4326)


radiation = pmax(radiation,4 ) # using 4 as a min value from boudary fit - else extrapolation 

radiation2=radiation^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = radiation2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

# b)  Make prediction using filled in values
# checking each and bounding where extrapolation occurs Alice


rain_early <- predictBL(model_pmm_est, log(data_pred$pmm_est))

rainE_min=predictBL(model_pmm_est, 6.125) # this x value is the largest the model was fitted to - see figure in report

rain_early = pmax(rain_early,rainE_min)
#rain_early[rain_early <rainE_min] = NA

rain_early2=rain_early^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = rain_early2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )
rain_veg <- predictBL(model_pmm_veg, data_pred$pmm_veg)

rainV_min=6 # this x value is the largest the model was fitted to - see figure in report

rain_veg = pmax(rain_veg,rainV_min)
#rain_veg[rain_veg <rainV_min] = NA

rain_veg2=rain_veg^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = rain_veg2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


rain_bulking <- predictBL(model_pmm_bulk, log(data_pred$pmm_bulk))

rainB_min=predictBL(model_pmm_bulk, 3.125) # this x value is the largest the model was fitted to - see figure in report

rain_bulking = pmax(rain_bulking,rainB_min)
#rain_bulking[rain_bulking <rainB_min] = NA

rain_bulking2=rain_bulking^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = rain_bulking2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


rain_maturation <- predictBL(model_pmm_mat, sqrt(data_pred$pmm_mat))

rainM_min=predictBL(model_pmm_mat, 5.5) # this x value is the largest the model was fitted to - see figure in report

rain_maturation = pmax(rain_maturation,rainM_min)
#rain_maturation[rain_maturation <rainM_min] = NA

rain_maturation2=rain_maturation^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = rain_maturation2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


spei_early <- predictBL(model_spei_est, data_pred$spei_est)

speiE_min=predictBL(model_spei_est, 0.85) # this x value is the largest the model was fitted to - see figure in report
spei_early = pmax(spei_early,speiE_min)
#spei_early[spei_early <speiE_min] = NA
spei_early2=spei_early^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = spei_early2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )






spei_veg <- predictBL(model_spei_veg, data_pred$spei_veg)
speiV_min=predictBL(model_spei_veg, 0.85) # this x value is the largest the model was fitted to - see figure in report
spei_veg = pmax(spei_veg,speiV_min)
#spei_veg[spei_veg <speiV_min] = NA
spei_veg2=spei_veg^2

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = spei_veg2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


spei_bulking <- predictBL(model_spei_bulk, data_pred$spei_bulk)

speiB_min=predictBL(model_spei_bulk, -0.85) # this x value is the largest the model was fitted to - see figure in report
spei_bulking = pmax(spei_bulking,speiB_min)
#spei_bulking[spei_bulking <speiB_min] = NA
spei_bulking2=spei_bulking^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = spei_bulking2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

spei_maturation <- predictBL(model_spei_mat, data_pred$spei_mat)

speiM_min=predictBL(model_spei_mat, -0.8) # this x value is the largest the model was fitted to - see figure in report
spei_maturation = pmax(spei_maturation,speiM_min)
#spei_maturation[spei_maturation <speiM_min] = NA
spei_maturation2=spei_maturation^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = spei_maturation2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

#################################################################################################


r <- rast("yield_from_clay_Nigeria.tif")

crs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

Clay_pV <- extract(r, vect(pts_sf))
Clay_p=Clay_pV$yield_clay

#Clay_p[which(Clay_p==0)]=NA

#Clay_min=predictBL(model_clay,2.7)

Clay_min=5.4113

Clay_p=pmax(Clay_p,Clay_min^2)
Clay_p=sqrt(Clay_p)
Clay_p2=Clay_p^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = Clay_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

r <- rast("yield_from_Sand_Nigeria.tif")
crs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

Sand_pV <- extract(r, vect(pts_sf))
Sand_p=Sand_pV$yield_sand

Sand_p[which(Sand_p==0)]=NA

#Sand_min=predictBL(model_sand,3.6)

Sand_min=3.996258

Sand_p=sqrt(pmax(Sand_p,Sand_min^2))

Sand_p2=Sand_p^2




ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = Sand_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

r <- rast("yield_from_pH_Nigeria.tif")
crs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

pH_pV <- extract(r, vect(pts_sf))
pH_p=pH_pV$yield_pH

pH_p[which(pH_p==0)]=NA

#pH_min=predictBL(model_pH,5)
pH_min=5.349788

pH_p=sqrt(pmax(pH_p,pH_min^2))

pH_p2=pH_p^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = pH_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


r <- rast("yield_from_soc_Nigeria.tif")
crs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

SOC_pV <- extract(r, vect(pts_sf))
SOC_p=SOC_pV$yield_soc

SOC_p[which(SOC_p==0)]=NA

#SOC_min=predictBL(model_SOC,3.65)

SOC_min=6.844721

SOC_p=sqrt(pmax(SOC_p,SOC_min^2))

SOC_p2=SOC_p^2



ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = SOC_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


r <- rast("yield_from_bd_Nigeria.tif")
#plcrs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

BD_pV <- extract(r, vect(pts_sf))
BD_p=BD_pV$yield_bd

#BD_p[which(BD_p==0)]=NA

#BD_min=predictBL(model_BD,1.45)
#BD_min=predictBL(model_BD,1.53)
BD_min=5.186661
BD_p=sqrt(pmax(BD_p,BD_min^2))

BD_p2=BD_p^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = BD_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )

r <- rast("yield_from_aw_Nigeria.tif")
#plcrs(r)
pts_df <- data.frame(
  lon =data_pred$long,
  lat = data_pred$lat
)
pts_sf <- st_as_sf(
  pts_df,
  coords = c("lon", "lat"),
  crs = 4326   # WGS84 lon/lat
)

PAW_pV <- extract(r, vect(pts_sf))
PAW_p=PAW_pV$yield_aw

#PAW_p[which(PAW_p==0)]=NA

PAW_min=predictBL(model_paw,0.205)
PAW_min=3.579946

PAW_p=sqrt(pmax(PAW_p,PAW_min^2))

PAW_p2=PAW_p^2


ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = PAW_p2), size = 1.75, alpha = 0.9) +
  
  scale_color_gradient(low = "red",high = "green",name = "Attainable Yield") +# Continuous color gradient
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(title = "Attainable Yield Map", x = "Longitude", y = "Latitude"
  )


dat_out <- min_by_indexNA(radiation = radiation,
                        rain_early = rain_early,
                        #rain_veg = rain_veg,
                        rain_bulking = rain_bulking,
                        rain_maturation = rain_maturation,
                        spei_early=spei_early,
                        spei_veg=spei_veg,
                        spei_bulking=spei_bulking,
                        spei_maturation=spei_maturation,
                       # Clay_p,
                        Sand_p,
                        pH_p,
                        SOC_p,
                        BD_p,
                        PAW_p)


#summary(dat_spei_veg_clean)

#We know pH has the smallest plateau s there is a chance when things are allocaled to Ph they should be 
#non bounded

dat_out <- dat_out %>%
  dplyr::mutate(
    vector_name = if_else(
      vector_name == "pH_p" & min_value >= 7.93967,
      "Unbounded",
      vector_name
    )
  )


#par(mar = c(10, 4, 4, 2))  # bottom, left, top, right
#barplot( table(dat_out$vector_name),las = 2, ylab = "Count",xlab = "",main = "Frequency of Names")
par(mar = c(9, 4, 4, 2))  # bottom, left, top, right
barplot(
  table(dat_out$vector_name),
  las = 2,
  ylab = "Count",
  xlab = "",
  main = "Frequency of Names",
  col = "skyblue"
)

data_pred <- cbind(data_pred,dat_out)
data_pred$min_value=data_pred$min_value^2

#plot(data_pred$spei_bulk,data_pred$min_value)
#length(data_pred$spei_bulk)
#length(data_pred$min_value)

# c.v) mapping the rainfall limited yield
Africa <- st_read("ng.shp") %>% st_transform(4326)
data_pred <- data_pred %>% filter(!is.na(long),!is.na(lat),!is.na(min_value))
pts <- st_as_sf(data_pred, coords = c("long", "lat"), crs = 4326)



# limiting yield-----------------

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = min_value), size = 2, alpha = 0.9) +
  scale_color_viridis_c(
#  scale_color_gradient(
#    low = "red",
#    high = "green",
    name = expression("Yield / t ha"^{-1})
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Predicted Yield by the Law of the Minimum: Soil and climate factors",
    x = "Longitude",
    y = "Latitude"
  )


mean(dat_out$min_value, na.rm=TRUE)^2

# limiting factor--------------------

# 1️⃣ Define all categories and factor levels explicitly
all_levels <- c(
  "radiation", "rain_early", "rain_veg", "rain_bulking", "rain_maturation",
  "spei_early", "spei_veg", "spei_bulking", "spei_maturation",
  "Clay_p", "SOC_p", "Sand_p2", "BD_p", "PAW_p", "pH_p", "Unbounded"
)

pts <- pts %>%
  mutate(vector_name = factor(vector_name, levels = all_levels))


my_colours <- c(
  # Radiation & Rain (blues)
  "radiation" = "#1f78b4",      # strong blue
  "rain_early" = "#6baed6",     # medium blue
  "rain_veg" = "#9ecae1",       # light blue
  "rain_bulking" = "#c6dbef",   # pale blue
  "rain_maturation" = "#084594",# dark blue
  
  # SPEI (greens)
  "spei_early" = "#31a354",     # strong green
  "spei_veg" = "#74c476",       # medium green
  "spei_bulking" = "#a1d99b",   # light green
  "spei_maturation" = "#e5f5e0",# pale green
  
  # Soil properties (orange → brown)
  "Clay_p" = "#e6550d",         # dark orange
  "SOC_p" = "#fd8d3c",          # medium orange
  "Sand_p2" = "#fdae6b",        # light orange
  "BD_p" = "#8c510a",           # brown
  "PAW_p" = "#d8b365",          # tan/orange-brown
  
  # Others
  "pH_p" = "#ffff99",           # yellow
  "Unbounded" = "#000000"       # black
)
# Plot
ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(
    data = pts,
    aes(color = vector_name),
    size = 1.75,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = my_colours,
    drop = FALSE,    
    name = "Limiting factor"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Limiting Factors Map",
    x = "Longitude",
    y = "Latitude"
  )


library(stringr)
pts2 <- pts %>%
  mutate(
    vector_simple = case_when(
      str_starts(vector_name, "rain_") ~ "Rainfall",
      str_starts(vector_name, "spei_") ~ "SPEI",
      TRUE ~ vector_name
    ),
    vector_simple = str_remove(vector_simple, "_p$")
  )

simple_cols <- c(
  "Radiation" = "#E69F00",   # amber / gold
  "Rainfall" = "#2A9D8F",   # strong blue
  "SPEI"     = "lightgreen",   # green
  "Clay"     = "#A63603",   # red-brown
  "SOC"      = "#5C4033",   # dark brown
  "Sand"     = "#F4A261",   # sandy orange
  "BD"       = "#8D6E63",   # grey-brown
  "PAW"      = "lightblue",   # teal
  "pH"       = "#F1C40F",   # yellow
  "Unbounded"= "#000000"    # black
)

ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(
    data = pts2,
    aes(color = vector_simple),
    size = 1.75,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = simple_cols,
    name = "Limiting factor"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Limiting Factors Map",
    x = "Longitude",
    y = "Latitude"
  )


# checking sow dates
ggplot() +
  geom_sf(data = Africa, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(
    data = pts,
    aes(color = planting_month), #vector_name
    size = 2,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "1"    = "grey",
      "2"    = "blue",
      "3"    = "green",
      "4"    = "purple",
      "5"    = "red",
      "6"    = "black",
      "7"    = "yellow",
      "8"    = "orange",
      "9"    = "pink",
      "10"    = "darkblue",
      "11"    = "darkgreen",
      "12"    = "lightblue"
    ),
    name = "Limiting factor"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Limiting Factors Map",
    x = "Longitude",
    y = "Latitude"
  )

aggregate(min_value ~ Country, data = data_pred, mean, na.rm = TRUE)



aggregate(min_value ~ Country, data = data_pred,FUN = function(x) c(mean = mean(x, na.rm = TRUE),n = sum(!is.na(x))))

























