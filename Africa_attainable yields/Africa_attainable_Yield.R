
library(readr)
library(dplyr)
library(BLA)
library(aplpack)
library(countrycode)
library(terra) # mapping
library(sf) # mapping
source("Extra_functions.R")


## 1) Read in data-----------------------------------------------------------------

data <- read_csv("03-merged_abbr_rainfall_tmin_tmax_(v3).csv")
names(data)

# 1) Data manipulation and Feature engineering ================================================================================

length(which(is.na(data$experiment_year)))# check for NA in the experiment year. To be used to aggregate climatic variables at site.
unique(data$experiment_year)

# The next step allows us to retain the only the first 4 characters for experiment years that have more than 4 digits e.g 2002-2003 will be 2002

data$experiment_year <- ifelse(
  nchar(data$experiment_year) > 4,
  substr(data$experiment_year, 1, 4),
  data$experiment_year
) 

data <- data %>% mutate(experiment_year = as.integer(experiment_year)) #converts the experimental year column to integer
years <- sort(unique(na.omit(data$experiment_year)))

length(which(is.na(data$planting_month)))# check for NA in the Planting month. PM is used as start point for aggregation of climatic variables
unique(data$planting_month) # check if the values for PM are sensible between 1 and 12

## a) Cumulative Precipitation ==========================================================

# The Cumulative rainfall is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), vegetative stage (2-4 months),
# bulking (5-8 months) and maturation (9-12 months).

rain_cols <- grep("^(?:pmm|ppm)_\\d{4}_\\d{1,2}$",
                  names(data),
                  value = TRUE,
                  ignore.case = TRUE)

if (length(rain_cols) == 0L) {
  data$pmm_1_2  <- NA_real_
  data$pmm_3_4  <- NA_real_
  data$pmm_5_8  <- NA_real_
  data$pmm_9_12 <- NA_real_
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
  pmm_mat <- as.matrix(data[rain_cols])
  storage.mode(pmm_mat) <- "double"
  
  ## 4) Precompute 12-month windows by (experiment_year, planting_month)
  data$experiment_year <- as.integer(
    ifelse(nchar(data$experiment_year) > 4,
           substr(data$experiment_year, 1, 4),
           data$experiment_year)
  )
  years  <- sort(unique(data$experiment_year))
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
  
  ## 5) Compute rainfall for 1–2, 3–4, 5–8, 9–12 months after planting
  keys <- paste0(
    data$experiment_year, "_",
    sprintf("%02d", as.integer(data$planting_month))
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
      yr <- data$experiment_year[i]
      pm <- as.integer(data$planting_month[i])
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
      
      # p1_2  <- sum_range(0, 1)   # months 1–2
      # p3_4  <- sum_range(2, 3)   # months 3–4
      # p5_8  <- sum_range(4, 7)   # months 5–8
      # p9_12 <- sum_range(8, 11)  # months 9–12
      
      p1_2  <- sum_range(0, 0)   # months 1
      p3_4  <- sum_range(1, 3)   # months 2–4
      p5_8  <- sum_range(4, 8)   # months 5–8
      p9_12 <- sum_range(9, 11)  # months 9–12
      
      c(p1_2, p3_4, p5_8, p9_12)
    },
    seq_len(nrow(data)), idx_list,
    SIMPLIFY = TRUE
  ))
  
  data$pmm_early  <- res_mat[, 1]
  data$pmm_vegetative  <- res_mat[, 2]
  data$pmm_bulking  <- res_mat[, 3]
  data$pmm_maturation <- res_mat[, 4]
}

# First exploration of plots of yield~rainfall by regions of Africa (East, West, South, central)
regions <- unique(data$region)
cols <- c("red", "blue", "green", "yellow", "black")
region_col <- cols[ match(data$region, regions) ]

# i) longitude~latitude
plot(data$decimal_longitude , data$decimal_latitude, col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# ii) cassava yield~longitude
plot(data$decimal_longitude , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iii) cassava yield~latitude
plot(data$decimal_latitude , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iv) cassava yield ~ rainfall (establishment)
plot(data$pmm_early , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# v) cassava yield ~ rainfall (vegetation)
plot(data$pmm_vegetative, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# vi) cassava yield ~ rainfall (bulking)
plot(data$pmm_bulking, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# vii) cassava yield ~ rainfall (maturation)
plot(data$pmm_maturation, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))


## b) Temperature Indices =============================================================

# Here we create Growing degree days, Heat stress degree days and chilling degree days.
# The GDD, HSDD and CDD is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), vegetative stage (2-4 months),
# bulking (5-8 months) and maturation (9-12 months).

data$planting_date <- as.Date(
  sprintf("%04d-%02d-01", data$experiment_year, data$planting_month)
)

## Identify daily Tmin/Tmax columns and common dates ---------------------------

tmin_cols  <- grep("^tmin_\\d{4}_\\d{2}_\\d{2}$", names(data), value = TRUE)
tmax_cols  <- grep("^tmax_\\d{4}_\\d{2}_\\d{2}$", names(data), value = TRUE)

tmin_dates <- as.Date(sub("^tmin_", "", tmin_cols), format = "%Y_%m_%d")
tmax_dates <- as.Date(sub("^tmax_", "", tmax_cols), format = "%Y_%m_%d")

common_dates <- sort(as.Date(intersect(tmin_dates, tmax_dates)))

if (length(common_dates) == 0L) {
  # no climatic overlap → everything NA
  data$GDD_early          <- NA_real_
  data$GDD_vegetative     <- NA_real_
  data$GDD_bulking        <- NA_real_
  data$GDD_maturation     <- NA_real_
  
  data$HSDD_early         <- NA_real_
  data$HSDD_vegetative    <- NA_real_
  data$HSDD_bulking       <- NA_real_
  data$HSDD_maturation    <- NA_real_
  
  data$CDD_early          <- NA_real_
  data$CDD_vegetative     <- NA_real_
  data$CDD_bulking        <- NA_real_
  data$CDD_maturation     <- NA_real_
  
} else {
  
  ## Align daily temperature matrices ----------------------------------
  
  tmin_common_cols <- paste0("tmin_", format(common_dates, "%Y_%m_%d"))
  tmax_common_cols <- paste0("tmax_", format(common_dates, "%Y_%m_%d"))
  
  tmin_mat <- as.matrix(data[tmin_common_cols]); storage.mode(tmin_mat) <- "double"
  tmax_mat <- as.matrix(data[tmax_common_cols]); storage.mode(tmax_mat) <- "double"
  tmean_mat <- (tmin_mat + tmax_mat) / 2
  
  ## Thresholds + stages (cassava) -------------------------------------
  
  Tbase_gdd   <- 10  # GDD base
  Tcrit_heat  <- 32  # heat-stress threshold 
  Tbase_chill <- 18 # cold-stress threshold 
  
  # stages in days after planting (DAP); upper bound is exclusive
  stages <- list(
    early       = c(0,   30),   # 0–60 DAP
    vegetative  = c(30,  120),  # 60–120 DAP
    bulking     = c(120, 240),  # 120–240 DAP
    maturation  = c(240, 365)   # 240–365 DAP
  )
  
  
  ## Row-wise calculation of GDD, HSDD, CDD by stage -------------------
  
  res <- mapply(
    function(i) {
      plant_date <- data$planting_date[i]
      
      if (is.na(plant_date)) {
        return(rep(NA_real_, 12L))
      }
      
      out <- numeric(12L)
      names(out) <- c(
        "GDD_early", "GDD_vegetative", "GDD_bulking", "GDD_maturation",
        "HSDD_early","HSDD_vegetative","HSDD_bulking","HSDD_maturation",
        "CDD_early", "CDD_vegetative", "CDD_bulking", "CDD_maturation"
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
    seq_len(nrow(data))
  )
  
  res <- t(res)
  
  ## Attach back to data -----------------------------------------------
  
  data$GDD_early          <- res[, "GDD_early"]
  data$GDD_vegetative     <- res[, "GDD_vegetative"]
  data$GDD_bulking        <- res[, "GDD_bulking"]
  data$GDD_maturation     <- res[, "GDD_maturation"]
  data$GDD_total <- rowSums(data[, c("GDD_early","GDD_vegetative","GDD_bulking","GDD_maturation")],na.rm = TRUE)
  
  data$HSDD_early         <- res[, "HSDD_early"]
  data$HSDD_vegetative    <- res[, "HSDD_vegetative"]
  data$HSDD_bulking       <- res[, "HSDD_bulking"]
  data$HSDD_maturation    <- res[, "HSDD_maturation"]
  data$HSDD_total <- rowSums(data[, c("HSDD_early","HSDD_vegetative","HSDD_bulking","HSDD_maturation")], na.rm = TRUE)
  
  data$CDD_early          <- res[, "CDD_early"]
  data$CDD_vegetative     <- res[, "CDD_vegetative"]
  data$CDD_bulking        <- res[, "CDD_bulking"]
  data$CDD_maturation     <- res[, "CDD_maturation"]
  data$CDD_total  <- rowSums(data[, c("CDD_early","CDD_vegetative","CDD_bulking","CDD_maturation")], na.rm = TRUE)
}

# First exploration of plots of yield~rainfall by regions of Africa (East, West, South, central)

# i) Yield~establishment (0-1 months)
plot(data$GDD_early , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# ii) Yield~vegetative (2-4 months)
plot(data$GDD_vegetative, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iii) Yield~bulking (5-8 months)
plot(data$GDD_bulking, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iv) Yield~maturation (9-12 months)
plot(data$GDD_maturation , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))
#plot(data$GDD_total, sqrt(data$cassava_fr_root_yld_tha), pch=16, col=region_col)

# v) Yield~establishment (0-1 months)
plot(data$HSDD_early, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16, xlim=c(0,20))
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# vi) Yield~vegetative (2-4 months)
plot(data$HSDD_vegetative, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topright", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# vii) Yield~bulking (5-8 months)
plot(data$HSDD_bulking , sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# viii) Yield~maturation (9-12 months)
plot(data$HSDD_maturation, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))
#plot(data$HSDD_total, sqrt(data$cassava_fr_root_yld_tha), pch=16, col=region_col)


## c) standard precipitation-evaporation index===========================================

# Determine the average spei index over the growing period of 12 months
# The spei is determined for the 4 critical stages of cassava growth
# Which include the i) establishment stage (0-1 month), vegetative stage (2-4 months),
# bulking (5-8 months) and maturation (9-12 months).

spei_cols  <- grep("^spei_\\d{4}-\\d{2}-\\d{1,2}$", names(data), value = TRUE)## Identify monthly SPEI columns with dates like spei_2007-10-2

if (length(spei_cols) == 0L) {
  data$spei_early  <- NA_real_
  data$spei_vegetative  <- NA_real_
  data$spei_bulking  <- NA_real_
  data$spei_maturation <- NA_real_
  
} else {
  
  ## Convert to Date
  spei_dates <- as.Date(sub("^spei_", "", spei_cols), format = "%Y-%m-%d")
  
  ## Order by date and align columns
  ord <- order(spei_dates)
  spei_dates <- spei_dates[ord]
  spei_cols  <- spei_cols[ord]
  
  ## SPEI matrix
  spei_mat <- as.matrix(data[spei_cols])
  storage.mode(spei_mat) <- "double"
  
  ## Clean experiment year (handle "2013-2016")
  data$experiment_year <- ifelse(
    nchar(data$experiment_year) > 4,
    substr(data$experiment_year, 1, 4),
    data$experiment_year
  )
  data$experiment_year <- as.integer(data$experiment_year)
  
  years  <- sort(unique(data$experiment_year))
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
  data$spei_early <- mapply(function(i, yr)
    get_period_mean(i, yr, 0, 1),
    seq_len(nrow(data)), data$experiment_year)
  
  data$spei_vegetative <- mapply(function(i, yr)
    get_period_mean(i, yr, 2, 4),
    seq_len(nrow(data)), data$experiment_year)
  
  data$spei_bulking <- mapply(function(i, yr)
    get_period_mean(i, yr, 5, 8),
    seq_len(nrow(data)), data$experiment_year)
  
  data$spei_maturation <- mapply(function(i, yr)
    get_period_mean(i, yr, 9, 12),
    seq_len(nrow(data)), data$experiment_year)
}

# First exploration of plots of yield~rainfall by regions of Africa (East, West, South, central)

# i) Yield~establishment (0-1 months)
plot(data$spei_early, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# ii) Yield~vegetation (2-4 months)
plot(data$spei_vegetative, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iii) Yield~bulking (5-8 months)
plot(data$spei_bulking, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# iv) Yield~maturation (9-12 months)
plot(data$spei_maturation, sqrt(data$cassava_fr_root_yld_tha), col=region_col, pch=16)
legend("topleft", legend=regions, pch=16, col=c("red", "blue", "green", "yellow"))

# 1) Boundary line model fitting  ==============================================

## a) Latitude-Yield ===========================================================

# The attainable yield are modeled based on the latitude. The distance away from 
# the equator, latitude 0. To check the influence of climate on cassava yields.

dat <- data.frame(x=data$decimal_latitude,y=data$cassava_fr_root_yld_tha) 
dat[dat == ""] <- NA # standardize missing  values (NA or blank spaces) to NA
dat  <- na.omit(dat) # omit all rows with the NA values
plot(dat, pch=16, col="grey")

## i) Summary statistics-----------
summastat(dat[,2]) # the plots do not look normally distributed. We try transforming
summastat(sqrt(dat[,2])) #looks better
dat$y <- sqrt(dat[,2])

## ii) Bagplot outlier detection and removal -----------
out <- bagplot(dat, show.whiskers=F, na.rm = TRUE) # identifies outliers
dat_clean<- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

## iii) Explore boundary structure in the data cloud ---------
expl_boundary(x=dat_clean[,1], y=dat_clean[,2], method = "Area") 

## iv) Define the boundary line model ----------------

trap_mirror <- function(x, a, b, c) { 
  # a=inflection point (distance from 0 where plateau starts), b=slope (positive),c=plateau value
  d <- abs(x) - a 
  y <- ifelse(d <= 0, c, c - b * d)
  return(y)
}

## iv) Initial start values for optimizing parameters of the bivariate distribution and the mirror model

# First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y

dist_lat <- c(mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), 
              sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs"))

plot(dat_clean, pch=16, col="grey")
startValues("trapezium") # Determine the start values for the mirror trapezium model i.e inflection point, slope and plateau value

# determine multiple start values
start_lat <- list(c(10, 0.42, 8.02, dist_lat), c(11, 0.44, 8.07, dist_lat), c(10.5, 0.31, 8.21, dist_lat))

## v) Determine the sigh value (Related to sd of measurement error), which a hyper-parameter of censored normal model

ble_extention(data = dat_clean, start = start_lat, sigh = sigh, model = "other", equation=trap_mirror) #determine sigh by maximum likelihood

## v) Fit the boundary model

model_lat <- cbvn_extention(data = dat_clean, start = start_lat, sigh = 0.3, model = "other", equation=trap_mirror, 
                            pch=16, col="grey", ylab=expression(bold("Yield / sqrt(t ha"^-1*")")), xlab=expression(bold("Latitude")),
                            main="Mirror boundary Model", optim.method = "Nelder-Mead")


## b) Rainfall-Yield Model ===========================================================

# we fit the boundary model for establishment (1-2 months), vegetative growth (2-4 months),
# root bulking (5-9 months) and maturation (9-12 months)

### b.1. Establishment =========================================================

dat_est <- data.frame(x=data$pmm_early, y=data$cassava_fr_root_yld_tha)

## i) Distribution properties of data-------

summastat(dat_est[,1])
summastat(log(dat_est[,1]))
summastat(dat_est[,2]) # the plots do not look normally distributed. We try transforming
summastat(sqrt(dat_est[,2])) #looks better

dat_est$x <- log(dat_est[,1])
dat_est$y <- sqrt(dat_est[,2])

## ii) Bagplot outlier identification and removal-----

out <- bagplot(dat_est, show.whiskers=F, na.rm = TRUE)
dat_est_clean <- rbind(out$pxy.bag, out$pxy.outer)

## iii) Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_pmm_est <- c(mean(dat_est_clean[,1], na.rm = T), mean(dat_est_clean[,2], na.rm = T), sd(dat_est_clean[,1], na.rm = T), 
                  sd(dat_est_clean[,2], na.rm = T), cor(dat_est_clean[,1],dat_est_clean[,2], use = "complete.obs"))

plot(dat_est_clean, pch=16, col="grey")
startValues("lp") # Determine the start values for the linear-plateau model i.e intercept, slope and plateau value
start_pmm_est <- list(c(25.45, -3.32, 8.17, dist_pmm_est), c(27.61, -3.67, 8.02, dist_pmm_est), c(27.52, -3.65, 7.89, dist_pmm_est))

## iv) Fit the boundary model---------

# we use the sigh from the yield~latitude model

model_pmm_est <- cbvn_extention(data = dat_est_clean, start = start_pmm_est, sigh = 0.3, model = "lp", 
                                  pch=16, col="grey",ylab=expression(bold("Yield / t ha"^-1)), 
                                  xlab=expression(bold("Rainfall (0-1 months)/mm")))

### b.2. vegetative growth (Canopy growth)===========================================

dat_veg <- data.frame(x=data$pmm_vegetative, y=data$cassava_fr_root_yld_tha)

## i) Distribution properties of data-------

summastat(dat_veg[,1])
summastat(log(dat_veg[,1]))
summastat(dat_veg[,2]) # The plots do not look normally distributed. We try transforming
summastat(sqrt(dat_veg[,2])) #looks better
dat_veg$y <- sqrt(dat_veg[,2])

## ii) Bagplot outlier identification and removal-----
out <- bagplot(dat_veg, show.whiskers=F, na.rm = TRUE)
dat_veg_pmm <- rbind(out$pxy.bag, out$pxy.outer)

## iii) Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_pmm_veg <- c(mean(dat_veg_pmm[,1], na.rm = T), mean(dat_veg_pmm[,2], na.rm = T), sd(dat_veg_pmm[,1], na.rm = T), 
                  sd(dat_veg_pmm[,2], na.rm = T), cor(dat_veg_pmm[,1],dat_veg_pmm[,2], use = "complete.obs"))

plot(dat_veg_pmm, pch=16, col="grey")
startValues("trapezium") # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_pmm_veg <- list(c(6.06, 0.011, 8.34,15.20,-0.011, dist_pmm_veg),c(5.87, 0.011, 8.29,15.04,-0.011, dist_pmm_veg),
                  c(5.69, 0.015, 8.35,14.31,-0.010, dist_pmm_veg), c(5.04, 0.015, 8.35,14.31,-0.010, dist_pmm_veg),
                  c(5.48, 0.018, 8.19,15.99,-0.011, dist_pmm_veg))

## iv) Fit the boundary model---------

# we use the sigh = 0.3 from the yield~latitude model

model_pmm_veg <- cbvn_extention(data = dat_veg_pmm, start = start_pmm_veg, sigh = 0.3, model = "trapezium", 
                                pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                                xlab=expression(bold("Rainfall (1-4 months)/mm")))


### b.3. Root Bulking ============================================================

dat_bulk <- data.frame(x=data$pmm_bulking, y=data$cassava_fr_root_yld_tha)

## i) Distribution properties of data-------

summastat(dat_bulk[,1]) 
summastat(log(dat_bulk[,1])) 
summastat(sqrt(dat_bulk[,1])) 
summastat(dat_bulk[,2]) # the plots do not look normally distributed. We try transforming
summastat(sqrt(dat_bulk[,2])) #looks better

dat_bulk$x <- log(dat_bulk[,1])
dat_bulk$y <- sqrt(dat_bulk[,2])

## ii) Bagplot outlier identification and removal-----

out <- bagplot(dat_bulk, show.whiskers=F, na.rm = TRUE)
dat_bulk_clean<- rbind(out$pxy.bag, out$pxy.outer)

## iii) Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_bulk <- c(mean(dat_bulk_clean[,1], na.rm = T), mean(dat_bulk_clean[,2], na.rm = T), sd(dat_bulk_clean[,1], na.rm = T), 
              sd(dat_bulk_clean[,2], na.rm = T), cor(dat_bulk_clean[,1],dat_bulk_clean[,2], use = "complete.obs"))

plot(dat_bulk_clean, pch=16, col="grey")
startValues("trapezium") # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_bulk <- list(c(6.06, 0.011, 8.34, dist_bulk), c(5.87, 0.011, 8.29, dist_bulk), c(5.48, 0.018, 8.19, dist_bulk))

## iv) Fit the boundary model---------

# we use the sigh = 0.3 from the yield~latitude model

model_pmm_bulk <- cbvn_extention(data = dat_bulk_clean, start = start_bulk, sigh = 0.3, model = "lp", pch=16, col="grey",
                                ylab=expression(bold("Yield / t ha"^-1)), xlab=expression(bold("Rainfall (4-9 months)/mm")), 
                                optim.method = "Nelder-Mead")

### b.4. Maturation===============================================================

dat_Mat <- data.frame(x=data$pmm_maturation, y=data$cassava_fr_root_yld_tha)

## i) Distribution properties of data-------
summastat(dat_Mat[,1])
summastat(sqrt(dat_Mat[,1]))
summastat(dat_Mat[,2]) # the plots do not look normally distributed. We try transforming
summastat(sqrt(dat_Mat[,2])) #looks better

dat_Mat$x <- sqrt(dat_Mat[,1])
dat_Mat$y <- sqrt(dat_Mat[,2])

## ii) Bagplot outlier identification and removal-----
out <- bagplot(dat_Mat, show.whiskers=F, na.rm = TRUE)
dat_Mat_clean<- rbind(out$pxy.bag, out$pxy.outer)

## iii) Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_Mat <- c(mean(dat_Mat_clean[,1], na.rm = T), mean(dat_Mat_clean[,2], na.rm = T), sd(dat_Mat_clean[,1], na.rm = T), 
              sd(dat_Mat_clean[,2], na.rm = T), cor(dat_Mat_clean[,1],dat_Mat_clean[,2], use = "complete.obs"))

plot(dat_Mat_clean, pch=16, col="grey")
startValues("trapezium")# Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_Mat <- list(c(2.68, 0.38, 8.01,12.56,-0.23,dist_Mat),c(1.98, 0.46, 8.17,12.65,-0.23,dist_Mat),c(2.39, 0.42, 8.04,12.30,-0.21, dist_Mat))

## iv) Fit the boundary model---------

# we use the sigh = 0.3 from the yield~latitude model

model_pmm_mat <- cbvn_extention(data = dat_Mat_clean, start = start_Mat, sigh = 0.3, model = "trapezium", 
                                 pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                                 xlab=expression(bold("Rainfall (9-12 months)/mm")))



