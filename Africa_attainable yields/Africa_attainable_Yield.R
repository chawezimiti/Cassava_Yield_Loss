
library(readr)
library(dplyr)
library(BLA)
library(aplpack)
library(countrycode)
library(terra) # mapping
library(sf) # mapping
source("Extra_functions.R")


## Read in data-----------------------------------------------------------------

#data <- read_csv("master_dataset_v14.csv")
data <- read_csv("03-merged_abbr_rainfall_tmin_tmax_(v2).csv")
names(data)

# 1) Data manipulation and Feature engineering ----------------------------------------------------------

# Select only African countries: applies only to the master data set

# data$continent <- countrycode(data$country,
#                               origin = "country.name",
#                               destination = "continent")
# 
# levels(as.factor(data$country))
# levels(as.factor(data$continent))
# levels(as.factor(data[which(data$continent=="Africa"),]$country))
# 
# data <- data[which(data$continent=="Africa"),]
# 
# Af <- vect("Africa_Countries.shp") # Tanzania shape file for extents
# Af <- project(Af, "EPSG:4326") #optional, convert to decimal degrees
# plot(Af)
# points(data$decimal_longitude, data$decimal_latitude, cex=0.5, pch=16, col="red")
# 
# plot(data$decimal_latitude,data$cassava_fr_root_yld_tha,  cex=0.5, pch=16, col="red")

# Here we aggregate some of the data according to the year for climatic variables 
# and by depth for soil properties.


which(is.na(data$experiment_year))# check for NA in the experiment year, to be used to aggregate climatic variables at site.
unique(data$experiment_year)

data$experiment_year <- ifelse(
  nchar(data$experiment_year) > 4,
  substr(data$experiment_year, 1, 4),
  data$experiment_year
) # this step allows us to retain the only the first 4 characters for experiment years that have more than 4 digits e.g 2002-2004 will be 2002

data <- data %>% mutate(experiment_year = as.integer(experiment_year)) #converts the experimental year column to integer
years <- sort(unique(na.omit(data$experiment_year)))

length(which(is.na(data$planting_month)))# check for NA in the experiment year
unique(data$planting_month) # check if the values are sensible between 1 and 12

## a) Cumulative rainfall----------------------------------------------------------

rain_cols <- grep("^(?:pmm|ppm)_\\d{4}_\\d{2}$",
                  names(data),
                  value = TRUE,
                  ignore.case = TRUE)

if (length(rain_cols) == 0L) {
  data$pmm_sum_12mo <- NA_real_
} else {
  ## 2) Parse year + month from column names
  rx <- "^(?:pmm|ppm)_(\\d{4})_(\\d{2})$"
  mm_yy <- do.call(
    rbind,
    regmatches(rain_cols, regexec(rx, rain_cols, ignore.case = TRUE))
  )
  # mm_yy columns: full match, year, month
  yr_int  <- as.integer(mm_yy[, 2])
  mon_int <- as.integer(mm_yy[, 3])
  
  # Use first of each month for date indexing
  rain_dates <- as.Date(sprintf("%04d-%02d-01", yr_int, mon_int))
  
  # Order both dates and column names chronologically
  ord <- order(rain_dates)
  rain_dates <- rain_dates[ord]
  rain_cols  <- rain_cols[ord]
  
  ## 3) Convert rainfall values into a numeric matrix
  pmm_mat <- as.matrix(data[rain_cols])
  storage.mode(pmm_mat) <- "double"
  
  ## 4) Precompute 12-month windows
  data$experiment_year <- as.integer(data$experiment_year)
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
  
  ## 5) Compute 12-month rainfall total for each row
  keys <- paste0(
    data$experiment_year, "_",
    sprintf("%02d", as.integer(data$planting_month))
  )
  idx_list <- idx_by_yrmon[keys]
  
  data$pmm_yr_data <- mapply(
    function(i, idx) {
      if (length(idx) == 0L) return(NA_real_)
      vals <- pmm_mat[i, idx]
      if (all(is.na(vals))) return(NA_real_)
      sum(vals, na.rm = TRUE)
    },
    seq_len(nrow(data)), idx_list
  )
}


plot(data$pmm_yr_data, data$cassava_fr_root_yld_tha)

#-------------new update

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
      
      p1_2  <- sum_range(0, 0)   # months 1–2
      p3_4  <- sum_range(1, 3)   # months 3–4
      p5_8  <- sum_range(4, 8)   # months 5–8
      p9_12 <- sum_range(9, 11)  # months 9–12
      
      c(p1_2, p3_4, p5_8, p9_12)
    },
    seq_len(nrow(data)), idx_list,
    SIMPLIFY = TRUE
  ))
  
  data$Establishment  <- res_mat[, 1]
  data$Canopy  <- res_mat[, 2]
  data$Bulking  <- res_mat[, 3]
  data$Maturation <- res_mat[, 4]
}

plot(data$Establishment, sqrt(data$cassava_fr_root_yld_tha))
plot(data$Canopy, sqrt(data$cassava_fr_root_yld_tha))
plot(data$Bulking, sqrt(data$cassava_fr_root_yld_tha))
plot(data$Maturation, sqrt(data$cassava_fr_root_yld_tha))

## b) Temperature Indices -------------------------------------------------------

# here we create Growing degree days, Heat stress degree days and chilling degree days

data$experiment_year  <- as.integer(data$experiment_year)
data$planting_month   <- as.integer(data$planting_month)

data$planting_date <- as.Date(
  sprintf("%04d-%02d-01", data$experiment_year, data$planting_month)
)

## 1. Identify daily Tmin/Tmax columns and common dates -------------------

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
  
  ## 2. Align daily temperature matrices ----------------------------------
  
  tmin_common_cols <- paste0("tmin_", format(common_dates, "%Y_%m_%d"))
  tmax_common_cols <- paste0("tmax_", format(common_dates, "%Y_%m_%d"))
  
  tmin_mat <- as.matrix(data[tmin_common_cols]); storage.mode(tmin_mat) <- "double"
  tmax_mat <- as.matrix(data[tmax_common_cols]); storage.mode(tmax_mat) <- "double"
  tmean_mat <- (tmin_mat + tmax_mat) / 2
  
  ## 3. Thresholds + stages (cassava) -------------------------------------
  
  Tbase_gdd   <- 10  # GDD base
  Tcrit_heat  <- 32  # heat-stress threshold (32 original)
  Tbase_chill <- 18 # cold-stress threshold (15 original)
  
  # stages in days after planting (DAP); upper bound is exclusive
  stages <- list(
    early       = c(0,   60),   # 0–60 DAP
    vegetative  = c(60,  120),  # 60–120 DAP
    bulking     = c(120, 240),  # 120–240 DAP
    maturation  = c(240, 365)   # 240–365 DAP
  )
  
  
  ## 4. Row-wise calculation of GDD, HSDD, CDD by stage -------------------
  
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
  
  ## 5. Attach back to data -----------------------------------------------
  
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


plot(data$GDD_early , sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$GDD_vegetative, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$GDD_bulking, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$GDD_maturation , sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$GDD_total, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")

plot(data$HSDD_early, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey", xlim=c(0,20))
plot(data$HSDD_vegetative, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$HSDD_bulking , sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$HSDD_maturation, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$HSDD_total, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")

plot(data$CDD_early, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$CDD_vegetative, sqrt(data$cassava_fr_root_yld_tha), pch=16, col="grey")
plot(data$CDD_bulking, sqrt(data$cassava_fr_root_yld_tha), xlim=c(0,20), pch=16, col="grey")
plot(data$CDD_maturation, sqrt(data$cassava_fr_root_yld_tha), xlim=c(0,20), pch=16, col="grey")
plot(data$CDD_total, sqrt(data$cassava_fr_root_yld_tha), xlim=c(0,20), pch=16, col="grey")


length(which(is.na(data$planting_month)==TRUE))/length(data$planting_month)
length(which(!is.na(data$planting_month)==TRUE))/length(data$planting_month)
length(which(is.na(data$experiment_year)==TRUE))/length(data$experiment_year)


## c) standard precipitation-evaporation index ---------------------------------

# Determine the average spei index over the growing period of 12 months

## 1) Identify monthly SPEI columns and parse their dates

spei_cols  <- grep("^spei_\\d{4}-\\d{2}-\\d{2}$", names(data), value = TRUE)

if (length(spei_cols) == 0L) {
  data$spei_mean_12mo <- NA_real_
} else {
  spei_dates <- as.Date(sub("^spei_", "", spei_cols), format = "%Y-%m-%d")
  
  ## Order by calendar time and align columns
  ord <- order(spei_dates)
  spei_dates <- spei_dates[ord]
  spei_cols  <- spei_cols[ord]
  
  ## 2) Build a matrix of monthly SPEI (force numeric)
  spei_mat <- as.matrix(data[spei_cols])
  storage.mode(spei_mat) <- "double"
  
  ## 3) Pre-compute indices for each (year, month) 12-month window
  data$experiment_year <- as.integer(data$experiment_year)
  years  <- sort(unique(data$experiment_year))
  months <- sprintf("%02d", 1:12)
  
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
  
  ## 4) Average SPEI over the 12-month window for each row
  keys <- paste0(
    data$experiment_year, "_",
    sprintf("%02d", as.integer(data$planting_month))  # robust if planting_month is char
  )
  idx_list <- idx_by_yrmon[keys]
  
  data$spei_yr_data <- mapply(
    function(i, idx) {
      if (length(idx) == 0L) return(NA_real_)
      vals <- spei_mat[i, idx]
      if (all(is.na(vals))) return(NA_real_)
      mean(vals, na.rm = TRUE)
    },
    seq_len(nrow(data)), idx_list
  )
}

plot(data$spei_yr_data, data$cassava_fr_root_yld_tha)


## c) Latitude-Yield -----------------------------------------------------------

# The attainable yield are modeled based on the latitude. The distance away from 
# the latitude 0. To check the influence of climate on cassava yields.

# Create a data frame and standardize the missing  values (either NA or blank spaces) to NA

dat <- data.frame(x=data$decimal_latitude,y=data$cassava_fr_root_yld_tha)
dat[dat == ""] <- NA
dat  <- na.omit(dat) # omit all rows with the NA values

plot(dat)

# Summary statistics------------------------------------------------------------

summastat(dat[,2]) # the plots do not look normally distributed. We try transforming
summastat(sqrt(dat[,2])) #looks better

dat$y <- sqrt(dat[,2])

## Bagplot outliers ------------------------------------------------------------

out <- bagplot(dat, show.whiskers=F, na.rm = TRUE)
dat_clean<- rbind(out$pxy.bag, out$pxy.outer)
plot(dat_clean)

{par(mfrow=c(1,2))
  plot(dat, ylim=c(0,10), pch=16)
  abline(h=max(dat_clean[,2]), col="red")
  plot(dat_clean, ylim=c(0,10), pch=16)
  par(mfrow=c(1,1))}



# Explore boundary presence ----------------------------------------------------

expl_boundary(x=dat_clean[,1], y=dat_clean[,2], method = "Area") # no evidence of boundary

# Fitting the boundary line model ----------------------------------------------

# a) trapezium model------------------------------------------------------------

# determine Start values for the boundary line model 

plot(dat_clean)
startValues("trapezium")

# start2 <- list(c(9.47, 0.20, 8.02, 14.44,-0.42, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),        
#                c(9.13, 0.19, 8.07, 14.81,-0.44, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),
#                c(9.10, 0.20, 8.11, 13.96,-0.31, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")))


start2 <- list(c(12.35, 0.41, 8.00, 13.33, -0.61, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),        
               c(11.46, 0.34, 8.03 , 13.35,-0.59, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),
               c(11.30, 0.31, 8.06, 13.47, -0.66, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")))
# determining of the standard deviation of measurement error

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6)
ble_extention(data = dat_clean, start = start2, sigh = sigh, model = "trapezium")

# model fitting 

par(mar=c(6,5,4,2)) # plot size parameter adjustment

cbvn_extention(data=dat_clean, start = start2, sigh = 0.6, model = "trapezium", pch=16, col="grey",
               ylab=expression(bold("Yield / t ha"^-1)), 
               xlab=expression(bold("Latitude")), main="trapezium Model")

## b) Mirror trapezium model --------------------------------------------------

# This function creates a trapezium function that is mirrored at the equator

trap_mirror <- function(x, a, b, c) {
  # x : predictor variable
  # a : inflection point (distance from 0 where plateau starts)
  # b : slope (positive number)
  # c : plateau value (maximum)
  
  d <- abs(x) - a                 # distance from plateau edge
  y <- ifelse(d <= 0, c, c - b * d)
  return(y)
}

plot(dat_clean)
startValues("trapezium") # use the trapezium to set the starting values (the slope and max value) for the trapezium mirror

start2 <- list(c(10, 0.42, 8.02, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),        
               c(11, 0.44, 8.07, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),
               c(10.5, 0.31, 8.21, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")))


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) #possible measurement error values

ble_extention(data = dat_clean, start = start2, sigh = sigh, model = "other", equation=trap_mirror) #determine measurement error by maximum likelihood


cbvn_extention(data = dat_clean, start = start2, sigh = 0.6, model = "other", equation=trap_mirror, pch=16, col="grey",
               ylab=expression(bold("Yield / t ha"^-1)), 
               xlab=expression(bold("Latitude")),
               main="Mirror Model")



