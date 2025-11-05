# This script is meant to carryout global yield gap analysis for cassava under 
# the CABI project. We employ the boundary line methodology (Webb, 1972) for this purpose.
# It used climatic, soil and management factors to determine the attainable yield 
# at country level as well as the most limiting factor.


# some useful information on cassava can be found here https://greg.app/cassava-lifecycle/

#------ Install and Load necessary libraries and Source files ------------------------------

library(devtools)
install_github("chawezimiti/BLA") # installs the development version of BLA


library(readr)
library(readxl)
library(dplyr)
library(BLA)
library(aplpack)
library(stringr)
library(terra) # mapping
library(sf) # mapping
library(randomcoloR)
source("Extra_functions.R") 

#---- 1. Data loading into R ===================================================

data <- read_csv("01-all_tanzania_daily_temp(v2).csv")

head(data)

#--- 2. Data manipulation and cleaning =======================================

# Here we aggregate some for the data according to the year for climatic variables 
# and by depth for soil properties. 


which(is.na(data$experiment_year))# check for NA in the experiment year, to be used to aggregate climatic variables at site.

unique(data$experiment_year)
data$experiment_year <- ifelse(data$experiment_year=="2013-2016", "2013", data$experiment_year)
data <- data %>% mutate(experiment_year = as.integer(experiment_year)) #converts the experimental year column to integer
years <- sort(unique(na.omit(data$experiment_year)))


# a) Growing Degree days

# This is calculated as sum(max((Tmim + Tmax)*0.5 -Tbase, 0)) over the growing period (12 months) from planting date 

tmin_cols  <- grep("^tmin_\\d{4}_\\d{2}_\\d{2}$", names(data), value = TRUE) # Identify tmin/tmax daily columns and parse their dates
tmax_cols  <- grep("^tmax_\\d{4}_\\d{2}_\\d{2}$", names(data), value = TRUE)

tmin_dates <- as.Date(sub("^tmin_", "", tmin_cols), format = "%Y_%m_%d")
tmax_dates <- as.Date(sub("^tmax_", "", tmax_cols), format = "%Y_%m_%d")

common_dates <- sort(as.Date(intersect(tmin_dates, tmax_dates))) # Keep only days with both Tmin and Tmax, and be sure they're sorted Dates


if (length(common_dates) == 0L) { # If no overlap, bail early
  data$GDD_planting_window <- NA_real_
} else {
  
  # Rebuild aligned column name vectors in the same order
  tmin_common_cols <- paste0("tmin_", format(common_dates, "%Y_%m_%d"))
  tmax_common_cols <- paste0("tmax_", format(common_dates, "%Y_%m_%d"))
  
  # Build a matrix of daily mean temperatures (force numeric)
  tmin_mat <- as.matrix(data[tmin_common_cols]); storage.mode(tmin_mat) <- "double"
  tmax_mat <- as.matrix(data[tmax_common_cols]); storage.mode(tmax_mat) <- "double"
  tmean_mat <- (tmin_mat + tmax_mat) / 2
  
  # Pre-compute indices for each (year, month) planting window
  data$experiment_year <- as.integer(data$experiment_year)
  years  <- sort(unique(data$experiment_year))
  months <- sprintf("%02d", 1:12)
  
  idx_by_yrmon <- vector("list", length = length(years) * 12L)
  names(idx_by_yrmon) <- as.vector(outer(years, months, paste, sep = "_"))
  
  for (yr in years) {
    for (m in months) {
      start <- as.Date(sprintf("%04d-%s-01", yr, m))
      end   <- as.Date(sprintf("%04d-%s-01", yr + 1L, m))  # exclusive
      idx   <- which(common_dates >= start & common_dates < end)
      idx_by_yrmon[[paste0(yr, "_", m)]] <- idx
    }
  }
  
  # Calculate GDD (thermal time) per row over its planting window
  
  Tbase <- 10 # base temperature for cassava (Paredes et al. 2025)
  
  keys     <- paste0(data$experiment_year, "_", sprintf("%02d", data$planting_month))
  idx_list <- idx_by_yrmon[keys]
  
  data$GDD_planting_window <- mapply(
    function(i, idx) {
      if (length(idx) == 0L) return(NA_real_)
      daily_gdd <- pmax(0, tmean_mat[i, idx] - Tbase)
      sum(daily_gdd, na.rm = TRUE)
    },
    seq_len(nrow(data)), idx_list
  )
}



# b)Rainfall aggregation-----------------------------------------------------

# The rainfall is determined as the cumulative rainfall for 12 months from planting date

## 1) Find rainfall columns of the form pmm_<mon>_<year>

rain_cols <- grep("^(?:pmm|ppm)_[a-z]{3}_\\d{4}$", names(data), value = TRUE, ignore.case = TRUE)

if (length(rain_cols) == 0L) {
  data$pmm_sum_12mo <- NA_real_
} else {
  ## 2) Parse month + year from column names and order them chronologically
  rx <- "^(?:pmm|ppm)_([a-z]{3})_(\\d{4})$"
  mm_yy <- do.call(rbind, regmatches(rain_cols, regexec(rx, rain_cols, ignore.case = TRUE)))
  # mm_yy columns: full match, month, year
  mon_chr <- tolower(mm_yy[, 2])
  yr_int  <- as.integer(mm_yy[, 3])
  
  mon_lut <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")
  mon_num <- match(mon_chr, mon_lut)
  
  rain_dates <- as.Date(sprintf("%04d-%02d-01", yr_int, mon_num))  # first of each month
  
  ord <- order(rain_dates)
  rain_dates <- rain_dates[ord]
  rain_cols  <- rain_cols[ord]
  
  ## 3) Build a matrix of monthly rainfall (force numeric)
  pmm_mat <- as.matrix(data[rain_cols])
  storage.mode(pmm_mat) <- "double"
  
  ## 4) Pre-compute indices for each (year, month) 12-month window
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
  
  ## 5) Sum rainfall over the 12-month window for each row
  keys <- paste0(
    data$experiment_year, "_",
    sprintf("%02d", as.integer(data$planting_month))  # robust if planting_month is character
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


# c) SPEI INDEX---------------------------------------------------------------

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



# d) Aggregation of Soil texture, SOC, BD, water content and pH by soil depth-------

# These soil properties are given for depths 0–5, 5–15, 15–30, 30–60, 60–100, 100-200.
# However, we assume that cassava roots go up-to 100m. A weighted average for each property 
# is determined and related to the yield.

weights_0_100 <- c(5, 10, 15, 30, 40)   # 0–5, 5–15, 15–30, 30–60, 60–100
total_depth <- sum(weights_0_100)       # = 100

data <- data %>%
  mutate(
    clay = {
      cols <- select(., matches("^clay_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    sand = {
      cols <- select(., matches("^sand_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    silt = {
      cols <- select(., matches("^silt_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    SOC = {
      cols <- select(., matches("^socgkg_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    BD = {
      cols <- select(., matches("^bd_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_10kpa = {
      cols <- select(., matches("^vwc_10kpa_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_33kpa = {
      cols <- select(., matches("^vwc_33kpa_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_1500kpa = {
      cols <- select(., matches("^vwc_1500kpa_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    soil_pH = {
      cols <- select(., matches("^soilgrids_ph_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      (rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth)
    }
  )


##---- 3. FITTING BOUNDARY LINE MODELS -----------------------------------------

#-- 3.1. Growing Degree Days ===========================================

dat_GDD <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","GDD_planting_window")]
dat_GDD <- dat_GDD[!duplicated(dat_GDD), ] # remove duplicate entries

x <- dat_GDD$GDD_planting_window
y <- dat_GDD$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------
dat_GDD <- data.frame(x=x, y=sqrt(y))
out <- bagplot(dat_GDD, show.whiskers=F, factor = 3.0)
dat_GDD <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------

plot(dat_GDD[,1], dat_GDD[,2], pch=16, col="grey")
startValues("lp")# determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(20.34, -0.002, 7.97, mean(dat_GDD[,1]), mean(dat_GDD[,2]), sd(dat_GDD[,1]), sd(dat_GDD[,2]), cor(dat_GDD[,1],dat_GDD[,2])),
               c(18.62, -0.002, 8.09, mean(dat_GDD[,1]), mean(dat_GDD[,2]), sd(dat_GDD[,1]), sd(dat_GDD[,2]), cor(dat_GDD[,1],dat_GDD[,2])),
               c(19.62, -0.002, 7.82, mean(dat_GDD[,1]), mean(dat_GDD[,2]), sd(dat_GDD[,1]), sd(dat_GDD[,2]), cor(dat_GDD[,1],dat_GDD[,2])),
               c(18.53, -0.002, 7.74, mean(dat_GDD[,1]), mean(dat_GDD[,2]), sd(dat_GDD[,1]), sd(dat_GDD[,2]), cor(dat_GDD[,1],dat_GDD[,2])),
               c(22.01, -0.003, 8.08, mean(dat_GDD[,1]), mean(dat_GDD[,2]), sd(dat_GDD[,1]), sd(dat_GDD[,2]), cor(dat_GDD[,1],dat_GDD[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error

ble_extention(data = dat_GDD, start = start2, sigh = sigh, model = "lp") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_GDD <- cbvn_extention(data=dat_GDD, start = start2, sigh = 0.4, model = "lp", pch=16, col="grey",
                            ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                            xlab = expression(bold("GDD / "^o*"C")))

akweight(model_GDD$AIC[1,1], model_GDD$AIC[2,1])# check the akaike weights for boundary and null model

x4 <- data$GDD_planting_window
x4[which(is.na(x4)==T)] <- mean(x4, na.rm = T)
GDD <- predictBL(model_GDD, x4)# predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

x_GDD <- seq(4000,6500,1)
y_GDD <- predictBL(model_GDD, x_GDD)
plot(dat_GDD[,1], dat_GDD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("GDD / "^o*"C")))
lines(x_GDD, y_GDD^2, col="red", lwd=2)  


#-- 3.2. Precipitation (Rainfall) ===============================================

dat_pmm <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","pmm_yr_data")]
dat_pmm <- dat_pmm[!duplicated(dat_pmm), ]# remove duplicate entries

x <- dat_pmm$pmm_yr_data
y <- dat_pmm$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_pmm <- data.frame(x=x, y=sqrt(y))
out <- bagplot(dat_pmm, show.whiskers=F, factor = 3.0, na.rm = TRUE)
dat_pmm <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values---------------------------------------------------

plot(dat_pmm[,1], dat_pmm[,2])
startValues("lp")# determine start values by clicking on the plot the points that make up the desired model.
# In this case for the linear-plateau, 2 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(7.72, 0.00, 7.72, mean(dat_pmm[,1]), mean(dat_pmm[,2]), sd(dat_pmm[,1]), sd(dat_pmm[,2]), cor(dat_pmm[,1],dat_pmm[,2])),
               c(7.82, 0.00, 7.82, mean(dat_pmm[,1]), mean(dat_pmm[,2]), sd(dat_pmm[,1]), sd(dat_pmm[,2]), cor(dat_pmm[,1],dat_pmm[,2])),
               c(8.00, 0.00, 8.00, mean(dat_pmm[,1]), mean(dat_pmm[,2]), sd(dat_pmm[,1]), sd(dat_pmm[,2]), cor(dat_pmm[,1],dat_pmm[,2])),
               c(7.42, 0.00, 7.42, mean(dat_pmm[,1]), mean(dat_pmm[,2]), sd(dat_pmm[,1]), sd(dat_pmm[,2]), cor(dat_pmm[,1],dat_pmm[,2])))


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
ble_extention(data = dat_pmm, start = start2, sigh = sigh, model = "lp")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_Precip_yr <- cbvn_extention(data=dat_pmm, start = start2, sigh = 0.2, model = "lp", pch=16, col="grey",
                                  ylab=expression(bold("Yield / "* sqrt(t~ha^-1))), 
                                  xlab=expression(bold("Rainfall/ mm yr"^-1)))

akweight(model_Precip_yr$AIC[1,1], model_Precip_yr$AIC[2,1])# check the akaike weights for boundary and null model

x3 <- data$pmm_yr_data
x3[which(is.na(x3)==T)] <- mean(x3, na.rm = T)
rainfall <- predictBL(model_Precip_yr, x3)# predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

x_pmm <- seq(550,2000,1)
y_pmm <- predictBL(model_Precip_yr, x_pmm)
plot(dat_pmm[,1], dat_pmm[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("Rainfall/ mm yr"^-1)))
lines(x_pmm, y_pmm^2, col="red", lwd=2)  


#-- 3.3 SPEI Index -------------------------------------------------------------

#- a) Summary statistics--------------------------------------------------------

dat_spei <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","spei_yr_data")]
dat_spei <- dat_spei[!duplicated(dat_spei), ] # remove duplicate entries

x <- dat_spei$spei_yr_data
y <- dat_spei$cassava_fr_root_yld_tha

summastat(x)
summastat(log(x+1)) # check for normality and transform if necessary
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------
dat_spei <- data.frame(x=log(x+1), y=sqrt(y))
out <- bagplot(dat_spei, show.whiskers=F, na.rm = TRUE)
dat_spei <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------

plot(dat_spei[,1], dat_spei[,2])
startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the lp, 2 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(9.06, 12.78, 7.77,15.39,-22.91, mean(dat_spei[,1]), mean(dat_spei[,2]), sd(dat_spei[,1]), sd(dat_spei[,2]), cor(dat_spei[,1],dat_spei[,2])),
               c(8.89, 12.68, 7.64,10.28,-7.79, mean(dat_spei[,1]), mean(dat_spei[,2]), sd(dat_spei[,1]), sd(dat_spei[,2]), cor(dat_spei[,1],dat_spei[,2])),
               c(11.12, 19.28, 7.68,17.54,-26.82, mean(dat_spei[,1]), mean(dat_spei[,2]), sd(dat_spei[,1]), sd(dat_spei[,2]), cor(dat_spei[,1],dat_spei[,2])),
               c(9.15, 13.00, 7.80,11.31,-10.55, mean(dat_spei[,1]), mean(dat_spei[,2]), sd(dat_spei[,1]), sd(dat_spei[,2]), cor(dat_spei[,1],dat_spei[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
ble_extention(data = dat_spei, start = start2, sigh = sigh, model = "trapezium") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) model fitting --------------------------------------------------------------

par(mar=c(6,5,4,2)) # adjusting ploting parameters

model_spei <- cbvn_extention(data=dat_spei, start = start2, sigh = 0.7, model = "trapezium", pch=16, col="grey",
                             ylab=expression(bold("Yield / "* sqrt(t~ha^-1))), 
                             xlab=expression(bold("SPEI INDEX")))

akweight(model_spei$AIC[1,1], model_spei$AIC[2,1]) # check the akaike weights for boundary and null model

x1 <- log(data$spei_yr_data+1)
x1[which(is.na(x1)==T)] <- mean(x1, na.rm = T)
spei <- predictBL(model_spei, x1) # predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------

x_spei<- seq(-0.35,0.6, 0.001)
y_spei <- predictBL(model_spei, x_spei)
plot(exp(dat_spei[,1])-1, dat_spei[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("SPEI INDEX")))
lines(exp(x_spei)-1, y_spei^2, col="red", lwd=2)  


#-- 3.4. Solar Radiation model =================================================

dat_SR <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","SR_annual_MJ_m2_tilted")]
dat_SR <- dat_SR[!duplicated(dat_SR), ]# remove duplicate entries

x <- dat_SR$SR_annual_MJ_m2_tilted
y <- dat_SR$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x)# check for normality and transform if necessary
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_SR <- data.frame(x=x, y=sqrt(y))
out <- bagplot(dat_SR, show.whiskers=F, factor = 3.0, na.rm=TRUE)
dat_SR <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------

plot(dat_SR[,1], dat_SR[,2])
startValues("trapezium")# determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(7.90, 0.00, 7.91, mean(dat_SR[,1]), mean(dat_SR[,2]), sd(dat_SR[,1]), sd(dat_SR[,2]), cor(dat_SR[,1],dat_SR[,2])),
               c(7.91, 0.00, 7.94, mean(dat_SR[,1]), mean(dat_SR[,2]), sd(dat_SR[,1]), sd(dat_SR[,2]), cor(dat_SR[,1],dat_SR[,2])),
               c(7.93, 0.00, 7.93, mean(dat_SR[,1]), mean(dat_SR[,2]), sd(dat_SR[,1]), sd(dat_SR[,2]), cor(dat_SR[,1],dat_SR[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error

ble_extention(data = dat_SR, start = start2, sigh = sigh, model = "lp") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_SR <- cbvn_extention(data=dat_SR, start = start2, sigh = 0.4, model = "lp", pch=16, col="grey",
                           ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                           xlab = expression(bold("Annual Solar Radiation / " * log(MJ~m^-2~year^-1))))

akweight(model_SR$AIC[1,1], model_SR$AIC[2,1])# check the akaike weights for boundary and null model

x4 <- log(data$SR_annual_MJ_m2_tilted)
x4[which(is.na(x4)==T)] <- mean(x4, na.rm = T)
SR <- predictBL(model_SR, x4)# predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

x_SR <- seq(5960,7400,1)
y_SR <- predictBL(model_SR, x_SR)
plot(dat_SR[,1], dat_SR[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("Annual Solar Radiation / MJm "^-2*" year"^-1)))
lines(x_SR, y_SR^2, col="red", lwd=2)  

#An annual total of 5,000–7,000 MJ m⁻² year⁻¹ is generally favorable cassava

#-- 3.5. Soil clay model =========================================================

dat_clay <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","clay")]
dat_clay <- dat_clay[!duplicated(dat_clay), ] # removes duplicated rows
dat_clay <- dat_clay[-which(dat_clay$clay==0),]

x <- dat_clay$clay
y <- dat_clay$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_clay <- data.frame(x=x, y=sqrt(y))
out <- bagplot(dat_clay, show.whiskers=F, factor = 3.0)
dat_clay <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------

plot(dat_clay[,1], dat_clay[,2])

startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(-5.27, 0.61, 7.91, 32.66, -0.58, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-7.07, 0.70, 7.68, 33.94, -0.58, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-7.10, 0.72, 7.74, 38.70, -0.67, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-5.92, 0.63, 7.67, 50.35, -0.91,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-5.32, 0.59, 8.04, 32.82, -0.57,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-5.12, 0.59, 7.70, 35.62, -0.63,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])))


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
# why only up-to 0.8 based on variogram result for upper limit

ble_extention(data = dat_clay, start = start2, sigh = sigh,model = "trapezium") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------


model_clay <- cbvn_extention(data=dat_clay, start = start2, sigh = 0.3, model = "trapezium", pch=16, col="grey",
                             ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                             xlab = expression(bold("Clay content / %")))

akweight(model_clay$AIC[1,1], model_clay$AIC[2,1]) # check the akaike weights for boundary and null model

x5 <- data$clay
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
clay <- predictBL(model_clay, x5) # predicts the boundary value for each value of x in the data set
clay <- ifelse(clay< 0,0, clay)

#- e) Plot boundary line on original scale--------------------------------------

x_clay<- seq(14,50,0.1)
y_clay <- predictBL(model_clay, x_clay)
plot(dat_clay[,1], dat_clay[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("Clay content / %")))
lines(x_clay, y_clay^2, col="red", lwd=2)  


#-- 3.6. Sand model ============================================================

dat_sand <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","sand")]
dat_sand <- dat_sand[!duplicated(dat_sand), ] # remove duplicate entries
dat_sand <- dat_sand[-which(dat_sand$sand==0), ]

x <- dat_sand$sand
y <- dat_sand$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_sand <- data.frame(x=x, y=sqrt(y))
out <- bagplot(dat_sand, show.whiskers=F, factor = 3.0)
dat_sand <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error----------------

plot(dat_sand[,1], dat_sand[,2])
startValues("lp")# determine start values by clicking on the plot the points that make up the desired model.


start2 <- list(c(2.00, 0.13, 7.99, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(2.25, 0.12, 7.95, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(2.29, 0.12, 8.04, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(2.97, 0.11, 7.99, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(1.41, 0.19, 7.73, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(1.09, 0.16, 7.83, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(2.35, 0.11, 7.88, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
               c(1.78, 0.13, 7.79, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
ble_extention(data = dat_sand, start = start2, sigh = sigh, model = "lp" ) # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_sand <- cbvn_extention(data=dat_sand, start = start2, sigh = 0.2, model = "lp", pch=16, col="grey",
                             ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                             xlab = expression(bold("sand content / %")))

akweight(model_sand$AIC[1,1], model_sand$AIC[2,1])# check the akaike weights for boundary and null model

x5 <- data$sand
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
sand <- predictBL(model_sand, x5) # predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

x_sand<- seq(35,80,0.1)
y_sand <- predictBL(model_sand, x_sand)
plot(dat_sand[,1], dat_sand[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("sand content / %")))
lines(x_sand, y_sand^2, col="red", lwd=2)  


#-- 3.7. Soil pH model ===========================================================

data_pH <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","soil_pH")]
data_pH <- data_pH[!duplicated(data_pH), ]# remove duplicate entries
data_pH <- data_pH[-which(data_pH$soil_pH==0),]

x <- data_pH$soil_pH
y <- data_pH$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_pH <- data.frame(x=x, y=sqrt(y))
out <- bagplot(data_pH, show.whiskers=F, factor = 3.0)
data_pH <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------------

plot(data_pH[,1], data_pH[,2])
startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(-24.96, 5.99, 7.92, 30.84, -3.54, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-27.89, 6.77, 7.77, 33.28, -3.92, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-21.48, 5.30, 7.68, 28.72, -3.26, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-22.13, 5.45, 7.60, 35.28, -4.25, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6)# vector of possible values for sd of measurement error
ble_extention(data = data_pH, start = start2, sigh = sigh, model = "trapezium") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_pH <- cbvn_extention(data=data_pH, start = start2, sigh = 0.4, model = "trapezium", pch=16, col="grey",
                           ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                           xlab = expression(bold("Soil pH")))

akweight(model_pH$AIC[1,1], model_pH$AIC[2,1]) # check the akaike weights for boundary and null model

x5 <- data$soil_pH
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
pH <- predictBL(model_pH , x5) # predicts the boundary value for each value of x in the data set
pH <- ifelse(pH< 0,0, pH)

#- e) Plot boundary line on original scale---------------------------------------

x_pH<- seq(4.7,7.3,0.01)
y_pH <- predictBL(model_pH, x_pH)
plot(data_pH[,1], data_pH[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("soil pH")))
lines(x_pH, y_pH^2, col="red", lwd=2)  

#-- 3.8. Soil organic carbon model =============================================

data_SOC <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","SOC")]
data_SOC <- data_SOC[!duplicated(data_SOC), ]# remove duplicate entries
data_SOC <- data_SOC[-which(data_SOC$SOC==0),]

x <- data_SOC$SOC
y <- data_SOC$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x)# check for normality and transform if necessary
summastat(log(x))# adds a very small number to x to allow for log transformation
summastat(sqrt(y))
plot(log(x),sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_SOC <- data.frame(x=log(x), y=sqrt(y))
out <- bagplot(data_SOC, show.whiskers=F, factor = 3.0, na.rm = T)
data_SOC <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values---------------------------------------------------

plot(data_SOC[,1], data_SOC[,2])
startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 2 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(-17.66, 12.45, 7.64, 25.88, -6.27, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
               c(-8.20, 7.74, 7.55, 28.62, -7.30, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
               c(-9.19, 8.33, 7.99, 25.83, -6.48, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
               c(-19.97, 14.31, 7.89, 29.75, -7.68, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
               c(-18.73, 13.64, 7.78, 28.30, -7.28, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
ble_extention(data = data_SOC, start = start2, sigh = sigh, model = "trapezium")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_SOC <- cbvn_extention(data=data_SOC, start = start2, sigh = 0.4, model = "trapezium", pch=16, col="grey",
                            ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                            xlab = expression(bold("SOC")))

akweight(model_SOC$AIC[1,1], model_SOC$AIC[2,1]) # check the akaike weights for boundary and null model

x5 <- log(data$SOC+0.000001)
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
SOC <- predictBL(model_SOC , x5) # predicts the boundary value for each value of x in the data set
SOC <- ifelse(SOC< 0,0, SOC )

#- e) Plot boundary line on original scale--------------------------------------

x_SOC<- seq(1.6,3.5,0.01)
y_SOC <- predictBL(model_SOC, x_SOC)
plot(exp(data_SOC[,1]), data_SOC[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("SOC")))
lines(exp(x_SOC), y_SOC^2, col="red", lwd=2)

#-- 3.9. Soil bulk density model ===============================================

data_BD <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","BD")]
data_BD <- data_BD[!duplicated(data_BD), ]# remove duplicate entries
data_BD <- data_BD[-which(data_BD$BD==0),]

x <- data_BD$BD/1000 # converts bulk density from kg/m3 to g/cm3
y <- data_BD$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_BD <- data.frame(x=x, y=sqrt(y))
#data_BD <- data_BD[-which(data_BD$x==0),]
out <- bagplot(data_BD, show.whiskers=F, factor = 3.0)
data_BD <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values---------------------------------------------------

plot(data_BD[,1], data_BD[,2])
startValues("trapezium")# determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(-73.96, 62.06, 7.95, 99.99, -67.02, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])),
               c(-102.47, 84.37, 7.75, 148.14, -100.91, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])),
               c(-91.79, 76.07, 7.87, 161.32, -110.32, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])),
               c(-106.01, 86.91, 7.87, 149.87, -102.10, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
ble_extention(data = data_BD, start = start2, sigh = sigh, model = "trapezium")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_BD <- cbvn_extention(data=data_BD, start = start2, sigh = 0.3, model = "trapezium", pch=16, col="grey",
                           ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                           xlab = expression(bold("BD g/kg")))

akweight(model_BD$AIC[1,1], model_BD$AIC[2,1]) # check the akaike weights for boundary and null model

x5 <- data$BD/1000
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
BD <- predictBL(model_BD , x5) # predicts the boundary value for each value of x in the data set
BD <- ifelse(BD< 0,0, BD )

#- e) Plot boundary line on original scale--------------------------------------

x_BD<- seq(1.25,1.45,0.001)
y_BD <- predictBL(model_BD, x_BD)
plot(data_BD[,1], data_BD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("BD g/cm"^3)))
lines(x_BD, y_BD^2, col="red", lwd=2)  


#-- 3.10. Water_holding capacity model==========================================

data_PAW <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","vwc_33kpa","vwc_1500kpa")]
data_PAW <- data_PAW[!duplicated(data_PAW), ] # remove duplicate entries
data_PAW$PAW <- data_PAW$vwc_33kpa - data_PAW$vwc_1500kpa # Calculates plant available water (PAW)
data_PAW <- data_PAW[-which(data_PAW$PAW <= 0),]

x <- data_PAW$PAW/1000 # converts PAW to cm3/cm3
y <- data_PAW$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x)  # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_PAW <- data.frame(x=x, y=sqrt(y))
out <- bagplot(data_PAW, show.whiskers=F, factor = 3.0, na.rm = T)
data_PAW <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error-----------------------

plot(data_PAW[,1], data_PAW[,2])
startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(1.25, 64.93, 7.66, 28.63, -103.61, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               c(-1.06, 101.23, 7.89, 35.20, -132.79, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               c(-0.98, 100.73, 7.59, 39.10, -150.39, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               c(-1.36, 105.42, 7.68, 33.97, -127.24, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               c(-0.34, 83.93, 7.77, 34.76,-130.93, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               c(-1.02, 97.83, 7.55, 26.44, -94.20, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
ble_extention(data = data_PAW, start = start2, sigh = sigh, model = "trapezium")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_paw <- cbvn_extention(data=data_PAW, start = start2, sigh = 0.2, model = "trapezium", pch=16, col="grey",
                            ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                            xlab = expression(bold("PAW cm"^-3*"cm"^-3)))

akweight(model_paw$AIC[1,1], model_paw$AIC[2,1]) # check the akaike weights for boundary and null model

x5 <- (data$vwc_33kpa - data$vwc_1500kpa)/1000
x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
paw <- predictBL(model_paw , x5) # predicts the boundary value for each value of x in the data set
paw <- ifelse(paw< 0,0, paw )

#- e) Plot boundary line on original scale--------------------------------------

x_paw<- seq(0.06,0.23,0.001)
y_paw <- predictBL(model_paw, x_paw)
plot(data_PAW[,1], data_PAW[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("PAW cm"^-3*"cm"^-3)))
lines(x_paw, y_paw^2, col="red", lwd=2)  


#-- 3.11. Disease ==============================================================

disease_data <- data.frame(x=as.factor(data$disease),y=data$cassava_fr_root_yld_tha)

# keeps only rows where x and y are non-missing and non-empty.
disease_data <- disease_data[!(is.na(disease_data$x) | is.na(disease_data$y) |
                                 disease_data$x == "" | disease_data$y == ""), ]

#- a) Outlier detection---------------------------------------------------------

disease_data <- do.call(rbind, lapply(split(disease_data, f=disease_data$x), outlier_remove, y=y, d=3))
rownames(disease_data) <- NULL

#- b) check for distribution of yield ~factor ----------------------------------

lapply(split(disease_data$y, f=disease_data$x), summastat,plot = TRUE)
disease_data <- data.frame(x=disease_data$x, y=sqrt(disease_data$y)) #transform yield 
lapply(split(disease_data$y, f=disease_data$x), summastat,plot = TRUE)

#- c) Initial starting values---------------------------------------------------

sigh<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
ble_Cat(x=disease_data$x,y=disease_data$y,sigh=sigh) # Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_disease <- ccnm(x=disease_data$x,y=disease_data$y,sigh=0.1, plot = TRUE)
akweight(model_disease $AIC[1,1], model_disease $AIC[2,1]) # check the akaike weights for boundary and null model

disease <- predictBL2(model_disease, data$disease) # predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

stripchart(disease_data$y^2 ~ disease_data$x, method = "jitter", jitter = 0.1, pch = 16, 
           col = "black", cex = 0.6, vertical = TRUE, 
           ylab = "", xlab = "Disease")
x_pos <- 1:length(model_disease[[2]][,2])
y_val <- model_disease[[2]][,2]^2
points(x_pos, y_val, pch = 16, col = "red") # add censoring points
segments(x_pos - 0.1, y_val, x_pos + 0.1, y_val, col = "red", lwd = 2)# Small horizontal lines through each point


#-- 3.12. Disease Severity =====================================================

severity_data <- data.frame(x=as.factor(data$severity),y=data$cassava_fr_root_yld_tha)

# keeps only rows where x and y are non-missing and non-empty.

severity_data <- severity_data[!(is.na(severity_data$x) | is.na(severity_data$y) |
                                   severity_data$x == "" | severity_data$y == ""), ]

#- a) Outlier detection---------------------------------------------------------

severity_data <- do.call(rbind, lapply(split(severity_data, f=severity_data$x), outlier_remove, y=y, d=3))
rownames(severity_data) <- NULL

#- b) check for distribution of yield ~factor ----------------------------------

lapply(split(severity_data$y, f=severity_data$x), summastat,plot = TRUE)
severity_data <- data.frame(x=severity_data$x, y=sqrt(severity_data$y)) #transform yield 
lapply(split(severity_data$y, f=severity_data$x), summastat,plot = TRUE)

#- c) Initial starting values---------------------------------------------------

sigh<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
ble_Cat(x=severity_data$x,y=severity_data$y,sigh=sigh) # Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_severity <- ccnm(x=severity_data$x,y=severity_data$y,sigh=0.1, plot = TRUE)
akweight(model_severity $AIC[1,1], model_severity $AIC[2,1]) # check the akaike weights for boundary and null model

disease_severity <- predictBL2(model_severity, data$severity) # predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

stripchart(severity_data$y^2 ~ severity_data$x, method = "jitter", jitter = 0.1, pch = 16, 
           col = "black", cex = 0.6, vertical = TRUE, 
           ylab = "", xlab = "Disease severity")
x_pos <- 1:length(model_severity[[2]][,2])
y_val <- model_severity[[2]][,2]^2
points(x_pos, y_val, pch = 16, col = "red") # add censoring points
segments(x_pos - 0.1, y_val, x_pos + 0.1, y_val, col = "red", lwd = 2)# Small horizontal lines through each point



##  4. Attainable Yield and the Most Limiting factor============================

# This is achieved by using the limfactor() function of the BLA R package. The inputs
# are the predicted boundary values from the various boundary line models, tmin, tmax, 
# rainfall, sand, pH, clay, SR, SOC, paw, BD determined earlier.

Limiting <- limfactor(sand, pH, clay, SR, SOC, paw, BD, GDD, disease_severity, spei)

# a) Attainable yield-----------------------------------------------------------

# The attainable yield is determined as the largest predicted yield from the 
# various fitted boundary line models. It can be extracted from the object of 
# limfactor() function.

att_yield <- Limiting[[2]]^2
paste("The attainable yield is ", round(att_yield,2)," t/ha")

# b) Most limiting factor-------------------------------------------------------
# The most limiting factor is determined according to the law of the minimum 
# by von Liebig using the function limfactor() above. so you can extract the 
# most limiting factor at each point in the data set. The proportions of the
# identified most limiting factors can be evaluated using the code below.

tab <- table(Limiting[[1]]$Lim_factor)
prop <- prop.table(tab)
prop_sorted <- sort(prop, decreasing = FALSE)
colfunc <- colorRampPalette(c("green", "yellow", "red", "darkred"))

par(mar=c(5,8,4,2))
barplot(
  prop_sorted,
  horiz = TRUE,   # flip axes
  col = colfunc(length(prop_sorted)),  space = 0.6, width = 0.2,
  xlab = expression(bold("Proportion of limiting factor")), ylab = "", xlim = c(0,1), las = 1)
box()



# View all the boundary line plots-----------------------------------------------
{
  par(mfrow=c(3,2))
  par(mar=c(5,7,4,2))
  
  plot(dat_GDD[,1], dat_GDD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("GDD / "^o*"C")))
  lines(x_GDD, y_GDD^2, col="red", lwd=2)  
  
  plot(dat_pmm[,1], dat_pmm[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("Rainfall/ mm yr"^-1)))
  lines(x_pmm, y_pmm^2, col="red", lwd=2)  
  
  plot(dat_spei[,1], dat_spei[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("SPEI INDEX")))
  lines(x_spei, y_spei^2, col="red", lwd=2)  
  
  
  plot(dat_SR[,1], dat_SR[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("Annual Solar Radiation / MJm "^-2*" year"^-1)))
  lines(x_SR, y_SR^2, col="red", lwd=2)  
  
  
  plot(dat_clay[,1], dat_clay[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("Clay content / %")))
  lines(x_clay, y_clay^2, col="red", lwd=2)  
  
  
  plot(dat_sand[,1], dat_sand[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("sand content / %")))
  lines(x_sand, y_sand^2, col="red", lwd=2)  
  
  plot(data_pH[,1], data_pH[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("soil pH")))
  lines(x_pH, y_pH^2, col="red", lwd=2)
  
  plot(data_SOC[,1], data_SOC[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("SOC")))
  lines(x_SOC, y_SOC^2, col="red", lwd=2) 
  
  plot(data_BD[,1], data_BD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("Soil Bulk Density/ g cm"^3)))
  lines(x_BD, y_BD^2, col="red", lwd=2) 
  
  plot(data_PAW[,1], data_PAW[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
       xlab=expression(bold("PAW cm"^-3*"cm"^-3)))
  lines(x_paw, y_paw^2, col="red", lwd=2)
  
  stripchart(disease_data$y^2 ~ disease_data$x, method = "jitter", jitter = 0.1, pch = 16, 
             col = "black", cex = 0.6, vertical = TRUE, 
             ylab = "", xlab = "Disease")
  x_pos <- 1:length(model_disease[[2]][,2])
  y_val <- model_disease[[2]][,2]^2
  points(x_pos, y_val, pch = 16, col = "red") # add censoring points
  segments(x_pos - 0.1, y_val, x_pos + 0.1, y_val, col = "red", lwd = 2)# Small horizontal lines through each point
  
  stripchart(severity_data$y^2 ~ severity_data$x, method = "jitter", jitter = 0.1, pch = 16, 
             col = "black", cex = 0.6, vertical = TRUE, 
             ylab = "", xlab = "Disease severity")
  x_pos <- 1:length(model_severity[[2]][,2])
  y_val <- model_severity[[2]][,2]^2
  points(x_pos, y_val, pch = 16, col = "red") # add censoring points
  segments(x_pos - 0.1, y_val, x_pos + 0.1, y_val, col = "red", lwd = 2)# Small horizontal lines through each point
  
  
  barplot(
    prop_sorted,
    horiz = TRUE,   # flip axes
    col = colfunc(length(prop_sorted)),  space = 0.6, width = 0.2,
    xlab = expression(bold("Proportion of limiting factor")), ylab = "", xlim = c(0,1), las = 1, cex.names=0.8)
  box()
  
  par(mfrow=c(1,1))
}


#--- 5. Mapping attainable yields and the most-limiting factors for Tanzania =========

## 5.1) Yield prediction by each factor of interest for Tanzania ----------------------------

bbox <- c(xmin = 24.90, xmax = 40.56, ymin = -15.38, ymax = 0.5)
depths <- c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm")
w <- c(5, 10, 15, 30, 40)# depth weights in cm (0–5, 5–15, 15–30, 30–60, 60–100)
tz <- vect("Tanzania.shp") # Tanzania shape file for extents
tz <- project(tz, "EPSG:4326") #optional, convert to decimal degrees

## a) Soil Clay-----------------------------------------------------------------

# load the clay raster files for each horizon

covs_clay <- sprintf("clay_%s_mean", depths)
base_clay <- "https://maps.isric.org/mapserv?map=/map/clay.map"

clay_rasters <- lapply(covs_clay, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_clay, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/10
})

clay_stack <- rast(clay_rasters)
names(clay_stack) <- depths

clay_pct <- app(clay_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average
clay_pct <- clamp(clay_pct, 0, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
names(clay_pct) <- "clay_pct"

# Function to Predict yields with parameters determined from boundary line model for clay

yield_fun <- function(clay_percent) {
  y1 <- (model_clay[[3]][1,1]+ model_clay[[3]][2,1]*clay_percent)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y2 <- (model_clay[[3]][4,1]+ model_clay[[3]][5,1]*clay_percent)
  y0 <- model_clay[[3]][3,1]
  
  y <- pmin(y0^2, y1^2,y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_clay <- app(clay_pct, yield_fun)  # pixel-wise apply
names(yield_clay) <- "yield_clay"

yield_clay <- clamp(yield_clay, 0, 120, values=TRUE)# clip yields to realistic range

# Mask map with Tanzania extents
tz <- project(tz, crs(clay_pct))# Make sure CRS matches the clay raster
yield_clay <- mask(crop(yield_clay, tz), tz)

# Clay_yield raster

# out_file <- "yield_from_clay_Tanzania.tif"
# writeRaster(
#   yield_clay, filename = out_file, overwrite = TRUE,
#   datatype = "FLT4S",  # 32-bit float
#   gdal = c("COMPRESS=LZW", "TILED=YES")
# )

# Quick look

plot(yield_clay, main="Predicted Yield (t/ha) from Clay (%)")


# b) Soil sand ----------------------------------------------------------------

# load the clay raster files for each horizon

covs_sand <- sprintf("sand_%s_mean", depths)
base_sand <- "https://maps.isric.org/mapserv?map=/map/sand.map"

sand_rasters <- lapply(covs_sand, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_sand, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/10
})

sand_stack <- rast(sand_rasters)
names(sand_stack) <- depths

sand_pct <- app(sand_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average
sand_pct <- clamp(sand_pct, 0, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
names(sand_pct) <- "sand_pct"

# Function to Predict yields with parameters determined from boundary line model for sand

yield_fun_sand <- function(sand_percent) {
  y1 <- (model_sand[[3]][1,1]+ model_sand[[3]][2,1]*sand_percent)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  #y2 <- (model_sand[[3]][4,1]+ model_sand[[3]][5,1]*sand_percent)
  y0 <- model_sand[[3]][3,1]
  
  y <- pmin(y0^2, y1^2)#, y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_sand <- app(sand_pct, yield_fun_sand)  # pixel-wise apply
names(yield_sand ) <- "yield_sand "

# Optional: clip yields to realistic range

yield_sand  <- clamp(yield_sand , 0, 120, values=TRUE)  # example upper cap

# Mask map with Tanzania extents

yield_sand <- mask(crop(yield_sand, tz), tz)

# sand_yield raster

out_file <- "yield_from_sand_Tanzania.tif"
writeRaster(
  yield_sand , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# Quick look

plot(yield_sand, main="Predicted Yield (t/ha) from sand (%)")


## c) Soil pH-------------------------------------------------------------------

covs_pH   <- sprintf("phh2o_%s_mean", depths)
base   <- "https://maps.isric.org/mapserv?map=/map/phh2o.map"

pH_rasters <- lapply(covs_pH, function(cov_id) {    # download each depth as pH (not pH*10)
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url) / 10  # convert from pH*10 to pH
})

pH_stack <- rast(pH_rasters)
names(pH_stack) <- depths

soilpH <- app(pH_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average

soilpH<- clamp(soilpH, 0, 14, values=TRUE)# Ensure plausible range and set nonsense to NA
names(soilpH) <- "pH"

# Function to Predict yields with parameters determined from boundary line model for soil pH

yield_fun_pH <- function(pH_percent) {
  y1 <- (model_pH[[3]][1,1]+ model_pH[[3]][2,1]*pH_percent)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y2 <- (model_pH[[3]][4,1]+ model_pH[[3]][5,1]*pH_percent)
  y0 <- model_pH[[3]][3,1]
  
  y <- pmin(y0^2, y1^2, y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_pH <- app(soilpH, yield_fun_pH)  # pixel-wise apply
names(yield_pH ) <- "yield_pH "

# Optional: clip yields to realistic range

yield_pH  <- clamp(yield_pH , 0, 120, values=TRUE)  # example upper cap

# Mask map with Tanzania extents

yield_pH <- mask(crop(yield_pH, tz), tz)

# soil pH_yield raster

# out_file <- "yield_from_pH_Tanzania.tif"
# writeRaster(
#   yield_pH , filename = out_file, overwrite = TRUE,
#   datatype = "FLT4S",  # 32-bit float
#   gdal = c("COMPRESS=LZW", "TILED=YES")
# )

# Quick look

plot(yield_pH, main="Predicted Yield (t/ha) from soil pH")


## d) Soil Organic Carbon ------------------------------------------------------

# load the clay raster files for each horizon

covs_soc  <- sprintf("soc_%s_mean", depths)
base_soc  <- "https://maps.isric.org/mapserv?map=/map/soc.map"

soc_rasters <- lapply(covs_soc, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_soc, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/10
})

soc_stack <- rast(soc_rasters)
names(soc_stack) <- depths

soc_pct <- app(soc_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average
soc_pct <- clamp(soc_pct, 0, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
names(soc_pct) <- "soc_pct"

# Function to Predict yields with parameters determined from boundary line model for SOC

yield_fun_soc <- function(soc) {
  y1 <- (model_SOC[[3]][1,1]+ model_SOC[[3]][2,1]*log(soc+0.0001))
  y1 <- ifelse( y1 < 0, 0, y1) # enforce non-negative and cap at a biological max if you like
  y2 <- (model_SOC[[3]][4,1]+ model_SOC[[3]][5,1]*log(soc+0.0001))
  y2 <- ifelse( y2 < 0, 0, y2) # enforce non-negative and cap at a biological max if you like
  y0 <- model_SOC[[3]][3,1]
  
  y <- pmin(y0^2, y1^2, y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_soc <- app(soc_pct, yield_fun_soc)  # pixel-wise apply
names(yield_soc ) <- "yield_soc "

# Optional: clip yields to realistic range

yield_soc  <- clamp(yield_soc , 0, 120, values=TRUE)  # example upper cap

yield_soc <- mask(crop(yield_soc, tz), tz)

# soc_yield raster

out_file <- "yield_from_soc_Tanzania.tif"
writeRaster(
  yield_soc , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# Quick look

plot(yield_soc, main="Predicted Yield (t/ha) from soc (%)")


## e) Soil bulk density --------------------------------------------------------

# load the clay raster files for each horizon

covs_bd   <- sprintf("bdod_%s_mean", depths)
base_bd   <- "https://maps.isric.org/mapserv?map=/map/bdod.map"

bd_rasters <- lapply(covs_bd, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_bd, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/100
})

bd_stack <- rast(bd_rasters)
names(bd_stack) <- depths

bd <- app(bd_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average
bd <- clamp(bd, 0, 2.5, values=TRUE)# Ensure plausible range and set nonsense to NA
names(bd) <- "bd"

# Function to Predict yields with parameters determined from boundary line model for BD

yield_fun_bd <- function(bd) {
  y1 <- (model_BD[[3]][1,1]+ model_BD[[3]][2,1]*bd)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y2 <- (model_BD[[3]][4,1]+ model_BD[[3]][5,1]*bd)
  y0 <- model_BD[[3]][3,1]
  
  y <- pmin(y0^2, y1^2, y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_bd <- app(bd, yield_fun_bd)  # pixel-wise apply
names(yield_bd ) <- "yield_bd "

# Optional: clip yields to realistic range

yield_bd  <- clamp(yield_bd , 0, 120, values=TRUE)  # example upper cap

yield_bd <- mask(crop(yield_bd, tz), tz)

# BD_yield raster

out_file <- "yield_from_bd_Tanzania.tif"
writeRaster(
  yield_bd , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# Quick look

plot(yield_bd, main="Predicted Yield (t/ha) from BD")


## f) Available Water ----------------------------------------------------------


# load the clay raster files for each horizon

## i) Field capacity

covs_w33 <- sprintf("wv0033_%s_mean", depths)
base_w33 <- "https://maps.isric.org/mapserv?map=/map/wv0033.map"

w33_rasters <- lapply(covs_w33, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_w33, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/1000
})

w33_stack <- rast(w33_rasters)
names(w33_stack) <- depths

w33 <- app(w33_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average


## ii) Wilting Point

covs_w1500 <- sprintf("wv1500_%s_mean", depths)
base_w1500 <- "https://maps.isric.org/mapserv?map=/map/wv1500.map"

w1500_rasters <- lapply(covs_w1500, function(cov_id) {
  url <- sprintf(paste0(
    "%s&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    "&COVERAGEID=%s",
    "&FORMAT=GEOTIFF_INT16",
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/4326",
    "&SUBSET=Long(%f,%f)&SUBSET=Lat(%f,%f)"),
    base_w1500, cov_id, bbox["xmin"], bbox["xmax"], bbox["ymin"], bbox["ymax"]
  )
  rast(url)/1000
})

w1500_stack <- rast(w1500_rasters)
names(w1500_stack) <- depths

w1500 <- app(w1500_stack, fun = function(x) sum(x * w) / sum(w)) # weighted average

## iii) available Water 

aw <- w33 - w1500
aw <- clamp(aw, 0, 1, values=TRUE)# Ensure plausible range and set nonsense to NA
names(aw) <- "aw"

# Function to Predict yields with parameters determined from boundary line model for PAW

yield_fun_aw <- function(aw) {
  y1 <- (model_paw[[3]][1,1]+ model_paw[[3]][2,1]*aw)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y2 <- (model_paw[[3]][4,1]+ model_paw[[3]][5,1]*aw)
  y0 <- model_paw[[3]][3,1]
  
  y <- pmin(y0^2, y1^2, y2^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_aw <- app(aw, yield_fun_aw)  # pixel-wise apply
names(yield_aw ) <- "yield_aw "

# Optional: clip yields to realistic range

yield_aw  <- clamp(yield_aw , 0, 120, values=TRUE)  # example upper cap

yield_aw <- mask(crop(yield_aw, tz), tz)

# PAW_yield raster

out_file <- "yield_from_aw_Tanzania.tif"
writeRaster(
  yield_aw , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

# Quick look

plot(yield_aw, main="Predicted Yield (t/ha) from aw (%)")



## 5.2) Law of Minimum prediction by factor -------------------------------------

# a) Minimum yield predicted

yield_min <- min(yield_clay, yield_sand, yield_pH, 
                 yield_soc, yield_bd, yield_aw, na.rm = TRUE)
names(yield_min) <- "yield_min_t_ha"


yield_min <- mask(crop(yield_min, tz), tz) # Mask and plot Quick plot

plot(yield_min, main = "Predicted Yield (t/ha) by law of Minimum")

out_file <- "yield_min_Tanzania.tif"
writeRaster(
  yield_min , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)


# b) Mapping the Identified Most limiting factor--------------------------------

x <- c(yield_clay, yield_sand, yield_pH, yield_soc, yield_bd, yield_aw)

names(x) <- sub("^yield_", "", names(x)) # Collects only the last part of the name
labs <- names(x)
n <- nlyr(x)

set.seed(123) # For reproducibility of random tie-breaking

# Function to pick minimum; return NA if min == 0

pick_min_index <- function(v) {
  if (all(is.na(v))) return(NA_real_)
  m <- min(v, na.rm = TRUE)
  if (m == 0) return(NA_real_)   # rule: if min == 0 → NA
  w <- which(v == m)
  if (length(w) == 1) w else sample(w, 1)
}


idx <- app(x, pick_min_index) # Apply function across all pixels

out <- as.factor(idx) # Convert to categorical raster and assign names
levels(out) <- data.frame(ID = seq_len(n), limiting = labs)

plot(out, main="Identified most limiting factor")


{par(mfrow=c(2,1))
  
  plot(yield_min, plg = list(title = "Predicted Yield (t/ha)", cex = 0.8), axes = FALSE,
       mar = c(2,2,2,10))
  
  plot(out, plg = list(title = "Most limiting Factor", cex = 1.3),  # legend title/size
       axes = FALSE, mar = c(2,2,1,10), col=c("#DD605C","#6AE653","#D74FC8","#CAE63E","#73DED1","#929685"))
  
  par(mfrow=c(1,1))}

# c). soil classes-------------------------------------------------------------

soil_class <- rast("soil_classes.tif")
soil_class <- mask(crop(soil_class, tz), tz)
soil_class <- as.factor(soil_class) # make it categorical

# Attach the soil-class names to the raster attribute table (RAT)
soil_names  <- c("Water", "Alisols", "Andosols", "Arenosols", "Calcisols", "Cambisols",
                 "Ferralsols", "Fluvisols", "Gleysols",  
                 "Leptosols", "Lixisols", "Luvisols", "Nitisols", "Phaeozems", 
                 "Planosols", "Plinthosols","Regosols", "Solonchaks", "Solonetz", "Vertisols") 


rat <- levels(soil_class)[[1]]
rat$soil_name <- soil_names     # add a label column
levels(soil_class)[[1]] <- rat[,c(1,3)]

set.seed(123) # for reproducibility                  
cols <- randomcoloR::distinctColorPalette(k = 20)

plot(soil_class, type = "classes", 
     col = cols, plg = list(title = "Soil class", cex = 0.8),  # legend title/size
     axes = FALSE, mar = c(2,2,1,10))

{par(mfrow=c(2,1))
  
  plot(yield_min, plg = list(title = "Predicted Yield (t/ha)", cex = 0.8), axes = FALSE,
       mar = c(2,2,2,10))
  
  plot(soil_class, type = "classes", 
       col = cols, plg = list(title = "Soil class", cex = 0.8),  # legend title/size
       axes = FALSE, mar = c(2,2,1,10))
  
  par(mfrow=c(1,1))}

#-------------------------------------------------------------------------------








