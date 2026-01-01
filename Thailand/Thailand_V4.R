# This script is meant to carryout global yield gap analysis for cassava under 
# the CABI project for Thailand. We employ the boundary line methodology (Webb, 1972) 
# for this purpose.
# It used climatic and soil factors to determine the attainable yield 
# at country level as well as the most limiting factor.


# some useful information on cassava can be found here https://greg.app/cassava-lifecycle/

#------ Install and Load necessary libraries and Source files ------------------------------

library(devtools)
## installs the development version of BLA

library(readr)
library(readxl)
library(dplyr)
library(BLA)
library(aplpack)
library(stringr)
library(terra) # mapping
library(sf) # mapping
library(MASS)
library(ggplot2)
source("Extra_functions.R") 

#---- 1. Data loading into R ===================================================

data <- read_csv("Thailand_soil_grid.csv")

names(data)

## --- 2. Locations ---------------------------------------------------------

Thailand <- st_read("thailand.shp") %>% st_transform(4326)
data$long <- data$decimal_longitude
data$lat <- data$decimal_latitude
pts <- st_as_sf(data, coords = c("long", "lat"), crs = 4326)


ggplot() +
  geom_sf(data = Thailand, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, size = 2, alpha = 0.9) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  xlim(97,106)+
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  )


#--- 2. Data manipulation and cleaning =======================================

# a) Aggregation of Soil texture, SOC, BD, water content and pH by soil depth-------

# These soil properties are given for depths 0–5, 5–15, 15–30, 30–60, 60–100, 100-200.
# However, we assume that cassava roots go up-to 100m. A weighted average for each property 
# is determined and related to the yield.

weights_0_100 <- c(5, 10, 15, 30, 40)   # 0–5, 5–15, 15–30, 30–60, 60–100
total_depth <- sum(weights_0_100)       # = 100

data <- data %>%
  mutate(
    clay = {
      cols <- dplyr::select(., matches("^clay_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    sand = {
      cols <- dplyr::select(., matches("^sand_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    SOC = {
      cols <- dplyr::select(., matches("^SOC_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    BD = {
      cols <- dplyr::select(., matches("^BD_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_33kpa = {
      cols <- dplyr::select(., matches("^VWC30kPa_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_1500kpa = {
      cols <- dplyr::select(., matches("^VWC1500kPa_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    soil_pH = {
      cols <- dplyr::select(., matches("^SoilpH_(0_5|5_15|15_30|30_60|60_100)$"))
      (rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth)
    }
  )


##---- 3. FITTING BOUNDARY LINE MODELS -----------------------------------------

#-- 3.1. Soil clay model =========================================================

dat_clay <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","clay")]
dat_clay <- dat_clay[!duplicated(dat_clay), ] # removes duplicated rows
dat_clay <- dat_clay[-which(dat_clay$clay==0),]

x <- dat_clay$clay/10
y <- dat_clay$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_clay <- data.frame(x=x, y=sqrt(y))
out <- bagplot(data_clay, show.whiskers=F, factor = 3.0, na.rm = TRUE)
data_clay <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

## iii) Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# a) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_clay <- c(mean(data_clay[,1], na.rm = T), mean(data_clay[,2], na.rm = T), sd(data_clay[,1], na.rm = T), 
                  sd(data_clay[,2], na.rm = T), cor(data_clay[,1],data_clay[,2], use = "complete.obs"))

# b) determine the initial model values
plot(data_clay, pch=16, col="grey")
startValues("trapezium",8) # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_clay <- list(c(-3.15, 0.545, 9.43,28.07,-0.45, dist_clay),
                   c(-15.79, 1.08, 9.36,23.49,-0.36, dist_clay),
                   c(-3.32, 0.52, 9.22,26.95,-0.43, dist_clay))

# c) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_clay, start = start_clay, sigh = sigh, model = "trapezium", optim.method = "Nelder-Mead")

## iv) Fit the boundary model---------

model_clay <- cbvn_extention(data = data_clay, start = start_clay, sigh = 0.4, model = "trapezium", 
                                pch=16, col="grey", ylab=expression(bold("Yield / " * sqrt(t~ha^-1))), 
                                xlab=expression(bold("Clay/%")))

akweight(model_clay$AIC[1,1],model_clay$AIC[2,1]) # no clear support for boundary line model. 50% support

#- e) Plot boundary line on original scale--------------------------------------

x_clay<- seq(19,50,0.01)
y_clay <- predictBL(model_clay, x_clay)
plot(data_clay[,1], data_clay[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("Clay/%")))
lines(x_clay, y_clay^2, col="red", lwd=2)


#-- 3.2. Sand model ============================================================

dat_sand <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","sand")]
dat_sand <- dat_sand[!duplicated(dat_sand), ] # remove duplicate entries
data_sand <- dat_sand[-which(dat_sand$sand==0), ]

x <- data_sand$sand/10
y <- data_sand$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(y) # check for normality and transform if necessary
summastat(sqrt(y))

data_sand <- data.frame(x=x, y=sqrt(y))

#- b Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# i) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_sand <- c(mean(data_sand[,1], na.rm = T), mean(data_sand[,2], na.rm = T), sd(data_sand[,1], na.rm = T), 
              sd(data_sand[,2], na.rm = T), cor(data_sand[,1],data_sand[,2], use = "complete.obs"))

# ii) determine the initial model values
plot(data_sand, pch=16, col="grey")
startValues("lp") # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_sand <- list(c(30.09, -0.41, 9.14, dist_sand),
                  c(31.56, -0.44, 9.08, dist_sand),
                  c(32.83, -0.45, 9.32, dist_sand))

# iii) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_sand, start = start_sand, sigh = sigh, model = "lp")

## iv) Fit the boundary model---------

model_sand <- cbvn_extention(data = data_sand, start = start_sand, sigh = 0.3, model = "lp", 
                            pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                            xlab=expression(bold("sand/%")))

akweight(model_sand$AIC[1,1],model_sand$AIC[2,1]) # no clear support for boundary line model. 50% support

#- v) Plot boundary line on original scale--------------------------------------

x_sand<- seq(18,60,0.01)
y_sand <- predictBL(model_sand, x_sand)
plot(data_sand[,1], data_sand[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("sand/%")))
lines(x_sand, y_sand^2, col="red", lwd=2)


#-- 3.3. Soil organic carbon model =============================================

data_SOC <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","SOC")]
data_SOC <- data_SOC[!duplicated(data_SOC), ]# remove duplicate entries
data_SOC <- data_SOC[-which(data_SOC$SOC==0),]

x <- data_SOC$SOC/10
y <- data_SOC$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x)# check for normality and transform if necessary
summastat(log(x))# adds a very small number to x to allow for log transformation
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_SOC <- data.frame(x=log(x), y=sqrt(y))
out <- bagplot(data_SOC, show.whiskers=F, factor = 3.0, na.rm = T)
data_SOC <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

# c Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# a) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_SOC <- c(mean(data_SOC[,1], na.rm = T), mean(data_SOC[,2], na.rm = T), sd(data_SOC[,1], na.rm = T), 
               sd(data_SOC[,2], na.rm = T), cor(data_SOC[,1],data_SOC[,2], use = "complete.obs"))

# b) determine the initial model values
plot(data_SOC, pch=16, col="grey")
startValues("schmidt",8) # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_SOC <- list(c(1, 2.38, 9.67, dist_SOC),
                   c(30, 2.51, 9.32, dist_SOC),
                   c(20, 2.28, 9.74, dist_SOC))

# c) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_SOC, start = start_SOC, sigh = sigh, model = "schmidt", optim.method = "Nelder-Mead")

## iv) Fit the boundary model---------

model_SOC <- cbvn_extention(data = data_SOC, start = start_SOC, sigh = 0.3, model = "schmidt", 
                             pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                             xlab=expression(bold("SOC/%")), optim.method = "Nelder-Mead")

akweight(model_SOC$AIC[1,1],model_SOC$AIC[2,1]) # no clear support for boundary line model. 50% support


#- e) Plot boundary line on original scale--------------------------------------

x_SOC<- seq(1.6,3.3,0.01)
y_SOC <- predictBL(model_SOC, x_SOC)
plot(exp(data_SOC[,1]), data_SOC[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("SOC/ %")))
lines(exp(x_SOC), y_SOC^2, col="red", lwd=2)


#-- 3.4. Soil Bulk density =============================================

data_BD <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","BD")]
data_BD <- data_BD[!duplicated(data_BD), ]# remove duplicate entries
data_BD <- data_BD[-which(data_BD$BD==0),]

x <- data_BD$BD/100 # converts bulk density from kg/m3 to g/cm3
y <- data_BD$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(sqrt(x)) # check for normality and transform if necessary
summastat(log(x)) 
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_BD <- data.frame(x=log(x), y=sqrt(y))
out <- bagplot(data_BD, show.whiskers=F, factor = 3.0, na.rm = TRUE)
data_BD <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

# c Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# a) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_BD <- c(mean(data_BD[,1], na.rm = T), mean(data_BD[,2], na.rm = T), sd(data_BD[,1], na.rm = T), 
              sd(data_BD[,2], na.rm = T), cor(data_BD[,1],data_BD[,2], use = "complete.obs"))

# b) determine the initial model values
plot(data_BD, pch=16, col="grey")
startValues("lp",8) # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_BD <- list(c(25.48, -44.68, 9.13, dist_BD),
                  c(23.81, -41.60, 8.85, dist_BD),
                  c(22.78, -38.20, 8.99, dist_BD))

# c) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_BD, start = start_BD, sigh = sigh, model = "lp")

## iv) Fit the boundary model---------

model_BD <- cbvn_extention(data = data_BD, start = start_BD, sigh = 0.3, model = "lp", 
                            pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                            xlab=expression(bold("BD/ log")))

akweight(model_BD$AIC[1,1],model_BD$AIC[2,1]) # no clear support for boundary line model. 50% support

#- e) Plot boundary line on original scale--------------------------------------

x_BD<- seq(0.25,0.5,0.001)
y_BD <- predictBL(model_BD, x_BD)
plot(exp(data_BD[,1]), data_BD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("BD/ g cm"^-3)))
lines(exp(x_BD), y_BD^2, col="red", lwd=2)


#-- 3.5. Soil pH =============================================

data_pH <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","soil_pH")]
data_pH <- data_pH[!duplicated(data_pH), ]# remove duplicate entries
data_pH <- data_pH[-which(data_pH$soil_pH==0),]

x <- data_pH$soil_pH/10 # converts bulk density from kg/m3 to g/cm3
y <- data_pH$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(y)
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_pH <- data.frame(x=x, y=sqrt(y))
out <- bagplot(data_pH, show.whiskers=F, factor = 3.0, na.rm = TRUE)
data_pH <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

# c Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# a) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_pH <- c(mean(data_pH[,1], na.rm = T), mean(data_pH[,2], na.rm = T), sd(data_pH[,1], na.rm = T), 
             sd(data_pH[,2], na.rm = T), cor(data_pH[,1],data_pH[,2], use = "complete.obs"))

# b) determine the initial model values
plot(data_pH, pch=16, col="grey")
startValues("inv-logistic",8) # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_pH <- list(c(6.17, 6.0, 9.13, dist_pH),
                 c(6.22, 5.8, 9.64, dist_pH),
                 c(6.30, 6.0, 9.78, dist_pH))

# c) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_pH, start = start_pH, sigh = sigh, model = "inv-logistic",optim.method = "Nelder-Mead")

## iv) Fit the boundary model---------

model_pH <- cbvn_extention(data = data_pH, start = start_pH, sigh = 0.2, model = "inv-logistic", 
                           pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                           xlab=expression(bold("pH")), optim.method = "Nelder-Mead")

akweight(model_pH$AIC[1,1],model_pH$AIC[2,1]) # no clear support for boundary line model. 50% support

#- e) Plot boundary line on original scale--------------------------------------

x_pH<- seq(5.5,6.5,0.001)
y_pH <- predictBL(model_pH, x_pH)
plot(data_pH[,1], data_pH[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("Soil pH")))
lines(x_pH, y_pH^2, col="red", lwd=2)


#-- 3.10. Water_holding capacity model==========================================

data_PAW <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","vwc_33kpa","vwc_1500kpa")]
data_PAW <- data_PAW[!duplicated(data_PAW), ] # remove duplicate entries
data_PAW$PAW <- data_PAW$vwc_33kpa - data_PAW$vwc_1500kpa # Calculates plant available water (PAW)
data_PAW <- data_PAW[-which(data_PAW$PAW <= 0),]

x <- data_PAW$PAW/1000 # converts PAW to cm3/cm3 from 10^-2*10
y <- data_PAW$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x)  # check for normality and transform if necessary
summastat(y)  
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

data_PAW <- data.frame(x=x, y=sqrt(y))
out <- bagplot(data_PAW, show.whiskers=F, factor = 3.0, na.rm = T)
data_PAW <- rbind(out$pxy.bag, out$pxy.outer)# removes outliers

# c Initial start values for optimizing parameters of the bivariate distribution and boundary model----

# a) First, determine the data distribution properties i.e mean of x and y, sd of x and y, and correlation of x and y
dist_PAW <- c(mean(data_PAW[,1], na.rm = T), mean(data_PAW[,2], na.rm = T), sd(data_PAW[,1], na.rm = T), 
             sd(data_PAW[,2], na.rm = T), cor(data_PAW[,1],data_PAW[,2], use = "complete.obs"))

# b) determine the initial model values
plot(data_PAW, pch=16, col="grey")
startValues("lp",8) # Determine the start values for the trapezium model i.e two intercept, two slopes and plateau value

start_PAW <- list(c(19.40, -69.33, 8.71, dist_PAW), #lp
                 c(22.47, -87.88, 9.01, dist_PAW),
                 c(21.30, -82.28, 8.94, dist_PAW))

# c) determine the standard deviation of measurement error
sigh <- c(0.1, 0.2, 0.3, 0.4, 0.5)
ble_extention(data = data_PAW, start = start_PAW, sigh = sigh, model = "lp")

## iv) Fit the boundary model---------

model_PAW <- cbvn_extention(data = data_PAW, start = start_PAW, sigh = 0.3, model = "lp", 
                           pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
                           xlab=expression(bold("PAW")))

akweight(model_PAW$AIC[1,1],model_PAW$AIC[2,1]) # no clear support for boundary line model. 50% support

#- e) Plot boundary line on original scale--------------------------------------

x_PAW<- seq(0.08,0.2,0.001)
y_PAW <- predictBL(model_PAW, x_PAW)
plot(data_PAW[,1], data_PAW[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
     xlab=expression(bold("PAW/ cm"^3*"cm"^3)))
lines(x_PAW, y_PAW^2, col="red", lwd=2)


#============= PREDICTIONS =====================================================


# 1. Crop Grid prediction ================================================================

data_pred <- read_csv("all_merged_thailand_pred_grid.csv")
names(data_pred)


# 1) Feature Engineering------------------------------------------------

weights_0_100 <- c(5, 10, 15, 30, 40)   # 0–5, 5–15, 15–30, 30–60, 60–100
total_depth <- sum(weights_0_100)       # = 100

data_pred <- data_pred %>%
  mutate(
    clay = {
      cols <- dplyr::select(., matches("^clay_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    sand = {
      cols <- dplyr::select(., matches("^sand_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    SOC = {
      cols <- dplyr::select(., matches("^socgkg_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    BD = {
      cols <- dplyr::select(., matches("^bd_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_33kpa = {
      cols <- dplyr::select(., matches("^VW33_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_1500kpa = {
      cols <- dplyr::select(., matches("^VW1500_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    soil_pH = {
      cols <- dplyr::select(., matches("^soilgrids_ph_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      (rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth)
    }
  )

# 2)  Make predictions avoiding extrapolation----------------------------------

# i) Clay---------------------------
clay_min <- min(data_clay[,1]) # the minimum from model fit
clay_max <- max(data_clay[,1]) # the maximum from model fit
clay_range <- pmin(pmax(data_pred$clay,clay_min),clay_max )# restricts the predictions to the model data range
clay <- predictBL(model_clay, clay_range)

# ii) sand--------------------------
sand_min <- min(data_sand[,1]) # the minimum from model fit
sand_max <- max(data_sand[,1]) # the maximum from model fit
sand_range <- pmin(pmax(data_pred$sand,sand_min),sand_max )# restricts the predictions to the model data range
sand <- predictBL(model_sand, sand_range)

# iii) soil bulk density--------------
BD_min <- min(data_BD[,1]) # the minimum from model fit
BD_max <- max(data_BD[,1]) # the maximum from model fit
BD_range <- pmin(pmax(log(data_pred$BD*0.001),BD_min),BD_max)# restricts the predictions to the model data range
BD <- predictBL(model_BD, BD_range)


# iv) soil organic carbon------------
SOC_min <- min(data_SOC[,1]) # the minimum from model fit
SOC_max <- max(data_SOC[,1]) # the maximum from model fit
SOC_range <- pmin(pmax(log(data_pred$SOC),SOC_min),SOC_max)# restricts the predictions to the model data range
SOC <- predictBL(model_SOC, SOC_range)


# v) soil pH------------------------
pH_min <- min(data_pH[,1]) # the minimum from model fit
pH_max <- max(data_pH[,1]) # the maximum from model fit
pH_range <- pmin(pmax(data_pred$soil_pH,pH_min),pH_max)# restricts the predictions to the model data range
pH <- predictBL(model_pH, pH_range)


# vi) Plant available water----------
PAW_min <- min(data_PAW[,1]) # the minimum from model fit
PAW_max <- max(data_PAW[,1]) # the maximum from model fit
PAW_range <- pmin(pmax(0.001*(data_pred$vwc_33kpa - data_pred$vwc_1500kpa),PAW_min),PAW_max)# restricts the predictions to the model data range
PAW <- predictBL(model_PAW, PAW_range)

# determine the most limiting

dat_out <- min_by_index(clay = clay,
                        sand = sand,
                        BD = BD,
                        SOC = SOC,
                        pH =pH,
                        PAW=PAW)


par(mar = c(10, 4, 4, 2))  # bottom, left, top, right
barplot( table(dat_out$vector_name),las = 2, ylab = "Count",xlab = "",
         main = "Frequency of identified most limiting factor", col="skyblue")

data_pred <- cbind(data_pred,dat_out)


# c.v) mapping the rainfall limited yield
Thailand <- st_read("thailand.shp") %>% st_transform(4326)
data_pred <- data_pred %>% filter(!is.na(long),!is.na(lat),!is.na(min_value))
pts <- st_as_sf(data_pred, coords = c("long", "lat"), crs = 4326)

# limiting yield-----------------

ggplot() +
  geom_sf(data = Thailand, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(data = pts, aes(color = min_value^2), size = 2, alpha = 0.9) +
  scale_color_gradient(
    low = "red",
    high = "green",
    name = "Yield (kg/ha)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(
    title = "",
    x = "Longitude",
    y = "Latitude"
  )

# limiting factor--------------------

ggplot() +
  geom_sf(data = Thailand, fill = "grey95", color = "grey60", linewidth = 0.3) +
  geom_sf(
    data = pts,
    aes(color = vector_name), #vector_name
    size = 2,
    alpha = 0.9
  ) +
  scale_color_manual(
    values = c(
      "clay"    = "red",
      "sand"    = "blue",
      "BD"    = "green",
      "SOC"    = "orange",
      "pH"    = "brown",
      "PAW"    = "black"
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
    title = "",
    x = "Longitude",
    y = "Latitude"
  )


# 2. Prediction across Thailand ==============================================

## 2.1) Yield prediction by each factor of interest ----------------------------

## a) Soil Clay-----------------------------------------------------------------

# load the clay raster files for each horizon

clay_0_5   <- rast("clay0_5.tif")
clay_5_15  <- rast("clay5_15.tif")
clay_15_30 <- rast("clay15_30.tif")
clay_30_60 <- rast("clay30_60.tif")
clay_60_100 <- rast("clay60_100.tif")

# Convert to clay content to % and then determine the depth weighted average %

w <- c(5, 10, 15,30,40)  # horizon thicknesses

clay_pct <- (clay_0_5/10*w[1] + clay_5_15/10*w[2] + clay_15_30/10*w[3]+ 
               clay_30_60/10*w[4]+ clay_60_100/10*w[5]) / sum(w)

# Avoid extrapolation. Only predict within model data range

clay_min <- min(data_clay[,1]) # the minimum from model fit
clay_max <- max(data_clay[,1]) # the maximum from model fit
clay_pct <- clamp(clay_pct, clay_min, clay_max, values=TRUE)# restricts the predictions to the model data range
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

# Mask map with Thailand extents and plot
th <- vect("Thailand.shp")
th <- project(th, crs(clay_pct))# Make sure CRS matches the clay raster
yield_clay <- mask(crop(yield_clay, th), th)
plot(yield_clay, main="Predicted Yield (t/ha) from Clay (%)")

# Clay_yield raster

out_file <- "Yield_clay_Thailand.tif"
writeRaster(
  yield_clay, filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## b) Soil sand----------------------------------------------------------------------

# load the sand raster files for each horizon
sand_0_5   <- rast("sand0_5.tif")
sand_5_15  <- rast("sand5_15.tif")
sand_15_30 <- rast("sand15_30.tif")
sand_30_60 <- rast("sand30_60.tif")
sand_60_100 <- rast("sand60_100.tif")

# Convert to clay content to % and then determine the depth weighted average %

w <- c(5, 10, 15,30,40)  # thicknesses

sand_pct <- (sand_0_5/10*w[1] + sand_5_15/10*w[2] + sand_15_30/10*w[3]+ 
               sand_30_60/10*w[4]+ sand_60_100/10*w[5]) / sum(w)

# Avoid extrapolation. Only predict within model data range

sand_min <- min(data_sand[,1]) # the minimum from model fit
sand_max <- max(data_sand[,1]) # the maximum from model fit
sand_pct <- clamp(sand_pct, sand_min, sand_max, values=TRUE)# Ensure plausible range and set nonsense to NA
names(sand_pct) <- "sand_pct"

# Function to Predict yields with parameters determined from boundary line model for sand

yield_fun_sand <- function(sand_percent) {
  y1 <- (model_sand[[3]][1,1]+ model_sand[[3]][2,1]*sand_percent)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y0 <- model_sand[[3]][3,1]
  y <- pmin(y0^2, y1^2)
  return(y)
}

# Compute the yield from the function and input raster
yield_sand <- app(sand_pct, yield_fun_sand)  # pixel-wise apply
names(yield_sand ) <- "yield_sand "

# Mask map with Thailand extents and plot
yield_sand <- mask(crop(yield_sand, th), th)
plot(yield_sand, main="Predicted Yield (t/ha) from sand (%)")

# sand_yield raster
out_file <- "Yield_Sand_Thailand.tif"
writeRaster(
  yield_sand , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)


## c) Soil Organic carbon-------------------------------------------------------

# load the sand raster files for each horizon
SOC_0_5   <- rast("SOC0_5.tif")
SOC_5_15  <- rast("SOC5_15.tif")
SOC_15_30 <- rast("SOC15_30.tif")
SOC_30_60 <- rast("SOC30_60.tif")
SOC_60_100 <- rast("SOC60_100.tif")

# Convert to clay content to % and then determine the depth weighted average %
w <- c(5, 10, 15,30,40)  # thicknesses

SOC_pct <- (SOC_0_5/10*w[1] + SOC_5_15/10*w[2] + SOC_15_30/10*w[3]+ 
              SOC_30_60/10*w[4]+ SOC_60_100/10*w[5]) / sum(w)

# Avoid extrapolation. Only predict within model data range
SOC_min <- exp(min(data_SOC[,1])) # the minimum from model fit
SOC_max <- exp(max(data_SOC[,1])) # the maximum from model fit
SOC_pct <- clamp(SOC_pct, SOC_min, SOC_max, values=TRUE)# Ensure plausible range and set nonsense to NA
names(SOC_pct) <- "SOC_pct"

# Function to Predict yields with parameters determined from boundary line model for sand

yield_fun_SOC <- function(soc_percent) {
  y1 <- model_sand[[3]][3,1] - model_SOC[[3]][1,1]* (log(soc_percent) - model_SOC[[3]][2,1])^2
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y <-  y1^2
  return(y)
}

# Compute the yield from the function and input raster
yield_SOC <- app(SOC_pct, yield_fun_SOC)  # pixel-wise apply
names(yield_SOC ) <- "yield_SOC "

# Optional: clip yields to realistic range
yield_SOC <- clamp(yield_SOC , 0, 120, values=TRUE)  # example upper cap

# Mask map with Tanzania extents
yield_SOC <- mask(crop(yield_SOC, th), th)

# Quick look
plot(yield_SOC, main="Predicted Yield (t/ha) from SOC (g/kg)")

# sand_yield raster
out_file <- "Yield_SOC_Thailand.tif"
writeRaster(
  yield_sand , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## d) Bulk density -------------------------------------------------------------

# load the bulk density raster files for each horizon

bd_0_5   <- rast("BD0_5.tif")
bd_5_15  <- rast("BD5_15.tif")
bd_15_30 <- rast("BD15_30.tif")
bd_30_60 <- rast("BD30_60.tif")
bd_60_100 <- rast("BD60_100.tif")

# Convert bulk density from cg/kg to kg/m3 and then determine the depth weighted average %
w <- c(5, 10, 15,30,40)  # thicknesses

bd_gcm <- (bd_0_5/100*w[1] + bd_5_15/100*w[2] + bd_15_30/100*w[3]+ 
             bd_30_60/100*w[4]+ bd_60_100/100*w[5]) / sum(w)


# Avoid extrapolation. Only predict within model data range
BD_min <- exp(min(data_BD[,1])) # the minimum from model fit
BD_max <- exp(max(data_BD[,1])) # the maximum from model fit
bd_gcm <- clamp(bd_gcm, BD_min, BD_max, values=TRUE)# Ensure plausible range and set nonsense to NA
names(bd_gcm) <- "bd_gcm"

# Function to Predict yields with parameters determined from boundary line model for BD

yield_fun_bd <- function(bd) {
  y1 <- (model_BD[[3]][1,1]+ model_BD[[3]][2,1]*log(bd))
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y0 <- model_BD[[3]][3,1]
  y <- pmin(y0^2, y1^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_bd <- app(bd_gcm, yield_fun_bd)  # pixel-wise apply
names(yield_bd ) <- "yield_bd "

# Mask to Thailand and plot the yield
yield_Bd <- mask(crop(yield_bd, th), th) # 
plot(yield_Bd, main=expression(bold("Predicted Yield (t/ha) from BD (g/cm"^3*")")))

# BD_yield raster
out_file <- "Yield_BD_Thailand.tif"
writeRaster(
  yield_Bd , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## e) Soil pH -------------------------------------------------------

# load the sand raster files for each horizon
pH_0_5   <- rast("pH0_5.tif")
pH_5_15  <- rast("pH5_15.tif")
pH_15_30 <- rast("pH15_30.tif")
pH_30_60 <- rast("pH30_60.tif")
pH_60_100 <- rast("pH60_100.tif")

# Convert to clay content to % and then determine the depth weighted average %
w <- c(5, 10, 15,30,40)  # thicknesses

pH <- (pH_0_5/10*w[1] + pH_5_15/10*w[2] + pH_15_30/10*w[3]+ 
             pH_30_60/10*w[4]+ pH_60_100/10*w[5]) / sum(w)

# Avoid extrapolation. Only predict within model data range
pH_min <- min(data_pH[,1]) # the minimum from model fit
pH_max <- max(data_pH[,1]) # the maximum from model fit
pH <- clamp(pH, pH_min, pH_max, values=TRUE)# Ensure plausible range and set nonsense to NA
names(pH) <- "pH"

# Function to Predict yields with parameters determined from boundary line model for sand
yield_fun_pH <- function(pH) {
  y1 <- model_pH[[3]][3,1]-(model_pH[[3]][3,1] / (1 + exp(model_pH[[3]][1,2]*(model_pH[[3]][1,1]- pH))))
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y <-  y1^2
  return(y)
}

# Compute the yield from the function and input raster
yield_pH <- app(pH, yield_fun_pH)  # pixel-wise apply
names(yield_pH ) <- "yield_pH "

# Mask map with Thailand extents and plot
yield_pH <- mask(crop(yield_pH, th), th)
plot(yield_pH, main="Predicted Yield (t/ha) from pH")

# sand_yield raster
out_file <- "Yield_pH_Thailand.tif"
writeRaster(
  yield_sand , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## f) Available Water Content --------------------------------------------------

# load the vwc_33 and vwc_1500 raster files for each horizon
# The data is Volumetric Water Content at  33kPa and 1500kPa suction in 10^-3 cm^3/cm^3 (0.1 v% or 1 mm/m)

vwc33_0_5   <- rast("VW33_0_5.tif")
vwc33_5_15  <- rast("VW33_5_15.tif")
vwc33_15_30 <- rast("VW33_15_30.tif")
vwc33_30_60 <- rast("VW33_30_60.tif")
vwc33_60_100 <- rast("VW33_60_100.tif")

vwc1500_0_5   <- rast("VW1500_0_5.tif")
vwc1500_5_15  <- rast("VW1500_5_15.tif")
vwc1500_15_30 <- rast("VW1500_15_30.tif")
vwc1500_30_60 <- rast("VW1500_30_60.tif")
vwc1500_60_100 <- rast("VW1500_60_100.tif")

# determine the PAW for each horizon

aw_0_5 <- vwc33_0_5 - vwc1500_0_5
aw_5_15 <- vwc33_5_15 - vwc1500_5_15
aw_15_30 <- vwc33_15_30 - vwc1500_15_30
aw_30_60 <- vwc33_30_60 - vwc1500_30_60
aw_60_100 <- vwc33_60_100 - vwc1500_60_100

# Determine the depth weighted average % and convert PAW to cm3/cm3
w <- c(5, 10, 15,30,40)  # thicknesses

aw <- (aw_0_5/1000*w[1] + aw_5_15/1000*w[2] + aw_15_30/1000*w[3]+
         aw_30_60/1000*w[4]+ aw_60_100/1000*w[5]) / sum(w)

# Avoid extrapolation. Only predict within model data range
PAW_min <- min(data_PAW[,1]) # the minimum from model fit
PAW_max <- max(data_PAW[,1])
aw <- clamp(aw, PAW_min, PAW_max, values=TRUE)# Ensure plausible range and set nonsense to NA
names(aw) <- "aw"

# Function to Predict yields with parameters determined from boundary line model for PAW

yield_fun_aw <- function(aw) {
  y1 <- (model_PAW[[3]][1,1]+ model_PAW[[3]][2,1]*aw)
  y1 <- ifelse( y1 < 0, 0, y1)# enforce non-negative and cap at a biological max if you like
  y0 <- model_PAW[[3]][3,1]
  y <- pmin(y0^2, y1^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_aw <- app(aw, yield_fun_aw)  # pixel-wise apply
names(yield_aw ) <- "yield_aw "

# Mask to Thailand and plot
yield_paw <- mask(crop(yield_aw, th), th) 
plot(yield_paw, main=expression(bold("Predicted Yield (t/ha) from PAW (cm"^3*"cm"^3*")")))

# PAW_yield raster

out_file <- "Yield_PAW_Thailand.tif"
writeRaster(
  yield_paw , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)

## 2.2) Law of Minimum prediction by factor -------------------------------------

# a) Minimum yield predicted

yield_min <- min(yield_clay, yield_sand, yield_pH, 
                 yield_SOC, yield_Bd, yield_paw, na.rm = TRUE)
names(yield_min) <- "yield_min_t_ha"

# Mask and plot Quick plot
yield_min <- mask(crop(yield_min, th), th) 
plot(yield_min, main = "Predicted Yield (t/ha) by law of Minimum")

out_file <- "yield_min_Thailand.tif"
writeRaster(
  yield_min , filename = out_file, overwrite = TRUE,
  datatype = "FLT4S",  # 32-bit float
  gdal = c("COMPRESS=LZW", "TILED=YES")
)


# b) Mapping the Identified Most limiting factor--------------------------------

x <- c(yield_clay, yield_sand, yield_pH, yield_SOC, yield_Bd, yield_paw)

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
  
  plot(out, plg = list(title = "Most limiting Factor", cex = 0.8),  # legend title/size
       axes = FALSE, mar = c(2,2,1,10), col=c("#DD605C","#6AE653","#D74FC8","#CAE63E","#73DED1","#929685"))
  
  par(mfrow=c(1,1))}
























