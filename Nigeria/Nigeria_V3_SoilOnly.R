# This script is meant to carryout global yield gap analysis for cassava under 
# the CABI project. We employ the boundary line methodology (Webb, 1972) for this purpose.
# It used climatic, soil and management factors to determine the attainable yield 
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
library(randomcoloR)
source("Extra_functions.R") 

#---- 1. Data loading into R ===================================================

#data <- read_csv("03-merged_abbr_rainfall_tmin_tmax_(v3).csv")
data <- read_csv("Nigeria_complete.csv")
head(data)


# 2. Soil stuff first 
# d) Aggregation of Soil texture, SOC, BD, water content and pH by soil depth-------

# These soil properties are given for depths 0–5, 5–15, 15–30, 30–60, 60–100, 100-200.
# However, we assume that cassava roots go up-to 100m. A weighted average for each property 
# is determined and related to the yield.

weights_0_100 <- c(5, 10, 15, 30, 40)   # 0–5, 5–15, 15–30, 30–60, 60–100
total_depth <- sum(weights_0_100)       # = 10
data <- data %>%
  mutate(
    clay = {
      cols <- dplyr::select(., dplyr::matches("^clay_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    sand = {
      cols <- dplyr::select(., dplyr::matches("^sand_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    silt = {
      cols <- dplyr::select(., dplyr::matches("^silt_(0_5|5_15|15_30|30_60|60_100)_pct$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    SOC = {
      cols <- dplyr::select(., dplyr::matches("^socgkg_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    BD = {
      cols <- dplyr::select(., dplyr::matches("^bd_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_10kpa = {
      cols <- dplyr::select(., matches("^vw10_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_33kpa = {
      cols <- dplyr::select(., dplyr::matches("^vw33_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    vwc_1500kpa = {
      cols <- dplyr::select(., dplyr::matches("^vw1500_(0_5|5_15|15_30|30_60|60_100)$"))
      rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth
    },
    soil_pH = {
      cols <- dplyr::select(., dplyr::matches("^soilgrids_ph_(0_5cm|5_15cm|15_30cm|30_60cm|60_100cm)$"))
      (rowSums(cols * weights_0_100, na.rm = TRUE) / total_depth)
    }
  )


#3. Fit boundary lines to soil variables
#-- 3.5. Soil clay model =========================================================

dat_clay <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","clay")]
dat_clay <- dat_clay[!duplicated(dat_clay), ] # removes duplicated rows
dat_clay <- dat_clay[-which(dat_clay$clay==0),]

x <- dat_clay$clay
y <- dat_clay$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(x) # check for normality and transform if necessary
summastat(log(x)) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_clay <- data.frame(x=log(x), y=sqrt(y))
out <- bagplot(dat_clay, show.whiskers=F, factor = 3.0)
dat_clay <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error -----------------

plot(dat_clay[,1], dat_clay[,2])

#startValues("trapezium") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.


start2 <- list(c(-22, 10.61, 7.33, 29.66, -6.88, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-22.4, 10.18, 7.68, 32.94, -7.58, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               c(-31.56, 13.97, 7.74, 27.270, -6.07, mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])), 
               c(-21.77, 10.13, 8.18, 27.35, -5.91,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])))#,
               #c(-5.32, 0.59, 8.04, 32.82, -0.57,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])),
               #c(-5.12, 0.59, 7.70, 35.62, -0.63,mean(dat_clay[,1]), mean(dat_clay[,2]), sd(dat_clay[,1]), sd(dat_clay[,2]), cor(dat_clay[,1],dat_clay[,2])))


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
# why only up-to 0.8 based on variogram result for upper limit

#ble_extention(data = dat_clay, start = start2, sigh = sigh,model = "trapezium") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------


model_clay <- cbvn_extention(data=dat_clay, start = start2, sigh = 0.3, model = "trapezium", pch=16, col="grey",cex=0.7,
                             ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                             xlab = expression(bold("Clay content / log %")))

akweight(model_clay$AIC[1,1], model_clay$AIC[2,1]) # check the akaike weights for boundary and null model

# x5 <- data$clay
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# clay <- predictBL(model_clay, x5) # predicts the boundary value for each value of x in the data set
# clay <- ifelse(clay< 0,0, clay)
# 
# #- e) Plot boundary line on original scale--------------------------------------
# 
# x_clay<- seq(14,50,0.1)
# y_clay <- predictBL(model_clay, x_clay)
# plot(dat_clay[,1], dat_clay[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("Clay content / %")))
# lines(x_clay, y_clay^2, col="red", lwd=2)  


#-- 3.6. Sand model ============================================================

dat_sand <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","sand")]
dat_sand <- dat_sand[!duplicated(dat_sand), ] # remove duplicate entries
dat_sand <- dat_sand[-which(dat_sand$sand==0), ]

x <- dat_sand$sand
y <- dat_sand$cassava_fr_root_yld_tha

#- a) Summary statistics--------------------------------------------------------

summastat(log(x)) # check for normality and transform if necessary
summastat(sqrt(y))

#- b) Outlier detection---------------------------------------------------------

dat_sand <- data.frame(x=log(x), y=sqrt(y))
out <- bagplot(dat_sand, show.whiskers=F, factor = 3.0)
dat_sand <- rbind(out$pxy.bag, out$pxy.outer) # removes outliers

#- c) Initial starting values for model and estimate of the sd of measurement error----------------

plot(dat_sand[,1], dat_sand[,2])
#startValues("trapezium")# determine start values by clicking on the plot the points that make up the desired model.


 start2 <- list(c(-41.4, 12.81, 8.11, 88.86, -19.79, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
                c(-26.25, 8.52, 7.45, 67.908, -14.134, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
                c(-41.9, 12.78, 8.383, 90.06, -19.73, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])))#,
#                c(2.97, 0.11, 7.99, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
#                c(1.41, 0.19, 7.73, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
#                c(1.09, 0.16, 7.83, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
#                c(2.35, 0.11, 7.88, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])),
#                c(1.78, 0.13, 7.79, mean(dat_sand[,1]), mean(dat_sand[,2]), sd(dat_sand[,1]), sd(dat_sand[,2]), cor(dat_sand[,1],dat_sand[,2])))





sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8) # vector of possible values for sd of measurement error
#ble_extention(data = dat_sand, start = start2, sigh = sigh, model = "trapezium" ) # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_sand <- cbvn_extention(data=dat_sand, start = start2, sigh = 0.3, model = "trapezium", pch=16, col="grey", cex=0.7,
                             ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                             xlab = expression(bold("sand content /log( %)")))

akweight(model_sand$AIC[1,1], model_sand$AIC[2,1])# check the akaike weights for boundary and null model

# x5 <- data$sand
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# sand <- predictBL(model_sand, x5) # predicts the boundary value for each value of x in the data set

#- e) Plot boundary line on original scale--------------------------------------

# x_sand<- seq(35,80,0.1)
# y_sand <- predictBL(model_sand, x_sand)
# plot(dat_sand[,1], dat_sand[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("sand content / %")))
# lines(x_sand, y_sand^2, col="red", lwd=2)  


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
#startValues("lp") # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(#c(-11.6, 3.4, 7.82, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-9.6, 2.9, 7.82, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-10.6, 2.8, 7.82, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])),
               c(-3.86, 2.09, 8.22, mean(data_pH[,1]), mean(data_pH[,2]), sd(data_pH[,1]), sd(data_pH[,2]), cor(data_pH[,1],data_pH[,2])))

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6)# vector of possible values for sd of measurement error
#ble_extention(data = data_pH, start = start2, sigh = sigh, model = "lp") # determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_pH <- cbvn_extention(data=data_pH, start = start2, sigh = 0.4, model = "lp", pch=20, col="grey",
                           ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                           xlab = expression(bold("Soil pH")))

akweight(model_pH$AIC[1,1], model_pH$AIC[2,1]) # check the akaike weights for boundary and null model

# x5 <- data$soil_pH
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# pH <- predictBL(model_pH , x5) # predicts the boundary value for each value of x in the data set
# pH <- ifelse(pH< 0,0, pH)
# 
# #- e) Plot boundary line on original scale---------------------------------------
# 
# x_pH<- seq(4.7,7.3,0.01)
# y_pH <- predictBL(model_pH, x_pH)
# plot(data_pH[,1], data_pH[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("soil pH")))
# lines(x_pH, y_pH^2, col="red", lwd=2)  

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
#startValues("schmidt",6) # determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 2 points are required. Determine multiple values and 
# insert them into the list as below.

#start2 <- list(c(10.46, 3.2, 7.77, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))#,
             #  c(43.20, 3.1, 7.75, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))#,
            #   c(36.19, 3.11, 8.71, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))#,
             #  c(-19.97, 14.31, 7.89, 29.75, -7.68, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
             #  c(-18.73, 13.64, 7.78, 28.30, -7.28, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))

#startValues("lp")
start2 <- list(c(13.66, -1.92, 7.88, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])),
               c(14.36, -2.04, 8.29, mean(data_SOC[,1]), mean(data_SOC[,2]), sd(data_SOC[,1]), sd(data_SOC[,2]), cor(data_SOC[,1],data_SOC[,2])))
               


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
#ble_extention(data = data_SOC, start = start2, sigh = sigh, model = "lp")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value


#- d) Model fitting-------------------------------------------------------------

model_SOC <- cbvn_extention(data=data_SOC, start = start2, sigh = 0.1, model = "lp", pch=16, col="grey",
                            ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                            xlab = expression(bold("SOC/log (%)")))



akweight(model_SOC$AIC[1,1], model_SOC$AIC[2,1]) # check the akaike weights for boundary and null model

# x5 <- log(data$SOC+0.000001)
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# SOC <- predictBL(model_SOC , x5) # predicts the boundary value for each value of x in the data set
# SOC <- ifelse(SOC< 0,0, SOC )
# 
# #- e) Plot boundary line on original scale--------------------------------------
# 
# x_SOC<- seq(1.6,3.5,0.01)
# y_SOC <- predictBL(model_SOC, x_SOC)
# plot(exp(data_SOC[,1]), data_SOC[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("SOC")))
# lines(exp(x_SOC), y_SOC^2, col="red", lwd=2)

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
#startValues("lp")# determine start values by clicking on the plot the points that make up the desired model.
# In this case for the trapezium, 4 points are required. Determine multiple values and 
# insert them into the list as below.

start2 <- list(c(48.5, -28.06, 7.95, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])),
               c(49.8, -28.96, 8.2, mean(data_BD[,1]), mean(data_BD[,2]), sd(data_BD[,1]), sd(data_BD[,2]), cor(data_BD[,1],data_BD[,2])))#,
           

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
#ble_extention(data = data_BD, start = start2, sigh = sigh, model = "lp")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_BD <- cbvn_extention(data=data_BD, start = start2, sigh = 0.3, model = "lp", pch=16, col="grey",
                           ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                           xlab = expression(bold("BD g/kg")))

akweight(model_BD$AIC[1,1], model_BD$AIC[2,1]) # check the akaike weights for boundary and null model

# x5 <- data$BD/1000
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# BD <- predictBL(model_BD , x5) # predicts the boundary value for each value of x in the data set
# BD <- ifelse(BD< 0,0, BD )
# 
# #- e) Plot boundary line on original scale--------------------------------------
# 
# x_BD<- seq(1.25,1.45,0.001)
# y_BD <- predictBL(model_BD, x_BD)
# plot(data_BD[,1], data_BD[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("BD g/cm"^3)))
# lines(x_BD, y_BD^2, col="red", lwd=2)  


#-- 3.10. Water_holding capacity model==========================================

data_PAW <- data[, c("decimal_longitude","decimal_latitude","cassava_fr_root_yld_tha","vwc_33kpa","vwc_1500kpa")]
data_PAW <- data_PAW[!duplicated(data_PAW), ] # remove duplicate entries
data_PAW$PAW <- data_PAW$vwc_33kpa - data_PAW$vwc_1500kpa # Calculates plant available water (PAW)
data_PAW <- data_PAW[which(data_PAW$PAW > 0),]

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

#start2 <- list(c(1.85, 56.78, 7.93, 22.11, -86.1, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
#              c(2.85, 39.23, 7.39, 18.8, -70.79, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])))#,
 #              c(1.72, 100.73, 7.59, 39.10, -150.39, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
              # c(-1.36, 105.42, 7.68, 33.97, -127.24, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
              # c(-0.34, 83.93, 7.77, 34.76,-130.93, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
              # c(-1.02, 97.83, 7.55, 26.44, -94.20, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])))

start2 <- list(c(-1.45, 85.78, 8.08, 22.11, -90.1, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2]))) #
               #c(2.85, 39.23, 7.39, 18.8, -70.79, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2]))) #,
               #c(1.72, 100.73, 7.59, 39.10, -150.39, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               # c(-1.36, 105.42, 7.68, 33.97, -127.24, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               # c(-0.34, 83.93, 7.77, 34.76,-130.93, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])),
               # c(-1.02, 97.83, 7.55, 26.44, -94.20, mean(data_PAW[,1]), mean(data_PAW[,2]), sd(data_PAW[,1]), sd(data_PAW[,2]), cor(data_PAW[,1],data_PAW[,2])))
               
               
sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)# vector of possible values for sd of measurement error
#ble_extention(data = data_PAW, start = start2, sigh = sigh, model = "trapezium")# determines the maximum likelihood for each  sd of measurement error value
# Select one with the largest likelihood value

#- d) Model fitting-------------------------------------------------------------

model_paw <- cbvn_extention(data=data_PAW, start = start2, sigh = 0.5, model = "trapezium", pch=16, col="grey",
                            ylab = expression(bold("Yield / " * sqrt(t~ha^-1))), 
                            xlab = expression(bold("PAW cm"^-3*"cm"^-3)))





akweight(model_paw$AIC[1,1], model_paw$AIC[2,1]) # check the akaike weights for boundary and null model

# x5 <- (data$vwc_33kpa - data$vwc_1500kpa)/1000
# x5[which(is.na(x5)==T)] <- mean(x5, na.rm = T)
# paw <- predictBL(model_paw , x5) # predicts the boundary value for each value of x in the data set
# paw <- ifelse(paw< 0,0, paw )
# 
# #- e) Plot boundary line on original scale--------------------------------------
# 
# x_paw<- seq(0.06,0.23,0.001)
# y_paw <- predictBL(model_paw, x_paw)
# plot(data_PAW[,1], data_PAW[,2]^2, pch=16, col="grey", ylab=expression(bold("Yield / t ha"^-1)), 
#      xlab=expression(bold("PAW cm"^-3*"cm"^-3)))
# lines(x_paw, y_paw^2, col="red", lwd=2)  


#--- 5. Mapping attainable yields and the most-limiting factors for Nigeria =========

## 5.1) Yield prediction by each factor of interest for Nigeria ----------------------------

bbox <- c(xmin = 2, xmax = 15, ymin = 4, ymax = 15)
depths <- c("0-5cm","5-15cm","15-30cm","30-60cm","60-100cm")
w <- c(5, 10, 15, 30, 40)# depth weights in cm (0–5, 5–15, 15–30, 30–60, 60–100)
tz <- vect("ng.shp") # Tanzania shape file for extents
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
#AEM Clamped to 1 and 100 so we can take logs 
clay_pct <- clamp(clay_pct, 1, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
names(clay_pct) <- "clay_pct"

clay_pct=log(clay_pct)

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

 out_file <- "yield_from_clay_Nigeria.tif"
 writeRaster(
   yield_clay, filename = out_file, overwrite = TRUE,
   datatype = "FLT4S",  # 32-bit float
   gdal = c("COMPRESS=LZW", "TILED=YES")
 )

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
sand_pct <- clamp(sand_pct, 1, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
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


sand_pct=log(sand_pct)

# Compute the yield from the function and input raster

yield_sand <- app(sand_pct, yield_fun_sand)  # pixel-wise apply
names(yield_sand ) <- "yield_sand "

# Optional: clip yields to realistic range

yield_sand  <- clamp(yield_sand , 0, 120, values=TRUE)  # example upper cap

# Mask map with Tanzania extents

yield_sand <- mask(crop(yield_sand, tz), tz)

# sand_yield raster

out_file <- "yield_from_sand_Nigeria.tif"
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
  #y2 <- (model_pH[[3]][4,1]+ model_pH[[3]][5,1]*pH_percent)
  y0 <- model_pH[[3]][3,1]
  
  y <- pmin(y0^2, y1^2)
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

 out_file <- "yield_from_pH_Nigeria.tif"
 writeRaster(
   yield_pH , filename = out_file, overwrite = TRUE,
   datatype = "FLT4S",  # 32-bit float
   gdal = c("COMPRESS=LZW", "TILED=YES")
 )

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
#AEM clamped to 1 - 100 as variate will be logged
soc_pct <- clamp(soc_pct, 1, 100, values=TRUE)# Ensure plausible range and set nonsense to NA
names(soc_pct) <- "soc_pct"
soc_pct <-log(soc_pct)

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

out_file <- "yield_from_soc_Nigeria.tif"
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
#  y2 <- (model_BD[[3]][4,1]+ model_BD[[3]][5,1]*bd)
  y0 <- model_BD[[3]][3,1]
  
  y <- pmin(y0^2, y1^2)
  return(y)
}

# Compute the yield from the function and input raster

yield_bd <- app(bd, yield_fun_bd)  # pixel-wise apply
names(yield_bd ) <- "yield_bd "

# Optional: clip yields to realistic range

yield_bd  <- clamp(yield_bd , 0, 120, values=TRUE)  # example upper cap

yield_bd <- mask(crop(yield_bd, tz), tz)

# BD_yield raster

out_file <- "yield_from_bd_Nigeria.tif"
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

out_file <- "yield_from_aw_Nigeria.tif"
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

plot(yield_min, main = "Predicted Yield (t/ha) by law of Minimum: Soil properties")

out_file <- "yield_min_Nigeria.tif"
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






