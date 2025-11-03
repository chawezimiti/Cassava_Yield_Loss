
library(readr)
library(BLA)
library(aplpack)
source("Extra_functions.R")


## Read in data-----------------------------------------------------------------

data <- read_csv("master_dataset_v14.csv")
head(data)


## Analysis---------------------------------------------------------------------

# The attainable yield are modeled based on the latitude. The distance away from 
# the latitude 0. To check the influence of climate on cassava yields.

# Create a data frame and standardize the missing  values (either NA or blank spaces) to NA

dat <- data.frame(x=data$decimal_latitude,y=data$cassava_fr_root_yld_tha)
dat <- data.frame(x,y)
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

start2 <- list(c(9.47, 0.20, 8.02, 14.44,-0.42, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),        
               c(9.13, 0.19, 8.07, 14.81,-0.44, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),
               c(9.32, 0.20, 8.21, 12.96,-0.31, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")))

# determining of the standard deviation of measurement error

sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6)
ble_extention(data = dat_clean, start = start2, sigh = sigh, model = "trapezium")

# model fitting 

par(mar=c(6,5,4,2)) # plot size parameter adjustment

cbvn_extention(data=dat_clean, start = start2, sigh = 0.3, model = "trapezium", pch=16, col="grey",
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
               c(15, 0.44, 8.07, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")),
               c(13, 0.31, 8.21, mean(dat_clean[,1], na.rm = T), mean(dat_clean[,2], na.rm = T), sd(dat_clean[,1], na.rm = T), sd(dat_clean[,2], na.rm = T), cor(dat_clean[,1],dat_clean[,2], use = "complete.obs")))


sigh <- c(0.1,0.2,0.3,0.4,0.5,0.6) #possible measurement error values

ble_extention(data = dat_clean, start = start2, sigh = sigh, model = "other", equation=trap_mirror) #determine measurement error by maximum likelihood


cbvn_extention(data = dat_clean, start = start2, sigh = 0.3, model = "other", equation=trap_mirror, pch=16, col="grey",
               ylab=expression(bold("Yield / t ha"^-1)), 
               xlab=expression(bold("Latitude")),
               main="Mirror Model")



