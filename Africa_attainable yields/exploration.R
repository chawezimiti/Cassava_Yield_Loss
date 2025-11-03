

plot(mask(crop(yield_min, tz), tz))


data <- data[data$country=="Tanzania",]

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha)


plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$urea_rate_kgha))

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$soil_p_mgkg))

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$disease))

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$severity))

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$incidence))

data <- data[-which(data$planting_month==2023),]
data <- data[-which(data$planting_month==2017),]
plot(as.factor(data$planting_month),data$cassava_fr_root_yld_tha, col=as.factor(data$incidence))


# Planting month-----------------------

plot(data$decimal_latitude,data$cassava_fr_root_yld_tha, col=as.factor(data$planting_month), pch=16)
plot(data$decimal_longitude, data$decimal_latitude,col=as.factor(data$planting_month), pch=16)
points(data$decimal_longitude, data$decimal_latitude,col=as.factor(data$planting_month), pch=16)

levels(as.factor(data$planting_month))


# Extract planting month and assign colours
cols <- rainbow(12)
month <- data$planting_month
col_vec <- cols[month]
# Make the scatter plot
par(mfrow=c(1,1))
plot(tz_dd)
points(data$decimal_longitude, data$decimal_latitude)

points(
  data$decimal_longitude, data$decimal_latitude,
  col = col_vec, pch = 16, cex = 1.3
)



# Add legend
legend(
  "topright",
  legend = month.abb,   # abbreviated month names
  col = cols,
  pch = 16, pt.cex = 1.2,
  title = "Month"
)

val <- data$cassava_fr_root_yld_tha
cols <- colorRampPalette(c("red", "yellow", "green"))(100)
col_index <- as.numeric(cut(val, breaks = 100, include.lowest = TRUE))
col_vec <- cols[col_index]
plot(tz_dd, main = "Cassava Fresh Root Yield (t/ha)")
points(
  data$decimal_longitude, data$decimal_latitude,
  col = col_vec, pch = 16, cex = 1.3
)
legend(
  "topright",
  legend = pretty(range(val, na.rm = TRUE), 6),
  col = colorRampPalette(c("red", "yellow", "green"))(6),
  pch = 16, pt.cex = 1.2,
  title = "Yield (t/ha)"
)
par(mfrow=c(1,1))

# disease---------------------

data$disease <- trimws(data$disease)
data$disease[data$disease == "" | is.na(data$disease)] <- "No_disease"
disease_data <- data.frame(x=as.factor(data$disease),y=data$cassava_fr_root_yld_tha)
disease_data <- do.call(rbind, lapply(split(disease_data, f=disease_data$x), outlier_remove, y=y, d=3))
rownames(disease_data) <- NULL

plot(as.factor(data$disease),data$cassava_fr_root_yld_tha) # original
plot(as.factor(disease_data$x),disease_data$y) # claen

plot(data$incidence,data$cassava_fr_root_yld_tha) # original


data$severity <- ifelse(
  data$disease == "No_disease",
  0,
  data$severity
)

## severity with boundary line----------------------------------------------------

severity_data <- data.frame(x=as.factor(data$severity),y=sqrt(data$cassava_fr_root_yld_tha))
severity_data <- do.call(rbind, lapply(split(severity_data, f=severity_data$x), outlier_remove, y=y, d=3))
rownames(severity_data) <- NULL


plot(data$severity,data$cassava_fr_root_yld_tha)# original
plot(severity_data$x,sqrt(severity_data$y))

x<-severity_data$x
y<-severity_data$y
dataxx <- data.frame(x=x, y=y)
datax <- dataxx[!(is.na(dataxx$x) | is.na(dataxx$y) | dataxx$x == "" | dataxx$y == ""), ]

ble_Cat(x=datax$x,y=datax$y,sigh=c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))
BL_Cat10_4(x=datax$x,y=datax$y,sigh=0.1, plot = TRUE)

plot(as.factor(data$severity),data$cassava_fr_root_yld_tha)# original
plot(as.factor(severity_data$x),severity_data$y)

#------------------------------

plot(data$decimal_longitude,data$cassava_fr_root_yld_tha)
plot(data$decimal_latitude,data$cassava_fr_root_yld_tha)

plot(data$decimal_longitude,data$decimal_latitude, pch=16)
plot(data$decimal_longitude,data$decimal_latitude, col=data$cassava_fr_root_yld_tha, pch=16)

library(ggplot2)

ggplot(data, aes(x = decimal_longitude, y = decimal_latitude,
                 color = cassava_fr_root_yld_tha)) +
  geom_point(size = 2) +
  scale_color_viridis_c(option = "plasma") +
  labs(
    title = "Cassava Fresh Root Yield (t/ha)",
    x = "Longitude", y = "Latitude", color = "Yield (t/ha)"
  ) +
  theme_minimal(base_size = 14)


data22 <- read_csv("master_dataset_v14.csv")

unique(data22$country)
hist(log(data$cassava_fr_root_yld_tha))


plot(tz)
points(3724587, -673793.2)
locator(1)

