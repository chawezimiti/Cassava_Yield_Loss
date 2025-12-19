predictBL_boot <- function(object, x, B = 1000, conf = 0.95) {
  
  # ---- parameter estimates ----
  theta_hat <- object$Parameters[,1]
  H         <- object$Hessian
  
  # ---- covariance ----
  Sigma <- tryCatch(
    solve(H),
    error = function(e) stop("Hessian is singular.")
  )
  
  # ---- bootstrap parameters ----
  theta_sim <- MASS::mvrnorm(B, mu = theta_hat, Sigma = Sigma)
  
  # ---- prediction function ----
  predict_theta <- function(theta) {
    tmp <- object
    tmp$Parameters[,1] <- theta
    predictBL(tmp, x)
  }
  
  # ---- bootstrap predictions ----
  Y <- apply(theta_sim, 1, predict_theta)
  # dimensions: length(x) × B
  
  # ---- summary statistics ----
  fit <- apply(Y, 1, mean, na.rm = TRUE)
  
  alpha <- 1 - conf
  lower <- apply(Y, 1, quantile, probs = alpha/2, na.rm = TRUE)
  upper <- apply(Y, 1, quantile, probs = 1 - alpha/2, na.rm = TRUE)
  
  data.frame(
    x     = x,
    fit   = fit,
    lower = lower,
    upper = upper
  )
}


rain_early2 <- predictBL_boot(model_pmm_est, log(data$Filled_pmm_early))
rain_veg2 <- predictBL_boot(model_pmm_veg, data$Filled_pmm_vegetative)
rain_bulking2 <- predictBL_boot(model_pmm_bulk, log(data$Filled_pmm_bulking))
rain_maturation2 <- predictBL_boot(model_pmm_mat, sqrt(data$Filled_pmm_maturation))

spei_early2 <- predictBL_boot(model_spei_est, data$spei_early)
spei_veg2 <- predictBL_boot(model_spei_veg, data$spei_vegetative)
spei_bulking2 <- predictBL_boot(model_spei_bulk, data$spei_bulking)
spei_maturation2 <- predictBL_boot(model_spei_mat, data$spei_maturation)


rain_early2 <- rain_early2[order(rain_early2$x), ]
rain_veg2 <- rain_veg2[order(rain_veg2$x), ]
rain_bulking2 <- rain_bulking2[order(rain_bulking2$x), ]
rain_maturation2 <- rain_maturation2[order(rain_maturation2$x), ]
spei_early2 <- spei_early2[order(spei_early2$x), ]
spei_veg2 <- spei_veg2[order(spei_veg2$x), ]
spei_bulking2 <-spei_bulking2[order(spei_bulking2$x), ]
spei_maturation2 <- spei_maturation2[order(spei_maturation2$x), ]

par(mfrow=c(1,1))

plot(rain_early2$x, rain_early2$fit, ylim=c(4.5,8.5))
lines(rain_early2$x, rain_early2$lower, col="red")
lines(rain_early2$x, rain_early2$upper, col="red")

plot(rain_veg2$x, rain_veg2$fit, ylim=c(4.5,8.5))
lines(rain_veg2$x, rain_veg2$lower, col="red")
lines(rain_veg2$x, rain_veg2$upper, col="red")

plot(rain_bulking2$x, rain_bulking2$fit, ylim=c(4.5,8.5))
lines(rain_bulking2$x, rain_bulking2$lower, col="red")
lines(rain_bulking2$x, rain_bulking2$upper, col="red")

plot(rain_maturation2$x, rain_maturation2$fit, ylim=c(4.5,8.5))
lines(rain_maturation2$x, rain_maturation2$lower, col="red")
lines(rain_maturation2$x, rain_maturation2$upper, col="red")

plot(spei_early2$x, spei_early2$fit, ylim=c(4.5,10))
lines(spei_early2$x, spei_early2$lower, col="red")
lines(spei_early2$x, spei_early2$upper, col="red")

plot(spei_veg2$x, spei_veg2$fit, ylim=c(4.5,8.5))
lines(spei_veg2$x, spei_veg2$lower, col="red")
lines(spei_veg2$x, spei_veg2$upper, col="red")

plot(spei_bulking2$x, spei_bulking2$fit, ylim=c(4.5,8.5))
lines(spei_bulking2$x, spei_bulking2$lower, col="red")
lines(spei_bulking2$x, spei_bulking2$upper, col="red")

plot(spei_maturation2$x, spei_maturation2$fit, ylim=c(4.5,10))
lines(spei_maturation2$x, spei_maturation2$lower, col="red")
lines(spei_maturation2$x, spei_maturation2$upper, col="red")







