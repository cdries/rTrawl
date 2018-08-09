## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

library(rTrawl)

## ----trawl intuition, echo=FALSE, out.width = "49%"----------------------

set.seed(2018)
n <- 25
mm <- 0
MM <- 10
lo <- 251
gr <- seq(mm, MM, length.out = lo)
x <- sort(runif(n, mm, MM))
y <- runif(n)
jumps <- rep(1, n)

x_ext <- c(x, x - log(y))
jumps_ext <- c(jumps, -jumps)
ord <- order(x_ext)
x_ext <- x_ext[ord]
jumps_ext <- jumps_ext[ord]

gr_ord <- (1:length(gr) + round(length(gr) / 2)) %% length(gr) + 1
ii <- gr_ord[1]

# plot Lévy seed
plot(x, y, ylim = c(0, 1), xlim = c(mm, MM),
     las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5,
     type = 'n', xlab = "time", ylab = "seed")

xtemp <- seq(0, gr[ii], length.out = 251)
ytemp <- exp(-(gr[ii] - xtemp))
polygon(c(xtemp, rev(xtemp)), c(rep(-0.005, length(xtemp)), rev(ytemp) + 0.005), 
        col = scales::alpha("red", 0.4), border = NA)
points(x[jumps > 0], y[jumps > 0], pch = 16)
points(x[jumps < 0], y[jumps < 0])
abline(v = gr[ii], lwd = 2)

# plot process values
xtemp <- c(0, x_ext[x_ext <= gr[ii]])
ztemp <- c(0, cumsum(jumps_ext[x_ext <= gr[ii]]))
plot(c(xtemp, gr[ii]), c(ztemp, ztemp[length(ztemp)]), ylim = c(0, 10), xlim = c(mm, MM),
     las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5,
     type = 's', lwd = 2, xlab = "time", ylab = "process value")
abline(v = gr[ii], lwd = 2)

## ----compute bases-------------------------------------------------------
sim_Poisson <- sim_trawl(list("levy_seed" = "Poisson", "levy_par" = 1.1, "b" = 1, "T0" = 0, "TT" = 10))
sim_Skellam <- sim_trawl(list("levy_seed" = "Skellam", "levy_par" = c(3, 2.5), "b" = 1, "T0" = 0, "TT" = 10))
sim_negBin <- sim_trawl(list("levy_seed" = "negBin", "levy_par" = c(2, 0.4), "b" = 1, "T0" = 0, "TT" = 10))
sim_DnegBin <- sim_trawl(list("levy_seed" = "DnegBin", "levy_par" = c(2, 0.4, 1.8, 0.5), "b" = 1, "T0" = 0, "TT" = 10))

## ----plot bases, echo=FALSE, out.width = "49%"---------------------------
plot(sim_Poisson$x_grid, sim_Poisson$p_grid, xlab = "time", ylab = "process value", main = "Poisson",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_Skellam$x_grid, sim_Skellam$p_grid, xlab = "time", ylab = "process value", main = "Skellam",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_negBin$x_grid, sim_negBin$p_grid, xlab = "time", ylab = "process value", main = "Negative Binomial",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_DnegBin$x_grid, sim_DnegBin$p_grid, xlab = "time", ylab = "process value", main = "Delta Negative Binomial",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

## ----pure trawl example--------------------------------------------------
sim_Poisson <- sim_trawl(list("levy_seed" = "Poisson", "levy_par" = 4, "trawl" =  "exp", 
                              "trawl_par" = 0.7, "b" = 0, "T0" = 0, "TT" = 30))
sim_Skellam <- sim_trawl(list("levy_seed" = "Skellam", "levy_par" = c(3, 2.5), "trawl" =  "gamma", 
                              "trawl_par" = c(1.1, 1.05), "b" = 0, "T0" = 0, "TT" = 30))
sim_negBin <- sim_trawl(list("levy_seed" = "negBin", "levy_par" = c(8, 0.7), "trawl" = "invGauss",
                             "trawl_par" = c(0.3, 0.6), "b" = 0, "T0" = 0, "TT" = 30))
sim_DnegBin <- sim_trawl(list("levy_seed" = "DnegBin", "levy_par" = c(3, 0.4, 2.8, 0.5), "trawl" = "gig",
                              "trawl_par" = c(0.1, 0.4, -0.7), "b" = 0, "T0" = 0, "TT" = 30))

## ----plot trawl example, echo=FALSE, out.width = "49%"-------------------
plot(sim_Poisson$x_grid, sim_Poisson$p_grid, xlab = "time", ylab = "process value", main = "Poisson",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_Skellam$x_grid, sim_Skellam$p_grid, xlab = "time", ylab = "process value", main = "Skellam",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_negBin$x_grid, sim_negBin$p_grid, xlab = "time", ylab = "process value", main = "Negative Binomial",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_DnegBin$x_grid, sim_DnegBin$p_grid, xlab = "time", ylab = "process value", main = "Delta Negative Binomial",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

## ----mvSkellam trawl example---------------------------------------------
trawl <- list("gamma", "exp")
trawl_par <- list(c(0.9, 1.8), c(0.5))
levy_par <- matrix(c(0.13, 0.23, 0.11), ncol = 1)
design_matrix <- matrix(c(1, 0, 0, 1, -1, -1), nrow = 2)
sim_mvSkellam <- sim_trawl(list("levy_seed" = "Skellam", "levy_par" = levy_par,
                                "trawl" = trawl, "trawl_par" = trawl_par,
                                "design_matrix" = design_matrix, "b" = c(0, 0),
                                "T0" = 0, "TT" = 3600, "observed_freq" = 1e-6), 
                           univariate = FALSE)

## ----plot mvSkellam trawl example, echo=FALSE, out.width = "49%"---------
plot(sim_mvSkellam$x_grid[[1]], sim_mvSkellam$p_grid[[1]], xlab = "time", ylab = "process value", main = "Skellam - component 1",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_mvSkellam$x_grid[[2]], sim_mvSkellam$p_grid[[2]], xlab = "time", ylab = "process value", main = "Skellam - component 2",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

## ----mvSkellamNegbin trawl example---------------------------------------
trawl <- list("gamma", "exp")
trawl_par <- list(c(0.9, 1.8), c(0.5))
levy_par <- matrix(c(1.5, 0.4, 0.3,
                     0.7, 0.4, NA,
                     0.6, 0.3, NA), ncol = 3, byrow = TRUE)
design_matrix <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 2)
sim_mvNegBin <- sim_trawl(list("levy_seed" = "negBin", "levy_par" = levy_par,
                               "trawl" = trawl, "trawl_par" = trawl_par,
                               "design_matrix" = design_matrix, "b" = c(0, 0),
                               "T0" = 15.3, "TT" = 1200, "observed_freq" = 1e-3), 
                          univariate = FALSE)

## ----plot mvSkellamNegbin trawl example, echo=FALSE, out.width = "49%"----
plot(sim_mvNegBin$x_grid[[1]], sim_mvNegBin$p_grid[[1]], xlab = "time", ylab = "process value", main = "Negative binomial - component 1",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
plot(sim_mvNegBin$x_grid[[2]], sim_mvNegBin$p_grid[[2]], xlab = "time", ylab = "process value", main = "Negative binomial - component 2",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

## ----fit univariate trawl example----------------------------------------
sim <- sim_trawl(list("levy_seed" = "Poisson", "levy_par" = 4, "trawl" =  "exp", 
                      "trawl_par" = 0.7, "b" = 0, "T0" = 0, "TT" = 3600))
sim$h <- 0.5
sim$trawl <- "exp"
sim$lag_max <- 3
sim$method <- "acf"
ft <- fit_trawl(sim)
print(ft$trawl_par)

## ----plot fit univariate trawl example, echo=FALSE, out.width = "50%"----
plot(sim$x_grid, sim$p_grid, xlab = "time", ylab = "process value", main = "Poisson",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

plot(1:10, acf_trawl(sim, sim$h, lag_max = 10), xlab = "lag", ylab = "autocorrelation", main = "ACF",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ft, sim$h, method = "acf", lag_max = 10), col = "blue", lwd = 2)

## ----fit univariate price trawl example----------------------------------
sim <- sim_trawl(list(levy_seed = "Skellam", levy_par = c(0.13, 0.11), b = 0.3, TT = 1200))
sim$h <- exp(seq(log(1e-2), log(60), length.out = 51))
sim$trawl <- "exp"
sim$include_b <- TRUE
sim$include_cum1 <- FALSE

sim$method <- "vs_C"
ftC <- fit_trawl(sim)
print(ft$trawl_par)

sim$method <- "vs_SY"
ftSY <- fit_trawl(sim)
print(ft$trawl_par)

## ----plot fit univariate price trawl example, echo=FALSE, out.width = "50%"----
plot(sim$x_grid, 100 + sim$p_grid, xlab = "time", ylab = "process value", main = "Skellam price process",
     type = 's', lwd = 2, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)

plot(log(sim$h), vs_trawl(sim, sim$h), xlab = "log(observation frequency)", 
     ylab = "Realized variance / observation frequency", 
     main = "Variance signature plot",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(log(sim$h), vs_trawl(ftC, sim$h, method = "vs_C"), col = "blue", lwd = 2)
lines(log(sim$h), vs_trawl(ftSY, sim$h, method = "vs_SY"), col = "red", lwd = 2, lty = 1)
legend("topright", c("vs_C", "vs_SY"), lwd = c(2, 2), col = c("blue", "red"), lty = c(1, 1))

## ----univariate step wise------------------------------------------------
sim <- sim_trawl(list("levy_seed" = "Skellam", "levy_par" = c(0.131, 0.130),
                      "trawl" = "gig", "trawl_par" = c(0.3, 0.45, -0.6),
                      "T0" = 72.03, "TT" = 75600, "observed_freq" = 1e-6, 
                      "b" = 0))

# fit the trawl parameters
sim$h <- 0.5
sim$trawl <- "gig"
sim$lag_max <- 3
sim$method <- "acf"
ft <- fit_trawl(sim) 
print(ft$trawl_par)

# fit the Skellam basis
lv_fit <- fit_levy(list("levy_seed" = "Skellam", "trawl" = "gig",
                        "trawl_par" = ft$trawl_par, "T0" = 72.03, "TT" = 75600,
                        "x_grid" = sim$x_grid, "p_grid" = sim$p_grid))  
print(lv_fit)

## ----multivariate Skellam------------------------------------------------
trawl <- list("gamma", "exp")
trawl_par <- list(c(0.9, 1.8), c(0.5))
levy_par <- matrix(c(0.13, 0.13, 0.23, 0.11, 0.05), ncol = 1)
design_matrix <- matrix(c(1, 0, 0, 1, -1, 1, 1, -1, 1, 0), nrow = 2)
sm <- sim_trawl(list("levy_seed" = "Poisson", "levy_par" = levy_par, "trawl" = trawl, "trawl_par" = trawl_par,
                     "design_matrix" = design_matrix, "b" = c(0, 0), 
                     "T0" = 35, "TT" = 24600, "observed_freq" = 1e-6),
                univariate = FALSE)

h <- 1
lag_max <- 2              
ft1 <- fit_trawl(list("T0" = 35, "TT" = 24600, "method" = "acf", "h" = h, "trawl" = trawl[[1]],
                      "x_grid" = sm$x_grid[[1]], "p_grid" = sm$p_grid[[1]],
                      "lag_max" = lag_max, "b" = 0), multi = 1)
print(ft1$trawl_par)
ft2 <- fit_trawl(list("T0" = 35, "TT" = 24600, "method" = "acf", "h" = h, "trawl" = trawl[[2]],
                      "x_grid" = sm$x_grid[[2]], "p_grid" = sm$p_grid[[2]],
                      "lag_max" = lag_max, "b" = 0), multi = 1)
print(ft2$trawl_par)

lv_fit <- fit_levy(list("levy_seed" = "Poisson", "trawl" = trawl, 
                        "trawl_par" = list(ft1$trawl_par, ft2$trawl_par), "T0" = 35, "TT" = 24600,
                        "x_grid" = sm$x_grid, "p_grid" = sm$p_grid, 
                        "design_matrix" = design_matrix))
print(lv_fit)

## ----multivariate negative binomial--------------------------------------
trawl <- list("gamma", "exp")
trawl_par <- list(c(0.9, 1.8), c(0.5))
levy_par <- matrix(c(1.5, 0.4, 0.3,
                     0.7, 0.4, NA,
                     0.6, 0.3, NA), ncol = 3, byrow = TRUE)
design_matrix <- matrix(c(1, 1, 1, 0, 0, 1), nrow = 2)
sm <- sim_trawl(list("levy_seed" = "negBin", "levy_par" = levy_par,
                     "trawl" = trawl, "trawl_par" = trawl_par,
                     "design_matrix" = design_matrix, "b" = c(0, 0),
                     "T0" = 0, "TT" = 25400, "observed_freq" = 1e-6),
                univariate = FALSE)

h <- 1       
lag_max <- 2       
ft1 <- fit_trawl(list("T0" = 0, "TT" = 25400, "method" = "acf", "h" = h, "trawl" = trawl[[1]],
                      "x_grid" = sm$x_grid[[1]], "p_grid" = sm$p_grid[[1]],
                      "lag_max" = lag_max, "b" = 0), multi = 1)
ft2 <- fit_trawl(list("T0" = 0, "TT" = 25400, "method" = "acf", "h" = h, "trawl" = trawl[[2]],
                      "x_grid" = sm$x_grid[[2]], "p_grid" = sm$p_grid[[2]],
                      "lag_max" = lag_max, "b" = 0), multi = 1)

lv_fit <- fit_levy(list("levy_seed" = "negBin", "trawl" = trawl,
                        "trawl_par" = list(ft1$trawl_par, ft2$trawl_par), "T0" = 0, "TT" = 25400,
                        "x_grid" = sm$x_grid, "p_grid" = sm$p_grid, 
                        "design_matrix" = design_matrix))

## ----xts, message=FALSE--------------------------------------------------
library(xts)
data(ESM8)

## ----ESM8 visual 1, echo=FALSE, out.width = "50%"------------------------
plot(ESM8, type = 's')
plot(diff(ESM8), pch = 16)

## ----ESM8 make object----------------------------------------------------
p_grid <- as.numeric(ESM8$Price) * 4
x_grid <- as.numeric(difftime(index(ESM8), as.POSIXct("2018-05-07 08:30:00", tz = "America/Chicago"), 
                              tz = "America/Chicago", units = "secs"))
T0 <- 0
TT <- 12600
trs <- list("x_grid" = x_grid, "p_grid" = p_grid, "T0" = T0, "TT" = TT)

## ----ESM8 estimate vs----------------------------------------------------
h <- exp(seq(log(1e-3), log(5 * 60), length.out = 51))

ftSY_exp <- fit_trawl(list("method" = "vs_SY", "trawl" = "exp", "x_grid" = x_grid, 
                           "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftSY_gamma <- fit_trawl(list("method" = "vs_SY", "trawl" = "gamma", "x_grid" = x_grid, 
                             "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftSY_invGauss <- fit_trawl(list("method" = "vs_SY", "trawl" = "invGauss", "x_grid" = x_grid, 
                                "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftSY_gig <- fit_trawl(list("method" = "vs_SY", "trawl" = "gig", "x_grid" = x_grid, 
                           "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))

ftC_exp <- fit_trawl(list("method" = "vs_C", "trawl" = "exp", "x_grid" = x_grid, 
                          "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftC_gamma <- fit_trawl(list("method" = "vs_C", "trawl" = "gamma", "x_grid" = x_grid, 
                            "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftC_invGauss <- fit_trawl(list("method" = "vs_C", "trawl" = "invGauss", "x_grid" = x_grid, 
                               "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))
ftC_gig <- fit_trawl(list("method" = "vs_C", "trawl" = "gig", "x_grid" = x_grid, 
                          "p_grid" = p_grid, "T0" = T0, "TT" = TT, "h" = h))

## ----ESM8 vs fittings, echo=FALSE, out.width = "49%"---------------------
plot(log(h), vs_trawl(trs, h), xlab = "log(observation frequency)", 
     ylab = "Realized variance / observation frequency", 
     main = "Variance signature plot - vs_SY",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(log(h), vs_trawl(ftSY_exp, h, method = "vs_SY"), col = "red", lwd = 2)
lines(log(h), vs_trawl(ftSY_gamma, h, method = "vs_SY"), col = "blue", lwd = 2)
lines(log(h), vs_trawl(ftSY_invGauss, h, method = "vs_SY"), col = "orange", lwd = 2)
lines(log(h), vs_trawl(ftSY_gig, h, method = "vs_SY"), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

plot(log(h), vs_trawl(trs, h), xlab = "log(observation frequency)", 
     ylab = "Realized variance / observation frequency", 
     main = "Variance signature plot - vs_C",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(log(h), vs_trawl(ftC_exp, h, method = "vs_C"), col = "red", lwd = 2)
lines(log(h), vs_trawl(ftC_gamma, h, method = "vs_C"), col = "blue", lwd = 2)
lines(log(h), vs_trawl(ftC_invGauss, h, method = "vs_C"), col = "orange", lwd = 2)
lines(log(h), vs_trawl(ftC_gig, h, method = "vs_C"), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

## ----ESM8 acf fittings, echo=FALSE, out.width = "49%"--------------------
h <- 1
plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=1 - vs_SY",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftSY_exp, h, method = "vs_SY", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftSY_gamma, h, method = "vs_SY", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftSY_invGauss, h, method = "vs_SY", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftSY_gig, h, method = "vs_SY", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=1 - vs_C",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftC_exp, h, method = "vs_C", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftC_gamma, h, method = "vs_C", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftC_invGauss, h, method = "vs_C", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftC_gig, h, method = "vs_C", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

h <- 0.1
plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=0.1 - vs_SY",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftSY_exp, h, method = "vs_SY", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftSY_gamma, h, method = "vs_SY", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftSY_invGauss, h, method = "vs_SY", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftSY_gig, h, method = "vs_SY", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=0.1 - vs_C",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftC_exp, h, method = "vs_C", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftC_gamma, h, method = "vs_C", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftC_invGauss, h, method = "vs_C", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftC_gig, h, method = "vs_C", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

h <- 0.01
plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=0.01 - vs_SY",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftSY_exp, h, method = "vs_SY", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftSY_gamma, h, method = "vs_SY", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftSY_invGauss, h, method = "vs_SY", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftSY_gig, h, method = "vs_SY", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

plot(1:10, acf_trawl(trs, h, dff = 1, lag_max = 10), xlab = "lag", ylab = "autocorrelation", 
     main = "ACF - h=0.01 - vs_C",
     pch = 16, las = 1, bty = 'n', cex.axis = 1.3, cex.lab = 1.5)
lines(1:10, acf_trawl(ftC_exp, h, method = "vs_C", dff = 1, lag_max = 10), col = "red", lwd = 2)
lines(1:10, acf_trawl(ftC_gamma, h, method = "vs_C", dff = 1, lag_max = 10), col = "blue", lwd = 2)
lines(1:10, acf_trawl(ftC_invGauss, h, method = "vs_C", dff = 1, lag_max = 10), col = "orange", lwd = 2)
lines(1:10, acf_trawl(ftC_gig, h, method = "vs_C", dff = 1, lag_max = 10), col = "purple", lwd = 2)
legend("topright", c("exp", "gamma", "invGauss", "gig"), lwd = rep(2, 4), 
       col = c("red", "blue", "orange", "purple"))

## ----ESM8 nonparametric levy par-----------------------------------------
ftSY_gig$levy_par
ftSY_gig$b

## ----ESM8 cumulants, echo=FALSE, results = 'asis'------------------------
library(knitr)
cums <- matrix(NA, nrow = 2, ncol = 5)
rownames(cums) <- c("second cumulant efficient process", "second cumulant trawl process")
colnames(cums) <- c("vs_SY", "vs_C (exp)", "vs_C (gamma)", "vs_C (invG)", "vs_C (gig)")
cums[1,] <- c(ftSY_gig$b * sum(ftSY_gig$levy_par[, 2]), ftC_exp$xi, ftC_gamma$xi, ftC_invGauss$xi, ftC_gig$xi)
cums[2,] <- c((1 - ftSY_gig$b) * sum(ftSY_gig$levy_par[, 2]), ftC_exp$omega, ftC_gamma$omega, ftC_invGauss$omega, ftC_gig$omega)
kable(cums, caption = "Estimated cumulants")

