### R script (avoiding use of lapply) to reproduce the results in:
###   Crippa A, Orsini N. Dose-response meta-analysis of differences in means. 
###   BMC medical research methodology. 2016 Aug 2;16(1):91.

## Packages required
#devtools::install_github("alecri/dosresmeta")
library(dosresmeta)
library(rms)
library(mvtnorm)
library(DoseFinding)
library(aod)

## Loading data
data("ari")
## description of the data and variable
help("ari")

## Auxiliary function to estimate target doses (to be implemented in dosresmeta pkg)
doseEff <- function(p, dose, Ep, trunc = FALSE){
  max <- max(Ep)
  EDmax <- dose[which.min(abs(Ep - max))]
  if (trunc == TRUE){
    if (EDmax == max(dose)) return(data.frame(p = NA, ED = NA, Ep = NA))      
  }
  ED <- apply(matrix(p), 1, function(x)
    dose[which.min(abs(Ep[dose < EDmax] - x * max(Ep)))])
  
  data.frame(p, ED, Ep = p * max(Ep))
}


## -----------------------------------------------------------------------------
## Illustrative example using the trial by Cutler et al.
cutler <- subset(ari, author == "Cutler 2006")
cutler

## Calculate mean differences, variances and covariance matrix
covar.smd(y = y, sd = sd, n = n, measure = "md", data = cutler)

## Dose-response model using restricted cubic spline
knots <- quantile(ari$dose, c(.1, .5, .9))
spl_ex <- dosresmeta(formula = y ~ rcs(dose, knots), sd = sd, n = n, 
                     covariance = "md", data = cutler)

## Summary of the model, print coefficients and covariance matrix
summary(spl_ex)
round(coef(spl_ex), 3)
round(vcov(spl_ex), 2)



## -----------------------------------------------------------------------------
## Main analysis

## Obtaining mean differences, variances, and (co)varinace matrices for the all the studies
ari$md <- c()
ari$vmd <- c()
Slist <- list()
for (i in 1:5){
  cov.md <- covar.smd(y, sd, n, "md", data = subset(ari, id == i))
  ari$md[ari$id == i] <- cov.md$y
  ari$vmd[ari$id == i] <- cov.md$v
  Slist[[i]] <- cov.md$S
}
## Data presented in Table 3
ari

## Dose-response meta-analysis
spl <- dosresmeta(formula = y ~ rcs(dose, knots), id = id, sd = sd, n = n, 
                  covariance = "md", data = ari)
## Summary of the model, print coefficients and covariance matrix
summary(spl)
round(coef(spl), 3)
round(vcov(spl), 2)

## Study-specific models (for Figure 1)
modi <- list()
for (i in 1:5){
  modi[[i]] <- dosresmeta(formula = y ~ rcs(dose, knots), sd = sd, n = n, 
                          covariance = "md", data = subset(ari, id == i))
}

## Figure 1
pdf("Figure1.pdf", width = 10, height = 6)
par(mfrow = c(2, 3), las = 1, bty = "n")
for (i in 1:5){
  newdata <- data.frame(dose = seq(0, max(ari$dose[ari$id == i]), length.out = 100))
  with(predict(modi[[i]], newdata), {
    matplot(newdata$dose, cbind(pred, ci.lb, ci.ub), type = "l", 
            ylim = c(-2, 25), xlim = c(0, 30), lty = c(1, 2, 3), col = "black",
            xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
  })
  with(subset(ari, id == i & dose != 0), 
              errbar(dose, md, md + 1.96*sqrt(vmd), md - 1.96*sqrt(vmd),
                     add = T, pch = 15, lty = 3, cap = .02))
  title(ari$author[ari$id == i & ari$dose == 0])
}
dev.off()

## Tabular prediction
newdata <- data.frame(dose = seq(0, 30, by = 5))
pred_md <- predict(spl, newdata = newdata, xref = 0, expo = FALSE)
round(pred_md, 2)

## Target doses with 'confidence interval'
p <- c(.5, .8, 1)
newdata <- data.frame(dose = seq(0, max(ari$dose), length.out = 5000))
edp <- with(predict(spl, newdata), 
            doseEff(p = p, dose = newdata$dose, Ep = pred, trunc = FALSE))
round(edp, 2)

set.seed(1234)
mvsample <- rmvnorm(10000, mean = coef(spl), vcov(spl))
results <- array(NA, dim = c(10000, 3, 3))
for (i in 1:10000){
  pred <- rcs(newdata$dose, knots) %*% mvsample[i, ]
  results[i, 1:3, 1:3] <- as.matrix(doseEff(p = p, dose = newdata$dose, Ep = pred, trunc = F))
}
round(quantile(results[, 1, 2], c(.025, .975), na.rm = T), 2)
round(quantile(results[, 2, 2], c(.025, .975), na.rm = T), 2)
round(quantile(results[, 3, 2], c(.025, .975), na.rm = T), 2)


## Figure 2
p <- seq(0, 1, .1)
newdata <- data.frame(dose = seq(0, max(ari$dose), length.out = 5000))
edp <- with(predict(spl, newdata),
            doseEff(p = p, dose = newdata$dose, Ep = pred, trunc = FALSE)
)
edp <- cbind(edp, Eprel = edp$Ep/max(edp$Ep))

pdf("Figure2.pdf", width = 8, height = 6)
par(mar = c(5, 4, 4, 4) + 1.5, mfrow = c(1, 1), las = 1, bty = "n")
with(predict(spl, newdata, xref = 0, expo = FALSE), {
  matplot(newdata$dose, cbind(pred, ci.lb, ci.ub), type = "l", 
          ylim = c(-2, 25), xlim = c(0, 30), lty = c(1, 2, 3), col = "black",
          xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
})
axis(side = 4, at = edp$Ep[seq(1, 11, 2)], lab = 100*edp$Eprel[seq(1, 11, 2)],
     pos = 32)
mtext("Relative Efficacy, %                      ", 
      side = 4, line = 2, las = 3, padj =  2.5)
w <- 1/ari$vmd[ari$vmd != 0]
with(subset(ari, dose!= 0), points(dose, md, pch = 1, cex = 2*w/max(w)))
dev.off()


## -----------------------------------------------------------------------------
## Sensitivity analysis

## 1) Location of knots
knmat <- combn(quantile(ari$dose, c(.1, .25, .5, .75, .9)), 3, simplify = T)[, -c(1, 10)]
#t(knmat)
modi_k <- list()
for (i in 1:8){
  modi_k[[i]] <- dosresmeta(formula = y ~ rcs(dose, knmat[, i]), sd = sd, n = n, 
                          covariance = "md", data = ari)
}

pdf("Figure3.pdf", width = 15, height = 6)
par(mfrow = c(1, 2))
newdata <- data.frame(dose = seq(0, max(ari$dose), length.out = 500))
par(mar = c(5, 4, 4, 4) + 1.5, las = 1, bty = "n")
with(predict(spl, newdata), {
  plot(newdata$dose, pred, type = "l", col = "white", ylim = c(0, 20), 
       xlim = c(0, 30), xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
})
for (i in 1:8){
  with(predict(modi_k[[i]], newdata), lines(newdata$dose, pred, lty = i))
}
legend(0, 21, apply(knmat, 2, function(k) paste("knots = ", paste(k, collapse = ", "))), 
       lty = 1:8, bty = "n")
axis(side = 4, at = edp$Ep[seq(1, 11, 2)], lab = 100*edp$Eprel[seq(1, 11, 2)], pos = 31)
mtext("Relative Efficacy, %                                   ",
      side = 4, line = 1, las = 3, padj =  2)


## 2) Other models:

## Frac Pol
pi <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
grid <- subset(expand.grid(p1 = pi, p2 = pi), p1 <= p2)
rownames(grid) <- seq(nrow(grid))
shift <- .5

modip <- list()
aics <- vector()
for (i in 1:nrow(grid)){
  modip[[i]] <- dosresmeta(formula = y ~ fracpol(dose, p = grid[i, ], shift = shift), id = id,
                           sd = sd, n = n, covariance = "md", data = ari, proc = "1stage")
  aics[i] <- AIC(modip[[i]])
}
(pf <- grid[which.min(aics), ])
mod_frac <- dosresmeta(formula = y ~ fracpol(dose, p = pf, shift = shift), id = id,
                       sd = sd, n = n, covariance = "md", data = ari)
## Emax model
emaxi <- list()
emax_coef <- matrix(nrow = 5, ncol = 2)
emax_vcov <- list()
for (i in 1:5){
  emaxi[[i]] <- with(subset(ari, id == i & dose != 0),
                     fitMod(dose, md, model = "emax", type = "general",
                            S = Slist[[i]], placAdj = TRUE, bnds = c(.001, max(dose))))
  emax_coef[i, ] <- coef(emaxi[[i]])
  emax_vcov[[i]] <- vcov(emaxi[[i]])
}
emax <- mvmeta(emax_coef, emax_vcov)
## Quadratic model
quadr <- dosresmeta(formula = y ~ dose + I(dose^2), id = id,
                    sd = sd, n = n, covariance = "md", data = ari)
## Piecewise with knot = 20.1
k <- 20.1
modi_ns <- list()
nscoef <- matrix(nrow = 5, ncol = 2)
nsvcov <- list()
for (i in 1:5){
  ind <- 1 + (max(ari$dose[ari$id == i]) >= k)
  modi_ns[[i]] <-  if (ind == 1){
    dosresmeta(formula = y ~ dose,
               sd = sd, n = n, covariance = "md", data = subset(ari, id == i))
  } else {
    dosresmeta(formula = y ~ dose + I((dose - k) * (dose > k)),
               sd = sd, n = n, covariance = "md", data = subset(ari, id == i))
  }
  nscoef[i, 1:ind] <- coef(modi_ns[[i]])
  nsvcov[[i]] <- matrix(NA, nrow = 2, ncol = 2)
  nsvcov[[i]][1:ind, 1:ind] <- vcov(modi_ns[[i]])
}
nspl <- mvmeta(nscoef ~ 1, nsvcov)
wald.test(vcov(nspl), coef(nspl), L = rbind(c(1, 1)))

## Predictions
newdata <- data.frame(dose = seq(0, max(ari$dose), length.out = 500))
predEmax <- emax(newdata$dose, 0, coef(emax)[1], coef(emax)[2])
predNs <- cbind(newdata$dose, (newdata$dose - k) * (newdata$dose > k)) %*% coef(nspl)

## Graphical comparison
par(mar = c(5, 4, 4, 4) + 1.5, las = 1, bty = "n")
with(predict(spl, newdata), {
  plot(newdata$dose, pred, type = "l", ylim = c(0, 20), xlim = c(0, 30),
       xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
})
with(predict(mod_frac, newdata), {
  lines(newdata$dose, pred, lty = 3)
})
with(predict(quadr, newdata), lines(newdata$dose, pred, lty = 2))
lines(newdata$dose, predEmax, lty = 3)
lines(newdata$dose, predNs, lty = 4)
w <- 1/ari$vmd[ari$vmd != 0]
with(subset(ari, dose!= 0), points(dose, md, pch = 1, cex = 2*w/max(w)))
legend(0, 21, c("Restricted cubic spline", "Fractional Polynomial", "Quadratic", 
                "Emax"),
       lty = 1:4, bty = "n")
axis(side = 4, at = edp$Ep[seq(1, 11, 2)], lab = 100*edp$Eprel[seq(1, 11, 2)],
     pos = 31)
mtext("Relative Efficacy, %                                   ", 
      side = 4, line = 1, las = 3, padj =  2)
dev.off()


## Table 4
## Note that there is a typo in the published article:
## the coefficients for the spline models are misreported: they are the coefficients
## of the piecewise linear models
tab <- matrix(nrow = 25, ncol = 5)
for (i in 1:5){
  tab[i, 1:2] <- coef(modi[[i]])
  tab[i, 3:5] <- c(vcov(modi[[i]]))[-2]
  tab[i + 5, 3:5] <- c(mod_frac$Si[[i]])[-2]
  tab[i + 10, 3:5] <- c(emax_vcov[[i]])[-2]
  tab[i + 15, 3:5] <- c(quadr$Si[[i]])[-2]
  tab[i + 20, 3:5] <- c(nsvcov[[i]])[-2]
}
tab[6:10, 1:2] <- mod_frac$bi
tab[11:15, 1:2] <- emax_coef
tab[16:20, 1:2] <- quadr$bi
tab[21:25, 1:2] <- nscoef
tab4 <- data.frame(rep(1:5, 5), tab)
rownames(tab4)[seq(1, 21, 5)] <- c("Restricted cubic splines", "Fractional Polynomials", 
                                   "Emax", "Quadratic",  "Piecewise linear")
colnames(tab4) <- c("id", "theta1", "theta2", "v11", "v12", "v22")
round(tab4, 4)
