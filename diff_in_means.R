### R script to reproduce the results in:
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
cov.md <- by(ari, ari$id, function(x) covar.smd(y, sd, n, "md", data = x))
ari$md <- unlist(lapply(cov.md, function(x) x$y))
ari$vmd <- unlist(lapply(cov.md, function(x) x$v))
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
modi <- lapply(split(ari, ari$id), function(x)
  dosresmeta(formula = y ~ rcs(dose, knots),
             sd = sd, n = n, covariance = "md", data = x)
)

## Figure 1
pdf("Figure1.pdf", width = 10, height = 6)
par(mfrow = c(2, 3), las = 1, bty = "n")
mapply(function(d, m){
  newdata <- data.frame(dose = seq(0, max(d$dose), length.out = 100))
  with(predict(m, newdata), {
    matplot(newdata$dose, cbind(pred, ci.lb, ci.ub), type = "l", 
            ylim = c(-2, 25), xlim = c(0, 30), lty = c(1, 2, 3), col = "black",
            xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
  })
  with(d[-1, ], errbar(dose, md, md + 1.96*sqrt(vmd), md - 1.96*sqrt(vmd), 
                       add = T, pch = 15, lty = 3, cap = .02))
  title(d$author[1])
}, split(ari, ari$id), modi)
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
results <- do.call("rbind", apply(mvsample, 1, function(y){
  pred <- rcs(newdata$dose, knots) %*% y
  doseEff(p = p, dose = newdata$dose, Ep = pred, trunc = F)
}))
by(results, results$p, function(x) 
  round(quantile(x$ED, c(.025, .975), na.rm = T), 2))


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
knlist <- combn(quantile(ari$dose, c(.1, .25, .5, .75, .9)), 3, simplify = F)[-c(1, 10)]
#do.call(rbind, knlist)
modi_k <- lapply(knlist, function(k)
  dosresmeta(formula = y ~ rcs(dose, k),
             sd = sd, n = n, covariance = "md", data = ari)
)

pdf("Figure3.pdf", width = 15, height = 6)
par(mfrow = c(1, 2))
newdata <- data.frame(dose = seq(0, max(ari$dose), length.out = 500))
par(mar = c(5, 4, 4, 4) + 1.5, las = 1, bty = "n")
with(predict(spl, newdata), {
  plot(newdata$dose, pred, type = "l", col = "white", ylim = c(0, 20), 
       xlim = c(0, 30), xlab = "Aripiprazole (mg/day)", ylab = "Mean Difference")
})
mapply(function(m, k){
  with(predict(m, newdata), lines(newdata$dose, pred, lty = k))
}, modi_k, 1:8)
legend(0, 21, lapply(knlist, function(k) 
  paste("knots = ", paste(k, collapse = ", "))), lty = 1:8, bty = "n")
axis(side = 4, at = edp$Ep[seq(1, 11, 2)], lab = 100*edp$Eprel[seq(1, 11, 2)], pos = 31)
mtext("Relative Efficacy, %                                   ",
      side = 4, line = 1, las = 3, padj =  2)


## 2) Other models:

## Frac Pol
pi <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)
grid <- subset(expand.grid(p1 = pi, p2 = pi), p1 <= p2)
rownames(grid) <- seq(nrow(grid))
shift <- .5
modip <- lapply(
  split(grid, seq(nrow(grid))), function(p)
    dosresmeta(formula = y ~ fracpol(dose, p = p, shift = shift), id = id,
               sd = sd, n = n, covariance = "md", data = ari, proc = "1stage")
)
(pf <- grid[which.min(sapply(modip, AIC)), ])
mod_frac <- dosresmeta(formula = y ~ fracpol(dose, p = pf, shift = shift), id = id,
                       sd = sd, n = n, covariance = "md", data = ari)
## Emax model
Slist <- lapply(cov.md, function(x) x$S)
emaxi <- mapply(function(d, S){
  fitMod(d$dose[-1], d$md[-1], model = "emax", type = "general",
         S = S, placAdj = TRUE, bnds = c(.001, max(d$dose)))
}, split(ari, ari$id), Slist, SIMPLIFY = FALSE)
emax <- mvmeta(do.call("rbind", lapply(emaxi, coef)), lapply(emaxi, vcov))
## Quadratic model
quadr <- dosresmeta(formula = y ~ dose + I(dose^2), id = id,
                    sd = sd, n = n, covariance = "md", data = ari)
## Piecewise with knot = 20.1
k <- 20.1
modi <- lapply(split(ari, ari$id), function(x)
  if (max(x$dose) < k){
    dosresmeta(formula = y ~ dose,
               sd = sd, n = n, covariance = "md", data = x)
  } else {
    dosresmeta(formula = y ~ dose + I((dose - k) * (dose > k)),
               sd = sd, n = n, covariance = "md", data = x)
  }
)
nscoef <- do.call("rbind", lapply(modi, coef))
nscoef[1:2, 2] <- NA
nsvcov <- lapply(modi, vcov)
nsvcov[1:2] <- lapply(nsvcov[1:2], function(s){
  m <- matrix(NA, nrow = 2, ncol = 2)
  m[1, 1] <- s
  m
})
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
coefi <- rbind(do.call("rbind", lapply(modi, coef)), 
               mod_frac$bi, do.call("rbind", lapply(emaxi, coef)),
               quadr$bi, nscoef)
vcovi <- rbind(do.call("rbind", lapply(modi, function(m) c(vcov(m))[-2])),
               do.call("rbind", lapply(mod_frac$Si, function(s) c(s)[-2])),
               do.call("rbind", lapply(lapply(emaxi, vcov), function(s) c(s)[-2])),
               do.call("rbind", lapply(quadr$Si, function(s) c(s)[-2])),
               do.call("rbind", lapply(nsvcov, function(s) c(s)[-2])))
tab4 <- data.frame(rep(1:5, 5), coef = coefi, vi = vcovi)
rownames(tab4)[seq(1, 21, 5)] <- c("Restricted cubic splines", "Fractional Polynomials", 
                                   "Emax", "Quadratic",  "Piecewise linear")
colnames(tab4) <- c("id", "theta1", "theta2", "v11", "v12", "v22")
round(tab4, 4)
