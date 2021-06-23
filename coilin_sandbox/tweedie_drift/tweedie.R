##-------------------------------------
## fit a drift term in a Tweedie model
## CM for discussion with Louise
## 22/06/2021
##-------------------------------------
library(icesDatras)
library(ggplot2); theme_set(theme_bw())
library(mgcv)
library(tweedie)

## Norway pout WoRMS ID
aphia_id <- 126444

##------------------------
## HH data - haul details
##------------------------

all_hh <- NULL

s <- "IE-IGFS"

available_years <- getSurveyYearList(s)
for(y in available_years){
    tmp_hh <- getDATRAS(record = "HH", s, y, quarters = 4)
    all_hh <- rbind(all_hh, tmp_hh)
}

## get valid tows
all_hh <- subset(all_hh, HaulVal == "V")

with(all_hh, plot(ShootLong, ShootLat))

## keeping all for now but can restrict areas etc

##-------------------------------
## HL data - biological sampling
##-------------------------------

all_hl <- NULL

for(y in available_years){
    tmp_hl <- getDATRAS(record = "HL", s, y, quarters = 4)
    tmp_hl <- subset(tmp_hl, Valid_Aphia == aphia_id)
    all_hl <- rbind(all_hl, tmp_hl)
}

table(all_hl$LngtCode) ## all in mm

## outlying sub-sampling factors?
## keeping all for now

## convert to mass
## using North Sea W = 0.00450 L^3.148 from
## https://www.fishbase.se/popdyn/LWRelationshipList.php?ID=1023&GenusName=Trisopterus&SpeciesName=esmarkii&fc=183

a <- 0.00450
b <- 3.148

## raised mass per length class per haul
all_hl$mass_g <- with(all_hl, HLNoAtLngt * SubFactor * a * (LngtClass / 10)^b)

## in kgs
all_hl$mass_kg <- all_hl$mass_g / 1e3

## sum to the haul level - ignoring swept area here
haul_hl <- aggregate(mass_kg ~ Year + HaulNo, FUN = sum, data = all_hl)

## bring in the zeros
vars2keep <- c("Year", "HaulNo")

all_dat <- merge(haul_hl, all_hh[, vars2keep], all = TRUE)

all_dat$mass_kg[is.na(all_dat$mass_kg)] <- 0

with(all_dat, plot(jitter(Year), mass_kg))

## fit a models
library(TMB)
compile("tweedie_ts.cpp")
dyn.load(dynlib("tweedie_ts"))

data <- list(y = all_dat$mass_kg,
             idx = all_dat$Year - min(all_dat$Year))
data$Y <- max(data$idx)

## Parameter initial guess
parameters <- list(logmu = rep(log(mean(data$y)), data$Y + 1),
                   logsdmu = log(0.1),
                   d = 0,
                   logphi = log(1),
                   logitp = qlogis(0.5))

## Fit model
obj <- MakeADFun(data,
                 parameters,
                 DLL = "tweedie_ts",
                 random = c("logmu"))
obj$fn()
obj$gr()
opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- sdreport(obj)
srep <- summary(rep)

## profile the drift term
prof <- tmbprofile(obj, "d")
plot(prof)
confint(prof)


p_hat <- 1.0001 + 0.9999 * plogis(srep["logitp", "Estimate"])
phi_hat <- exp(srep["logphi", "Estimate"])
mu_hat <- exp(srep["logmu", "Estimate"])

library(cplm)
fit <- cpglm(mass_kg ~ -1 + factor(Year), data = all_dat, link = "log")
fit2 <- cpglm(mass_kg ~ -1 + factor(Year), data = all_dat, link = "log")

tmp <- summary(fit)$coefficients
rownames(tmp) <- NULL
pred_R <- data.frame(Year = available_years,
                     est = exp(tmp[, "Estimate"]),
                     lwr = exp(tmp[, "Estimate"] - 2 * tmp[, "Std. Error"]),
                     upr = exp(tmp[, "Estimate"] + 2 * tmp[, "Std. Error"]),
                     mod = "R Tweedie")

## add own
tmp <- srep[grep("logmu", rownames(srep)), ]
pred_own <- data.frame(Year = available_years,
                       est = exp(tmp[, "Estimate"]),
                       lwr = exp(tmp[, "Estimate"] - 2 * tmp[, "Std. Error"]),
                       upr = exp(tmp[, "Estimate"] + 2 * tmp[, "Std. Error"]),
                       mod = "Own C++ Tweedie with drift")

all_pred <- rbind(pred_R, pred_own)

pdf("Tweedie_drift_model_example.pdf", height = 7, width = 8)
ggplot(all_pred, aes(x = Year, y = est, colour = mod)) +
    geom_line(size = 1) +
    geom_point() +
    geom_line(aes(y = lwr), lty = 2) +
    geom_line(aes(y = upr), lty = 2) +
    scale_colour_manual("Model", values = c("blue", "slategrey")) +
    theme(legend.position = "bottom") +
    xlab("Year") +
    ylab("Mean CPUE (kg/tow)") +
    ggtitle("Norway pout")
##
plot(prof, xlab = "Drift", ylab = "Profile likelihood", bty = "l", main = "Norway pout drift term profile likelihood")
abline(v = 0, lty = 2)
dev.off()
