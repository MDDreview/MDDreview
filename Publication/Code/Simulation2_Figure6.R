## based on simulations in Simulation2-EffectSizeFilter_Figure5.R: nonsig dataset
## >> run Simulation2-EffectSizeFilter_Figure5.R first

##-----
# Effect sizes passing regulatory thresholds:
# two filters: 1. alpha = 0.05, 2. threshold ranging from 30 to 100

outeffect <- data.frame(
  indicator = NA,
  threshold = NA,
  effectsize = NA
)

ind <- list("pCI", "pMDE", "pMDD")

for (i in ind){
  for (j in 30:100){
    sub <- nonsig[nonsig[, i] <= j, ]
    temp <- data.frame(
      indicator = rep(i, nrow(sub)),
      threshold = rep(j, nrow(sub)),
      effectsize = sub$effectsize
    )
    outeffect <- rbind(outeffect, temp)
  }
}
outeffect <- na.omit(outeffect)

pCIpass <- outeffect[outeffect$indicator == "pCI", ] 
pMDEpass <- outeffect[outeffect$indicator == "pMDE", ] 
pMDDpass <- outeffect[outeffect$indicator == "pMDD", ] 


library(dplyr)
summ.pCI <- pCIpass %>%
  group_by(threshold) %>%
  summarise(min = min(effectsize), q25 = quantile(effectsize, 0.25), median = median(effectsize), 
            q75 = quantile(effectsize, 0.75), max = max(effectsize), N = n(), 
            mean = mean(effectsize), sd = sd(effectsize))
summ.pMDE <- pMDEpass %>%
  group_by(threshold) %>%
  summarise(min = min(effectsize), q25 = quantile(effectsize, 0.25), median = median(effectsize), 
            q75 = quantile(effectsize, 0.75), max = max(effectsize), N = n(),
            mean = mean(effectsize), sd = sd(effectsize))
summ.pMDD <- pMDDpass %>%
  group_by(threshold) %>%
  summarise(min = min(effectsize), q25 = quantile(effectsize, 0.25), median = median(effectsize), 
            q75 = quantile(effectsize, 0.75), max = max(effectsize), N = n(),
            mean = mean(effectsize), sd = sd(effectsize))


##########-------------------

pdf("PLOTS/Figure4_effectfilter.pdf", height = 4, width = 7)

## plot max, quartiles, min
layout(matrix(c(1,2,3), 1, 3, byrow = TRUE), width = c(1.65,1,1))

## plot1: pCI (9)
op <- par(mar = c(3,0,1,1), oma = c(0.5,8,0,0), cex = 0.8)

plot(summ.pCI$max ~ summ.pCI$threshold, type ="n", ylim = c(0,1.03), xlim = rev(c(30,100)),
     las = 2, xlab = "", ylab = "", las = 1, bty ="l", yaxt = "n", xaxt = "n")
axis(side = 1)
axis(side = 2, las = 2, labels = T, at = seq(0, 1, by =0.1), line = 4, outer = T)
mtext("Real effect size", side = 2, line = 7)
xcor = c(summ.pCI$threshold, rev(summ.pCI$threshold))
ycor = c(smooth.spline(summ.pCI$max ~ summ.pCI$threshold)$y, rev(smooth.spline(summ.pCI$min ~ summ.pCI$threshold)$y))
polygon(xcor, ycor, border =NA, col = adjustcolor("slategray1", alpha.f = 0.2))
ycor = c(smooth.spline(summ.pCI$q75 ~ summ.pCI$threshold)$y, rev(smooth.spline(summ.pCI$q25 ~ summ.pCI$threshold)$y))
polygon(xcor, ycor, border =NA, col = "slategray1")
lines(smooth.spline(summ.pCI$max ~ summ.pCI$threshold), lty = 3, col = "darkblue")
lines(smooth.spline(summ.pCI$q75 ~ summ.pCI$threshold), lty = 2, col = "darkblue")
lines(smooth.spline(summ.pCI$median ~ summ.pCI$threshold), lty = 1, lwd =2, col = "darkblue")
lines(smooth.spline(summ.pCI$q25 ~ summ.pCI$threshold), lty = 2, col = "darkblue")
lines(smooth.spline(summ.pCI$min ~ summ.pCI$threshold), lty = 3, col = "darkblue")
text(100, 1.03, "pCI filter", col = "darkblue", pos = 4)

# add effectsizes passing p-value (> 0.05): min, 1st quartile, median, 3rd quartile and max
abline(h = max(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = min(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = median(nonsig$effectsize), lty = 1, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.25), lty = 2, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.75), lty = 2, col = "grey60")

mtext(c("p-value filter","max", "75%", "median", "25%", "min"), 
      side = 2, line = 0.5, las = 2, cex = 0.8, adj = 1,
      at = c(1.03,max(nonsig$effectsize)-0.01, quantile(nonsig$effectsize, 0.75), median(nonsig$effectsize),
             quantile(nonsig$effectsize, 0.25), min(nonsig$effectsize)))


## plot2: pMDE (8)
op <- par(mar = c(3,0,1,1))

plot(summ.pMDE$max ~ summ.pMDE$threshold, type ="n", ylim = c(0,1.03), xlim = rev(c(30,100)),
     las = 2, xlab = "", ylab = "", las = 1, bty ="l",
     yaxt = "n", xaxt = "n")
axis(side = 1)
mtext("Threshold level", side = 1, line = 2.5)
#axis(side = 4, las = 2, pos = c(30, 1), labels = F, at = seq(0, 1, by =0.1))
ycor = c(smooth.spline(summ.pMDE$max ~ summ.pMDE$threshold)$y, rev(smooth.spline(summ.pMDE$min ~ summ.pMDE$threshold)$y))
polygon(xcor, ycor, border =NA, col = adjustcolor("darkseagreen2", alpha.f = 0.2))
ycor = c(smooth.spline(summ.pMDE$q75 ~ summ.pMDE$threshold)$y, rev(smooth.spline(summ.pMDE$q25 ~ summ.pMDE$threshold)$y))
polygon(xcor, ycor, border =NA, col = "darkseagreen2")
lines(smooth.spline(summ.pMDE$max ~ summ.pMDE$threshold), lty = 3, col = "darkgreen")
lines(smooth.spline(summ.pMDE$q75 ~ summ.pMDE$threshold), lty = 2, col = "darkgreen")
lines(smooth.spline(summ.pMDE$median ~ summ.pMDE$threshold), lty = 1, lwd =2, col = "darkgreen")
lines(smooth.spline(summ.pMDE$q25 ~ summ.pMDE$threshold), lty = 2, col = "darkgreen")
lines(smooth.spline(summ.pMDE$min ~ summ.pMDE$threshold), lty = 3, col = "darkgreen")
text(100, 1.03, "pMDE filter", col = "darkgreen", pos = 4)

# add effectsizes passing p-value (> 0.05): min, 1st quartile, median, 3rd quartile and max
abline(h = max(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = min(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = median(nonsig$effectsize), lty = 1, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.25), lty = 2, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.75), lty = 2, col = "grey60")

## plot3: pMDD (7)
plot(summ.pMDD$max ~ summ.pMDD$threshold, type ="n", ylim = c(0,1.03), xlim = rev(c(30,100)),
     las=2, xlab = "", 
     ylab = "", las = 1, bty ="l", yaxt = "n", xaxt = "n")
axis(side = 1)
#axis(side = 4, las = 2, pos = c(30, 1), labels = T, at = seq(0, 1, by =0.1))
ycor = c(smooth.spline(summ.pMDD$max ~ summ.pMDD$threshold)$y, rev(smooth.spline(summ.pMDD$min ~ summ.pMDD$threshold)$y))
polygon(xcor, ycor, border =NA, col = adjustcolor( "mistyrose2", alpha.f = 0.2))
ycor = c(smooth.spline(summ.pMDD$q75 ~ summ.pMDD$threshold)$y, rev(smooth.spline(summ.pMDD$q25 ~ summ.pMDD$threshold)$y))
polygon(xcor, ycor, border =NA, col = "mistyrose2")
lines(smooth.spline(summ.pMDD$max ~ summ.pMDD$threshold), lty = 3, col = "darkred")
lines(smooth.spline(summ.pMDD$q75 ~ summ.pMDD$threshold), lty = 2, col = "darkred")
lines(smooth.spline(summ.pMDD$median ~ summ.pMDD$threshold), lty = 1, lwd =2, col = "darkred")
lines(smooth.spline(summ.pMDD$q25 ~ summ.pMDD$threshold), lty = 2, col = "darkred")
lines(smooth.spline(summ.pMDD$min ~ summ.pMDD$threshold), lty = 3, col = "darkred")
text(100, 1.03, "pMDD filter", col = "darkred", pos = 4)

# add effectsizes passing p-value (> 0.05): min, 1st quartile, median, 3rd quartile and max
abline(h = max(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = min(nonsig$effectsize), lty = 3, col = "grey60")
abline(h = median(nonsig$effectsize), lty = 1, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.25), lty = 2, col = "grey60")
abline(h = quantile(nonsig$effectsize, 0.75), lty = 2, col = "grey60")

par(op)

dev.off()

########--------------------------------------------------------------------
  ## plot: prop of nonsig experiments passing the secondary filter:
# plot(summ.pCI$N/nrow(nonsig) ~ summ.pCI$threshold, type = "n", bty ="l", ylim = c(0,1))
# lines(smooth.spline(summ.pCI$N/nrow(nonsig) ~ summ.pCI$threshold), col = "darkblue", lwd = 2)
# lines(smooth.spline(summ.pMDE$N/nrow(nonsig) ~ summ.pMDE$threshold), col = "darkgreen", lwd = 2)
# lines(smooth.spline(summ.pMDD$N/nrow(nonsig) ~ summ.pMDD$threshold), col = "darkred", lwd = 2)
# 
# plot(summ.pCI$N ~ rev(summ.pCI$threshold), type = "n", bty ="l", ylim = c(0,max(summ.pCI$N)),
#      xlim = rev(c(30,100)),
#      ylab = "Sample size", xlab = "Threshold", las = 1)
# lines(smooth.spline(summ.pCI$N ~ summ.pCI$threshold), col = "darkblue", lwd = 2)
# lines(smooth.spline(summ.pMDE$N ~ summ.pMDE$threshold), col = "darkgreen", lwd = 2)
# lines(smooth.spline(summ.pMDD$N ~ summ.pMDD$threshold), col = "darkred", lwd = 2)
# 
# 
# ##################---------------------------------
# ## plot mean +-SD
# op <- par(las = 1)
# 
# plot(summ.pCI$mean ~ summ.pCI$threshold, type = "n", bty = "l", 
#      ylab = "Real effect size", xlab = "Threshold",
#      ylim = c(0,1), xlim = rev(c(30,100)))
# lines(smooth.spline(summ.pCI$mean ~ summ.pCI$threshold), col = "darkblue", lwd = 2)
# lines(smooth.spline(summ.pCI$mean + summ.pCI$sd ~ summ.pCI$threshold), col = "darkblue", lty = 2)
# lines(smooth.spline(summ.pCI$mean - summ.pCI$sd ~ summ.pCI$threshold), col = "darkblue", lty = 2)
# 
# lines(smooth.spline(summ.pMDE$mean ~ summ.pMDE$threshold), col = "darkgreen", lwd = 2)
# lines(smooth.spline(summ.pMDE$mean + summ.pMDE$sd ~ summ.pMDE$threshold), col = "darkgreen", lty = 2)
# lines(smooth.spline(summ.pMDE$mean - summ.pMDE$sd ~ summ.pMDE$threshold), col = "darkgreen", lty = 2)
# 
# lines(smooth.spline(summ.pMDD$mean ~ summ.pMDD$threshold), col = "darkred", lwd = 2)
# lines(smooth.spline(summ.pMDD$mean + summ.pMDD$sd ~ summ.pMDD$threshold), col = "darkred", lty = 2)
# lines(smooth.spline(summ.pMDD$mean - summ.pMDD$sd ~ summ.pMDD$threshold), col = "darkred", lty = 2)
# 
# abline(h = mean(nonsig$effectsize), lty = 1, col = "grey60")
# abline(h = mean(nonsig$effectsize) + sd(nonsig$effectsize), lty = 2, col = "grey60")
# abline(h = mean(nonsig$effectsize) - sd(nonsig$effectsize), lty = 2, col = "grey60")
# 
# par(op)
# 
