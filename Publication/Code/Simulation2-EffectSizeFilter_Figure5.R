source("./Publication/MDD.R", local = T)

## 
# II. Data simulation: rnbinom to generate overdispersed count data, 
# data transformation according to Brock et al. (2015)

size = 1000

experiments = data.frame(
  effectsize = runif(n = size),
  theta = runif(size, max = 50,min = 1),
  sampleSize = sample(3:10, size = size, replace = T),
  abund = sample(10:100, size, replace=T)
)
experiments$variance = experiments$theta * experiments$abund

out <- data.frame(
  p = rep(NA,size),
  pMDD = rep(NA,size),
  pMDE = rep(NA,size),
  pCI = rep(NA,size)
)


for (i in 1:size)
{
  N <- experiments$sampleSize[i]
  theta <- experiments$theta[i]
  abund <- experiments$abund[i]
  effectsize <- experiments$effectsize[i]
  control <- rnbinom(n = N, mu = abund, size = abund/(theta-1))
  tmean <- abund * (1 - effectsize)
  treatment <- rnbinom(n = N, mu = tmean, size = tmean/(theta-1))
  transc <- log((2*control)+1)
  transt <- log((2*treatment)+1)
  #tests and calculations on transformed data:
  test <- t.test(transc, transt, alternative = "greater",var.equal = T) #for p-values
  mdd.ln <- MDD(N1 = N, N2 = N, variance1 = (var(transc)+var(transt))/2,
                alpha=0.05, two.sided = F, var.equal = T)
  postpower <- power.t.test(n = N, delta = NULL, sd = sd(c(transc-mean(transc),transt-mean(transt))), 
                             sig.level = 0.05, alternative = "one.sided", type = "two.sample", power = 0.8)
  upperCI <- mean(transc)-mean(transt) + 
    qt(0.05/1, df = 2*N-2, lower.tail = FALSE) * sqrt(var(c(transc-mean(transc),transt-mean(transt)))) * sqrt(2/N)
  #backtransformation
  mdd.abu <- (exp(mean(transc))-exp(mean(transc) - mdd.ln$mdd))/2
  mde.abu <- (exp(mean(transc))-exp(mean(transc) - postpower$delta))/2
  upperCI.abu <- (exp(mean(transc))-exp(mean(transc) - upperCI))/2
  #relation to control mean
  pMDD <- 100 * mdd.abu / ((exp(mean(transc))-1)/2)
  pMDE <- 100 * mde.abu / ((exp(mean(transc))-1)/2)
  pCI <- 100 * upperCI.abu / ((exp(mean(transc))-1)/2)
  out[i,1] <- test$p.value
  out[i,2] <- pMDD
  out[i,3] <- pMDE
  out[i,4] <- pCI
}

#------------------------------------

## PLOT pCI, pMDE pMDD vs. real effectsize:

palette(c("grey60", "black", "red"))

all <- cbind(experiments, out)
nonsig <- all[all$p >= 0.05,]

ylim = max(max(all$pCI), max(all$pMDD), max(all$pMDE)) + 5

op <- par(mfrow=c(2, 3), oma = c(1, 1, 3, 1), mar = c(4, 5, 1, 1), pch = 1, bty = "l", 
          cex.lab = 1.5, las = 1)

# (A-C) All experiments:
dotcol <- as.numeric(all$effectsize * 100 <= all$pCI)
plot(all$pCI ~ all$effectsize, ylim = c(0, ylim),  
     xlab = "", col = dotcol + 1, ylab = "pCI", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 5)
fit1 <- lm(all$pCI ~ all$effectsize)
abline(fit1, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "A", cex = 1.5)

dotcol <- as.numeric(all$effectsize * 100 <= all$pMDE)
plot(all$pMDE ~ all$effectsize, ylim = c(0, ylim),  
     xlab = "", col = dotcol + 1, ylab = "pMDE", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 5)
fit2 <- lm(all$pMDE ~ all$effectsize)
abline(fit2, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "B", cex = 1.5)

dotcol <- as.numeric(all$effectsize * 100 <= all$pMDD)
plot(all$pMDD ~ all$effectsize, ylim = c(0, ylim),  
     xlab = "", col = dotcol + 1, ylab = "pMDD", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 5)
fit3 <- lm(all$pMDD ~ all$effectsize)
abline(fit3, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "C", cex = 1.5)


# (D-F) only non-significant experiments:
dotcol <- as.numeric(nonsig$effectsize * 100 <= nonsig$pCI)
plot(nonsig$pCI ~ nonsig$effectsize, ylim = c(0, ylim), xlim = c(0, 1),
     xlab = "", col = dotcol + 1, ylab = "pCI", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 2)
fit1 <- lm(nonsig$pCI ~ nonsig$effectsize)
abline(fit1, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "D", cex = 1.5)

dotcol <- as.numeric(nonsig$effectsize * 100 <= nonsig$pMDE)
plot(nonsig$pMDE ~ nonsig$effectsize, ylim = c(0, ylim),  xlim = c(0, 1),
     xlab = "Real effect size", col = dotcol + 1, ylab = "pMDE", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 2)
fit2 <- lm(nonsig$pMDE ~ nonsig$effectsize)
abline(fit2, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "E", cex = 1.5)

dotcol <- as.numeric(nonsig$effectsize * 100 <= nonsig$pMDD)
plot(nonsig$pMDD ~ nonsig$effectsize, ylim = c(0, ylim),  xlim = c(0, 1),
     xlab = "", col = dotcol + 1, ylab = "pMDD", pch = dotcol, 
     xaxt = "n", yaxt = "n", cex = 0.8)
axis(2, at = c(0, 50, 100, 150, 200))
axis(1, at = c(0, 0.5, 1))
abline(0, 100, col = "grey40", lwd = 2, lty = 2)
fit3 <- lm(nonsig$pMDD ~ nonsig$effectsize)
abline(fit3, col = "darkred", lwd = 2, lty = "solid")
text(0.1, ylim, "F", cex = 1.5)

par(op)