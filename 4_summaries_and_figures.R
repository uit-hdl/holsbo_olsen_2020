library(dplyr)

load("data/metastasis.rda")
load("data/genewise_logistic.rda")

posterior <- extract(fit)

pscreen <- mean(covariables$status == "screening")

# this integrates out detection method by 
# p(theta) = p(screen)*p(theta|screen) + p(interval)*p(theta|interval)
N <- nrow(posterior$alpha[,,1])
n_scr <- round(pscreen*N)
n_int <- round((1-pscreen)*N)

i_scr <- sample(N, n_scr)
i_int <- sample(N, n_int)

alpha2 <- rbind(posterior$alpha[i_scr,,1], posterior$alpha[i_int,,2])
beta2 <- rbind(posterior$beta[i_scr,,1], posterior$beta[i_int,,2])

# converts from logit to probability scale
logistic <- function(x) 1/(1 + exp(-x))

abs_effect2 <- plyr::aaply(1:ncol(fchange), 1, function(i) {
  logistic(alpha2[, i] + beta2[, i]*0.1) - logistic(alpha2[, i])
})
abs_effect2 <- t(abs_effect2)
colnames(abs_effect2) <- colnames(fchange)

# calculate sign probability
sign_prob2 <- plyr::aaply(1:ncol(abs_effect2), 1, function(i) {
  x <- abs_effect2[, i]
  ifelse(median(x) > 0, mean(x > 0), mean(x < 0))
})

ord <- order(sign_prob2, decreasing = T)
errp <- 1 - sign_prob2[ord]
fdp <- plyr::aaply(1:length(errp), 1, function(i) {
  mean(errp[1:i])
})
gnames <- colnames(abs_effect2)[ord]

topgenes2 <- data.frame(gene=gnames, excess_risk=signif(aaply(abs_effect2, 2, median)[ord], 3), `p(increased_risk)`=signif(sign_prob2[ord], 3), FDP=signif(fdp, 2))
#topgenes <- topgenes %>% arrange(desc(excess_risk))
write.csv(topgenes2[1:100, ], file="top100single.csv")
head(topgenes2, n=100)

hist(sign_prob2, nclass=60, border=NA, col="black")

# plot that compares MLE with bayesian estimate
MLE_estimates <- plyr::aaply(fchange, 2, function(x) {
  sset <- covariables$status=="screening"
  dat <- data.frame(metastasis=as.factor(covariables$metastasis), x=(x-mean(x))/sd(x))

  glm_fit <- glm(metastasis~x, data=dat, family=binomial, subset = sset)
  b1 <- coef(glm_fit)[2]
  a1 <- coef(glm_fit)[1]
  p1 <- logistic(a1 + 0.1*b1) - logistic(a1)

  glm_fit <- glm(metastasis~x, data=dat, family=binomial, subset = !sset)
  b2 <- coef(glm_fit)[2]
  a2 <- coef(glm_fit)[1]
  p2 <- logistic(a2 + 0.1*b2) - logistic(a2)

  pscreen*p1 + (1-pscreen)*p2
}, .progress="text")

shrunk_estimates <- plyr::aaply(abs_effect2, 2, mean)
mx <- max(abs(c(MLE_estimates, shrunk_estimates)))

hist(abs(MLE_estimates - shrunk_estimates))
plot(MLE_estimates, shrunk_estimates)
abline(0,1)

pdf(file="shrinkage_size.pdf", width=15, height=4.5)
plot(c(-mx, mx), 0:1, type="n", bty="n", yaxt="n", ylab="", xlab="Point estimate",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
# points(MLE_estimates, rep(0,length(MLE_estimates)))
# points(shrunk_estimates, rep(1,length(MLE_estimates)))
segments(MLE_estimates, rep(0,length(MLE_estimates)), shrunk_estimates, rep(1,length(MLE_estimates)), lwd=2)
text(-.0165, c(.05, .95), labels=c("MLE", "Bayesian"), cex = 1.5)
dev.off()

### Geneset analyses
library(geneset)   # https://github.com/3inar/geneset

c2 <- load_gmt("data/genesets/c2.all.v7.0.symbols.gmt")
c5bp <- load_gmt("data/genesets/c5.bp.v7.0.symbols.gmt")
hallm <- load_gmt("data/genesets/h.all.v7.0.symbols.gmt")

#exctract those genes from gene sets that overlap with our genes
c2 <- gsintersect(c2, colnames(abs_effect2))
c2 <- gsfilter(c2, min=5)
c5bp <- gsintersect(c5bp, colnames(abs_effect2))
c5bp <- gsfilter(c5bp, min=5)
hallm <- gsintersect(hallm, colnames(abs_effect2))
hallm <- gsfilter(hallm, min=5)

bothsets <- gset(c(c2$names, c5bp$names, hallm$names), 
                 c(c2$descriptions, c5bp$descriptions, hallm$descriptions), 
                 c(c2$genesets, c5bp$genesets, hallm$genesets))

sets <- bothsets$genesets
absabs <- sign_prob2
names(absabs) <- colnames(abs_effect2)

# calculates average sign probability in the genes that are in a given gene set
check <- plyr::laply(sets, function(gs) {
  inside <- which(names(absabs) %in% gs)
  inset <- absabs[inside]
  outset <- absabs[-inside]    # also calculates average sign probability outside the set, but this is not used
  
  c(mean(inset), mean(outset))
}, .progress = "text")
# check <- t(check)

pr_sets <- plyr::aaply(check, 1, function(x) x[1] - x[2])

ord <- order(pr_sets, decreasing = T)
geneset_ranking <- data.frame(set_name=bothsets$names[ord], avg_p_sign=signif(check[ord, 1], 3))


write.csv(geneset_ranking[1:50, ], file="top50sets.csv")

# extract top sets
topsets <- bothsets[order(pr_sets, decreasing = T)[1:50]]
represented <- unique(unlist(topsets$genesets))
insets <- plyr::laply(topsets$genesets, function(st) { represented %in% st })
conts <- colSums(insets)
names(conts) <- represented

# find genes that occur in several top sets
in_several_gsets <- sort(names(conts[conts > 1]))
locs <- list()

for (gn in in_several_gsets) {
  for (i in 1:length(topsets)) {
    st <- topsets[i]
    if (gn %in% st$genesets[[1]]) { 
      if (is.null(locs[[gn]])) locs[[gn]] <- as.vector(st$names)
      else locs[[gn]] <- sort(unlist(c(locs[[gn]], list(st$names))))
    }
  }
}

dframe <- data.frame()
for (i in 1:length(locs)) {
  nm <- names(locs)[i]
  gsv <- locs[[i]]
  lt <- length(gsv)
  gsv <- sort(gsv)
  gsv <- paste0(gsv, collapse=", ")
  
  dframe <- rbind(dframe, data.frame(gene=nm, num=lt, sets=gsv))
  
}

dframe <- dframe %>% arrange(desc(num))
write.csv(dframe, file="genes_in_several_sets.csv", row.names = F)

# produce plots of all gene sets
strip <- function(dd, i) {
  current <- dd[, i]
  seqs <- seq(.05, .45, by=.05)
  seqs <- c(seqs, 1-rev(seqs))
  quants <- quantile(current, probs=seqs)
  greys <-  rev(grey.colors(length(seqs)/2))
  greys <- RColorBrewer::brewer.pal(9, "BuGn")
  for (j in 1:(length(seqs)/2)) {
    xx <- c(rep(quants[j], 2), rep(quants[length(quants) + 1 - j], 2))
    yy <- i-1+c(1,1.5, 1.5, 1)
    polygon(x=xx, y=yy, border=NA, density=NA, col=greys[j])
  }
  
  #segments(mean(current), i, mean(current), i+.5, lwd=2)
}

# helper function for plots
genestrips <- function(pathway, dd) {
  gene_names <- as.character(pathway$genesets[[1]])
  dd <- dd[, gene_names]
  print("ACADVL" %in% colnames(dd))
  print(gene_names)
  print(colnames(dd))
  setname <- pathway$names[1]
  omar <- par()$mar
  ls <- max(plyr::aaply(colnames(dd), 1, nchar))
  par(mar=c(5, ls, 4, 2) + 0.1)
  quants <- plyr::aaply(dd, 2, quantile, probs=c(.05, .95))
  plot(range(quants), c(1, ncol(dd)+1), type="n", bty="n", yaxt="n", ylab="", main=setname,
       xlab="excess risk")
  
  dorder <- plyr::aaply(dd, 2, function(x) mean(x > 0))
  dd <- dd[, order(dorder)]
  
  
  for (i in 1:ncol(dd)) strip(dd, i)
  abline(v=0, col=rethinking::col.alpha("black"), lwd=2)
  mtext(side=2, at=1:ncol(dd)+.5, text = colnames(dd), las=2)
  par(mar=omar)
}

for (i in 1:length(topsets)) {
  pathway <- topsets[i]
  size <- length(pathway$genesets[[1]])
  pdf(file=paste0("pathway_plots/", pathway$names, ".pdf"), width=8, height=.2*size + 2)
  genestrips(pathway, abs_effect2)
  dev.off()
}

top100 <- topgenes2[1:100, ]
top100down <- top100[top100$excess_risk < 0, ]
top100up <- top100[top100$excess_risk > 0, ]
nrow(top100down)
nrow(top100up)
in_several <- dframe

other_sets <- gset(gsnames = c("Upregulated out of top 100",
                               "Downregulated out of top100",
                               "Occur in more than one top gene set"),
                   gsdesc = rep("", 3),
                   gslist = list(top100up$gene,
                                 top100down$gene,
                                 in_several$gene))

for (i in 1:length(other_sets)) {
  pw <- other_sets[i]
  size <- length(pw$genesets[[1]])
  pdf(file=paste0(pw$names, ".pdf"), width=8, height=.2*size + 2)
  genestrips(pw, abs_effect2)
  dev.off()
}


hist(abs_effect2[, sample(ncol(abs_effect2), 1)])
library(ggplot2)
library(ggridges)



## prior predictive distribution plots
sd_beta <- exp(rnorm(100000, 0, .2));
prior_mu <- rnorm(100000, 0, 0.1) 
prior_a <- rnorm(100000, -1, 1)
prior_b <- rnorm(100000, prior_mu, sd_beta)

logit <- function(p) log(p) - log(1-p)
logistic <- function(lt) 1/(1 + exp(-lt))

pp_er <- logistic(prior_a + 0.1*prior_b) - logistic(prior_a)

pdf(file="prior_predictive.pdf", width=8*3, height=4.5)
par(mfrow=c(1, 3), cex=1.5)
hist(prior_a, breaks = 100, prob=T, col="darkgrey", border = "white", main=expression("Prior predictive density for "*alpha[g][s]),
     xlab=expression(alpha[g][s]),
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
     abline(v=log(.25) - log(.75), lwd=2)
     legend("topleft", legend="logit(.25)", lty=1, lwd=2, bty = "n")
hist(prior_b, breaks = 100, prob=T, col="darkgrey", border = "white", main=expression("Prior predictive density for "*beta[g][s]),
     xlab=expression(beta[g][s]),
     ylab="",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
hist(pp_er, breaks = 100, prob=T, col="darkgrey", border = "white", main=expression("Prior predictive density for "*rho[g][s]),
     xlab=expression(rho[g][s]),
     ylab="",
     cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
dev.off()

mean(abs(pp_er) > .05)

quantile(pp_er, probs = c(.05, .2, .8, .95))
##           5%         20%         80%         95% 
##   -0.03109530 -0.01386715  0.01429831  0.03228762 
