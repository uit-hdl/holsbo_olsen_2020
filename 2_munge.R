library(dplyr); library(tidylog)
load("data/raw/metastasis.rda")

################################################################################
############### Fix background info, remove useless obs. #######################
################################################################################
covariables <- co_variables %>% filter(case_control=="case")
summary(covariables)

covariables_superset <- cov_superset
summary(cov_superset)

# Some missing subtypes, set to unknown
covariables$subtype[is.na(covariables$subtype)] <- "Unknown"
covariables_superset$subtype[is.na(covariables_superset$subtype)] <- "Unknown"

# Some missing BMI values; removing these observations
covariables <- covariables %>% filter(!is.na(BMI))
covariables_superset <- covariables_superset %>% filter(!is.na(BMI))

# some negative followup days (by very little); keeping these
covariables$followup[which(covariables$followup < 0)]

# some missing HRT status, removing
covariables <- covariables %>% filter(!is.na(HRT))
covariables_superset <- covariables_superset %>% filter(!is.na(HRT))

# some status "attended" (as opposed screening, interval, clinical). Removing.
covariables <- covariables %>% filter(status != "attended")
covariables_superset <- covariables_superset %>% filter(status != "attended")

# remove clinical
covariables <- covariables %>% filter(status != "outside")
covariables_superset <- covariables_superset %>% filter(status != "outside")

# some missing control in expression values, removing
covariables <- covariables %>% filter(is.element(match_labnr, rownames(gene_exp)))

# define new 0/1 metastasis variable
meta <- rep("Unknown", nrow(covariables))
meta[covariables$metastasis %in% c("A", "B", "C")] <- "Yes"
meta[covariables$metastasis == "0"] <- "No"

covariables$PN[meta == "Unknown"]
## [1] 0 Y 0 1 Y 0 X
# X and Y also unknown

meta[meta=="Unknown" & covariables$PN=="0"] <- "No"
meta[meta=="Unknown" & covariables$PN=="1"] <- "Yes"

covariables <- covariables %>% mutate(metastasis_creg=metastasis, metastasis=meta) %>% 
  filter(metastasis != "Unknown")

meta <- rep("Unknown", nrow(covariables_superset))
meta[covariables_superset$metastasis %in% c("A", "B", "C")] <- "Yes"
meta[covariables_superset$metastasis == "0"] <- "No"

covariables_superset$PN[meta == "Unknown"]
## [1] <NA> X    X    X    X    X    Y    X    X    X 
# X and Y also unknown

meta[meta=="Unknown" & covariables_superset$PN=="0"] <- "No"
meta[meta=="Unknown" & covariables_superset$PN=="1"] <- "Yes"

covariables_superset <- covariables_superset %>% mutate(metastasis_creg=metastasis, metastasis=meta) %>% 
  filter(metastasis != "Unknown")

covariables_controls <- co_variables %>% dplyr::filter(match_labnr %in% covariables$labnr)

summary(covariables_controls)

rm(cov_superset, meta)

################################################################################
########################### Data reduction #####################################
################################################################################
gexpression <- gene_exp[c(covariables$labnr, covariables$match_labnr), ]
rm(gene_exp)

deviations <- plyr::aaply(gexpression, 2, sd)
means <- plyr::aaply(gexpression, 2, mean)
meds <- plyr::aaply(gexpression, 2, median)

plot(density(log2(means/deviations)))
abline(v=log2(20))
axis(1, at=log2(1:100), labels=(1:100))

plot(means, abs(means/deviations), col=ifelse(means/deviations > 20, "black", "grey"), pch=20)
plot(means, deviations, col=ifelse(means/deviations > 20, "black", "grey"), pch=20)

# remove those genes where mean/sd > 20
remove <- means/deviations > 20; sum(remove)
gexpression <- gexpression[, !remove]
dim(gexpression)

rownames(co_variables) <- co_variables$labnr

# used to inform prior
smk <- co_variables[rownames(gexpression), ]$smoking
smkx <- colMeans(gexpression[smk=="Yes", ]) - colMeans(gexpression[smk=="No", ])
sds <- plyr::aaply(gexpression, 2, sd)

hrt <- co_variables[rownames(gexpression), ]$HRT
hrtx <- colMeans(gexpression[hrt=="Yes", ]) - colMeans(gexpression[hrt=="No", ])

ratios_smk <- (smkx/sds)
ratios_hrt <- (hrtx/sds)

pooledx <- c(smkx, hrtx)
pooledrat <- c(ratios_smk, ratios_hrt)


# fold change is log2(expression in case) - log2(expression in control)
fchange <- gexpression[covariables$labnr, ] - gexpression[covariables$match_labnr, ]
dim(fchange)

save(fchange, covariables, covariables_superset, covariables_controls, file="data/metastasis.rda")
