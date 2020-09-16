library(nowac)     # internal nowac package that extracts data from questionnaires
library(stringr)
library(dplyr)

data("samples")
head(samples)

prospective <- samples %>% filter(str_detect(Project, "breast cancer run"))
  
######################### load EXPRESSION DATA ######################################
load("censored")
ex_run_1 <- list(expression=exprs, background=background10, lab_info=labInfo, negative_controls=negCtrl, lumi_batch=data)

load("censored")
ex_run_2 <- list(expression=exprs, background=background10, lab_info=labInfo, negative_controls=negCtrl, lumi_batch=data)

load("censored")
ex_run_3 <- list(expression=exprs, background=background10, lab_info=labInfo, negative_controls=negCtrl, lumi_batch=data)

save(ex_run_1, ex_run_2, ex_run_3, file="data/raw/prospective_expression.rda")
rm(exprs, background10, labInfo, negCtrl, data, samples)

# need expression, dates
tmp_1 <- ex_run_1$background %>% dplyr::select(labnr = labnr, lpnr = LPNR, 
                                        diagnosis_date = dxdate, 
                                        sample_date = BPROVEDATO) %>%
                                        mutate(followup = diagnosis_date - sample_date)
tmp_2 <- ex_run_2$background %>% dplyr::select(labnr = labnr, lpnr = LPNR, 
                                        diagnosis_date = dxdate, 
                                        sample_date = BPROVEDATO) %>%
                                        mutate(followup = diagnosis_date - sample_date)
tmp_3 <- ex_run_3$background %>% dplyr::select(labnr = labnr, lpnr = LPNR, 
                                        diagnosis_date = dxdate, 
                                        sample_date = BPROVEDATO) %>%
                                        mutate(followup = diagnosis_date - sample_date)

# only 9 observations outside of run 1 have followup within 2 years,
# we leave these out to avoid thinking too hard about batch effects
sum(na.omit(tmp_1$followup <= 2*365)) # 250 case
sum(na.omit(tmp_2$followup <= 2*365)) # 8
sum(na.omit(tmp_3$followup <= 2*365)) # 1

followup <- rbind(tmp_1,tmp_2,tmp_3)
rm(tmp_1, tmp_2, tmp_3)

rm(ex_run_2, ex_run_3)


########################  load CANCER REGISTRY DATA #################################
cancer_registry <- read.csv("censored",  sep = ",") %>% filter(IO_NR %in% prospective$LPNR)
save(cancer_registry, file="data/raw/cancer_registry.rda")

rownames(cancer_registry) <- cancer_registry$IO_NR
c_variables <- cancer_registry[as.character(prospective$LPNR), ] %>% dplyr::select(IO_NR, DS1, ICD101, DIAGNOSEDA, P_T1, P_N1, P_M1, ER1, PR1, HER21, INSITU1, METASTASE1)

# RECODE ER
ER <- c_variables$ER1
ER[ER==1] <- -1  # neg
ER[ER==98] <- -3 # missing
ER[!ER %in% c(-1, -3, NA)] <- -2 # pos
ER <- factor(-ER, levels=1:3, labels=c("Negative", "Positive", "Missing"))

# RECODE PR
PR <- c_variables$PR1
PR[PR %in% c(1, 12645)] <- -1 # neg
PR[PR==98] <- -3 # missing
PR[!PR %in% c(-1, -3, NA)] <- -2 # pos
PR <- factor(-PR, levels=1:3, labels=c("Negative", "Positive", "Missing"))

# RECODE HER2
HER2 <- c_variables$HER21
HER2 <- factor(HER2, levels=1:3, labels=c("Missing", "Negative", "Positive"))

# molecular subtypes
luminal_a <- (ER == "Positive" | PR == "Positive") & HER2 == "Negative"
luminal_b <- (ER == "Positive" | PR == "Positive") & HER2 == "Positive"
triple_negative <- (ER == "Negative" & PR == "Negative") & HER2 == "Negative"
HER2_positive <- (ER == "Negative" & PR == "Negative") & HER2 == "Positive"
unknown <- (ER == "Missing" | PR == "Missing") | HER2 == "Missing" 

summary(luminal_a)
summary(luminal_b)
summary(triple_negative)
summary(HER2_positive)
summary(unknown)

subtype <- rep(NA, length(prospective))
subtype[luminal_a] <- "Luminal A"
subtype[luminal_b] <- "Luminal B"
subtype[triple_negative] <- "Triple Negative"
subtype[HER2_positive] <- "HER2 Positive"
subtype[unknown] <- "Unknown"

subtype[is.na(subtype) & !(is.na(ER) & is.na(PR) & is.na(HER2))] <- "Unknown"

subtype <- as.factor(subtype)
summary(subtype)

length(which(luminal_b + luminal_a + triple_negative + HER2_positive > 1))

cancer_registry <- c_variables %>% dplyr::select(lpnr=IO_NR, PT=P_T1, PN=P_N1, PM=P_M1, 
                                          metastasis=METASTASE1, ICD10=ICD101, diagnostic_certainty=DS1) %>% 
  mutate(HER2=HER2, ER=ER, PR=PR, subtype=subtype)

rm(c_variables, ER, PR, HER2, HER2_positive, luminal_a, luminal_b, subtype, triple_negative, unknown)

##  detection method
blood <- read.table("censored", na.strings="",sep="\t",header=TRUE) # "status" in this one
save(blood, file="data/raw/extra_variables.rda")

questionnaire <- blood %>% dplyr::select(labnr=labnr, case_control=Case_ctrl, 
                                  smoking=BROYK,  HRT=BHRT, status=status, 
                                  age=age_at_blood, age_diagnosis=age_at_dx) %>% 
  mutate(smoking=factor(smoking, levels=0:1, labels=c("Yes", "No")), 
         HRT=factor(HRT, levels=0:1, labels=c("Yes", "No"))) 

rm(blood)

# PARITY
parity <- read.csv("censored") %>% filter(labnr %in% prospective$LabNumber)
rownames(parity) <- parity$labnr
parity <- parity[as.character(prospective$LabNumber), ] %>% transmute(labnr, project, parity=antallbarn)
save(parity, file="data/raw/parity.rda")
  
  
#########################  load 2page QUESTIONNAIRE DATA ###################################
two_page <- nowac::getTwoPage(prospective)
rownames(two_page) <- two_page$LPNR
q_variables <- two_page[as.character(prospective$LPNR), ] %>% dplyr::select(LPNR, Case_ctrl, match_labnr, FAAR, BPROVEDATO, BROYK, BVEKT, BHOYDE, BHRT)

# BMI  kg/m^2
q_variables <- q_variables %>% transmute(BMI=BVEKT/(BHOYDE/100)^2, lpnr=LPNR, match_labnr, case_control=Case_ctrl, height=BHOYDE, weight=BVEKT)

save(two_page, file="data/raw/two_page_questionnaire.rda")
rm(two_page)


summary(q_variables)
summary(followup)
summary(questionnaire)
summary(cancer_registry) 
summary(parity)

predictors <- q_variables %>% left_join(followup) %>% 
left_join(questionnaire) %>% 
  left_join(cancer_registry) %>% 
  left_join(parity) %>% mutate(case_control=as.factor(case_control), 
                             followup=as.numeric(followup),
                             labnr=as.character(labnr), 
                             match_labnr=as.character(match_labnr)) 

rm(q_variables, followup, questionnaire, cancer_registry, parity, prospective)

# Only breast cancers and healthy controls
predictors <- predictors %>% filter(ICD10=="C50" | is.na(ICD10))



########################## processing of expression values  #####################################
exprs <- exprs(ex_run_1$lumi_batch)
negative_c <- ex_run_1$negative_controls

# background correction of expression values
# nowaclean package contains wrappers to bioconductor functionality 
# commonly used in nowac: https://github.com/3inar/nowaclean
exprs_corrected <- nowaclean::corrected(exprs, t(negative_c))
exprs_log <- log2(exprs_corrected)      

bat <- ex_run_1$lumi_batch
exprs(bat) <- exprs_log

normed <- nowaclean::normalized(bat)        # quantile normalization
filterd <- nowaclean::filtered(normed, pval = 0.01, fval=.15, verbose=T)   # filter out genes below detection threshold
aggd <- nowaclean::probe_aggregated(filterd, verbose=T)     # aggregate separate probes into one value
symbols <- nowaclean::gene_names(aggd)                      # translate probe ids to gene symbols

rm_these <- names(which(symbols[, "Symbol"] == ""))     # remove genes with missing symbol

# extract questionnaire data that corresponds to gene expression, randomize sample ids
two_yrs <- which(predictors$followup <= 2*365)

keepers <- c(predictors$labnr[two_yrs], predictors$match_labnr[two_yrs])

co_variables <- predictors %>% filter(labnr %in% sampleNames(aggd)) %>% 
  filter(match_labnr %in% sampleNames(aggd)) %>% filter(labnr %in% keepers) %>% dplyr::select(-lpnr)
cov_superset <- predictors %>% filter(case_control=="case") %>% 
   filter(!(labnr %in% keepers)) %>% dplyr::select(-lpnr) %>%
  dplyr::select(-labnr, -match_labnr)

labnos <- sample(co_variables$labnr)

newnames <- paste0("meta_", as.character(1:length(keepers)))
names(newnames) <- keepers

gene_exp <- exprs(aggd)[!is.element(rownames(exprs(aggd)), rm_these), labnos]
colnames(gene_exp) <- newnames[colnames(gene_exp)]
rownames(gene_exp) <- symbols[rownames(gene_exp), "Symbol"]

co_variables$labnr <- newnames[co_variables$labnr]
co_variables$match_labnr <- newnames[co_variables$match_labnr]

gene_exp <- t(gene_exp)

save(gene_exp, co_variables, cov_superset, file = "metastasis.rda")


