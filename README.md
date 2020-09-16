# Metastatic breast cancer and pre-diagnostic blood gene expression profiles – The Norwegian Women and Cancer (NOWAC) Post-genome Cohort
This repo contains the data processing and analysis scripts for the manuscript Holsbø and Olsen, 2020. It comprises the following files:

1. `1_gene_expression_prep.R` holds the initial processing of gene expression data (normalization, transformation, etc.,) and the extraction of questionnaire data.
1. `2_munge.R` holds further processing of these data, mostly filtering out of observations that do not fit the analysis.
1. `3a_fit_model.R` does the setting up and execution of the model specified in `3b_model_specification.stan`
1. `3b_model_specification.stan` defines a genewise logistic regression model, see manuscript for details.
1. `4_summaries_and_figures.R` holds the exploration, summaries, and plotting of the above model.

Some path names have been censored. This repo is intended for documentation only.
