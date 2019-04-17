# Load some useful packages
library(tidyverse)  # Set of packages including ggplot2
library(caret)  # Used for knn and cross-validation
library(car)  # Calculates VIF to check for multicollinearity
library(ResourceSelection)  # Hosmer-Lemeshow goodness-of-fit test
library(survey)  # Build glm using weights

# Bring the chms_2018 dataset into R
study_data = read_csv("chms_2018.csv")

# Make SMK_12 and CLC_SEX a factor
study_data$SMK_12 = as.factor(study_data$SMK_12)
study_data$CLC_SEX = as.factor(study_data$CLC_SEX)

# Choose seeds so that results from imputing censored data is fixed.
# These numbers were drawn from a Uniform distribution with boundaries (1, 10^9). This will not be truly random, but will suffice for our purposes
seed1 = 336741157
seed2 = 874515716


#----------------------------------------------------------------------------------------------------------------
# Estimate censored data

# The values for LAB_BCD and LAB_BHG have been set to 999.5 when the patients recorded values are too low to be measured (ie. they are below the LOD). 
# First understand how common this is
BCD_LOD = !is.na(study_data$LAB_BCD) & study_data$LAB_BCD == 999.5 
sum(BCD_LOD)  # 999.5 appears 54 times in LAB_BCD
BHG_LOD = !is.na(study_data$LAB_BHG) & study_data$LAB_BHG == 999.5 
sum(BHG_LOD)  # 999.5 appears 599 times in LAB_BHG

# The values below the LOD are replaced with randomly drawn elements from the interval (0, LOD); this is called jittering
set.seed(seed1)
study_data$LAB_BCD[!is.na(study_data$LAB_BCD) & study_data$LAB_BCD == 999.5] = runif(sum(BCD_LOD), 0, 0.71)
set.seed(seed2)
study_data$LAB_BHG[!is.na(study_data$LAB_BHG) & study_data$LAB_BHG == 999.5] = runif(sum(BHG_LOD), 0, 2.1)

# Confirm that all values below the LOD have been replaced.
sum(study_data$LAB_BCD == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BCD
sum(study_data$LAB_BHG == 999.5, na.rm = TRUE)  # 999.5 appears 0 times in LAB_BHG


#----------------------------------------------------------------------------------------------------------------
# Use K-Nearest-Neighbors (KNN) to estimate missing values
# Check how many missing values there are in this data set
sum(is.na(study_data))  # 175 missing values

# Return variables that contain missing values, as well as their respective number of missing values
n_missing = colSums(is.na(study_data)) # number of missing values in each column of data set
matrix(c(names(study_data)[n_missing > 0], n_missing[n_missing > 0]), nrow = 2, byrow = TRUE)

# Are there observations that contain multiple missing values?
multiple_missing = rowSums(is.na(study_data))  # number of missing values in each row of data set
study_data[multiple_missing >= 2, seq(2, 8)] # return observations that have multiple missing values

# To simplify the following estimation, will estimate the LAB_BHG values for these two observations with the 
# sample mean of LAB_BHG. This insures that every observation has no more than 1 missing value
study_data$LAB_BHG[multiple_missing >= 2] = mean(study_data$LAB_BHG, na.rm = TRUE) 

multiple_missing = rowSums(is.na(study_data)) 
sum(multiple_missing > 1) # No observations have more than one missing value now

# Use cross-validation to find the best value of K to use in KNN 
ctrl <- trainControl(method="cv", number = 5) 

# Using HWMDBMI as the response variable, find the best value of K using as the training set only the observations 
# that have no missing values
knn_HWMDBMI <- train(HWMDBMI ~ as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE+LAB_BCD+LAB_BHG, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
# Estimate the missing values in HWMDBMI
study_data$HWMDBMI[is.na(study_data$HWMDBMI)] = predict(knn_HWMDBMI, study_data[is.na(study_data$HWMDBMI),])

# Repeat the above steps for the other variables with missing data
knn_LAB_BCD <- train(LAB_BCD ~ HWMDBMI+as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE+LAB_BHG, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BCD[is.na(study_data$LAB_BCD)] = predict(knn_LAB_BCD, study_data[is.na(study_data$LAB_BCD),])

knn_LAB_BHG <- train(LAB_BHG ~ LAB_BCD+HWMDBMI+as.factor(HIGHBP)+SMK_12+CLC_SEX+CLC_AGE, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$LAB_BHG[is.na(study_data$LAB_BHG)] = predict(knn_LAB_BHG, study_data[is.na(study_data$LAB_BHG),])

knn_SMK_12 <- train(SMK_12 ~ LAB_BHG+LAB_BCD+HWMDBMI+as.factor(HIGHBP)+CLC_SEX+CLC_AGE, data = study_data, subset = (multiple_missing == 0), method = "knn", preProcess = c("center","scale"), trControl = ctrl, tuneLength = 30)
study_data$SMK_12[is.na(study_data$SMK_12)] = predict(knn_SMK_12, study_data[is.na(study_data$SMK_12),])

sum(is.na(study_data))  # 0 missing values


#----------------------------------------------------------------------------------------------------------------
# Do some EDA

# Check correlations between continuous explanatory variables (CLC_AGE, HWMDBMI, LAB_BCD, and LAB_BHG)
cor(study_data[, 4:8])
pairs(study_data[2:8])

# Variables are not linearly correlated. Thus, we can maybe ignore other explanatory variables when plotting log odds 
# against a single variable

eda = study_data[, 1:9] %>%
  mutate(CLC_AGE_CAT = cut_number(CLC_AGE, 30), HWMDBMI_CAT = cut_number(HWMDBMI, 30), LAB_BCD_CAT = cut_number(LAB_BCD, 20), LAB_BHG_CAT = cut_number(LAB_BHG, 30))

# Group by categorical CLC_AGE and for each group, calculate log odds and sum of survey weights
eda1 = eda %>%
  group_by(CLC_AGE_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))), weights = sum(WGT_FULL))

# Plot logit(HIGHBP) against our categorical CLC_AGE, with size proportional to the relative weight
ggplot(eda1, aes(CLC_AGE_CAT, logit)) +
  geom_point(aes(size = weights)) +
  xlab("Age") + 
  ylab("Log Odds of Mean Hypertension") + 
  ggtitle("Relation Between Logit and Age is Linear") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1), axis.title.x = element_text(size = 15)) +
  scale_size_continuous("Survey Weights")

# No higher order terms for CLC_AGE are needed

# Repeat for HWMDBMI
eda2 = eda %>%
  group_by(HWMDBMI_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))), weights = sum(WGT_FULL))

# Plot logit(HIGHBP) against our categorical HWMDBMI, with size proportional to the relative weight
ggplot(eda2, aes(HWMDBMI_CAT, logit)) +
  geom_point(aes(size = weights)) +
  xlab("Body Mass Index") + 
  ylab("Log Odds of Mean Hypertension") + 
  ggtitle("Relation Between Logit and BMI Possibly Quadratic") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1), axis.title.x = element_text(size = 15)) +
  scale_size_continuous("Survey Weights")

# A quadratic term for HWMDBMI should be added

# Repeat for LAB_BCD
eda3 = eda %>%
  group_by(LAB_BCD_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))), weights = sum(WGT_FULL))

# Plot logit(HIGHBP) against our categorical HWMDBMI, with size proportional to the relative weight
ggplot(eda3, aes(LAB_BCD_CAT, logit)) +
  geom_point(aes(size = weights)) +
  xlab("Concentration of Blood Cadmium (nmol/L)") + 
  ylab("Log Odds of Mean Hypertension") + 
  ggtitle("Relation Between Logit and Blood Cadmium Possibly Quadratic") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1), axis.title.x = element_text(size = 15)) +
  scale_size_continuous("Survey Weights")

# A quadratic term for LAB_BCD should be added

# Repeat for LAB_BHG
eda4 = eda %>%
  group_by(LAB_BHG_CAT) %>%
  summarize(logit = log(mean(HIGHBP-1)/(1 - mean(HIGHBP-1))), weights = sum(WGT_FULL))

# Plot logit(HIGHBP) against our categorical HWMDBMI, with size proportional to the relative weight
ggplot(eda4, aes(LAB_BHG_CAT, logit)) +
  geom_point(aes(size = weights)) +
  xlab("Concentration of Blood Mercury (nmol/L)") + 
  ylab("Log Odds of Mean Hypertension") + 
  ggtitle("Relation Between Logit and Blood Mercury Possibly Quadratic") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.y = element_text(size = 15)) +
  theme(axis.text.x = element_text(angle = -30, hjust = 0, vjust = 1), axis.title.x = element_text(size = 15)) +
  scale_size_continuous("Survey Weights")

# A quadratic term for LAB_BHG could be added


#----------------------------------------------------------------------------------------------------------------
# Build three models each using a different link function

# Uses log odds link function
model_logit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'logit'), study_data)

# Uses inverse CDF of standard Normal distribution link function
model_probit = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'probit'), study_data)

# Uses complimentary log-log link function
model_cloglog = glm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, binomial(link = 'cloglog'), study_data)


#----------------------------------------------------------------------------------------------------------------
# Hosmer-Lemeshow statistic

# See that the minimum estimated probabilities that HIGHBP=2 are all greater than 10%
min(fitted.values(model_logit))
min(fitted.values(model_probit))
min(fitted.values(model_cloglog))

# The minimum estimated probabilities that HIGHBP=1 are also greater than 10%
max(fitted.values(model_logit))
max(fitted.values(model_probit))
max(fitted.values(model_cloglog))

# Calculate the Hosmer-Lemeshow statistic multiple times for each of the three models
pvalues_logit = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_logit$y, y = fitted.values(model_logit))["p.value",]
pvalues_probit = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_probit$y, y = fitted.values(model_probit))["p.value",]
pvalues_cloglog = sapply(X = seq(10, 60, 1), FUN = hoslem.test, x = model_cloglog$y, y = fitted.values(model_cloglog))["p.value",]

pvalues_logit = unlist(pvalues_logit)
pvalues_probit = unlist(pvalues_probit)
pvalues_cloglog = unlist(pvalues_cloglog)

# Plot distribution of Hosmer-Lemeshow statistics
ggplot(mapping = aes(x = factor(rep(c("Logit", "Probit", "Cloglog"),  each = 51), levels = c("Logit", "Probit", "Cloglog")),  y = c(pvalues_logit, pvalues_probit, pvalues_cloglog))) + 
  geom_boxplot() +
  xlab("Link Functions") + 
  scale_y_continuous("P-values", c(0.2, 0.4, 0.6, 0.8), limits = c(0, 1)) +
  ggtitle("Different Link Functions have Similiar P-values for Hosmer-Lemeshow Statistic") +
  theme(plot.title = element_text(hjust = 0.5))


#----------------------------------------------------------------------------------------------------------------
# Continue model adequacy checking

# Calculate the Generalized Variance Inflation Factors (GVIF) for each predictor
vif(model_logit) 

# Build residual plot
ggplot(mapping = aes(model_logit$linear.predictors, residuals(model_logit, "pearson"))) + 
  geom_point() +
  xlab("Predicted Log Odds Ratio") + ylab("Pearson Residuals") +
  ggtitle("Residual Plot with No Replicates is Uninformative") +
  theme(plot.title = element_text(hjust = 0.5))


#----------------------------------------------------------------------------------------------------------------
# Create a few functions that will be useful when doing hypothesis tests below

# For the given model, returns a dataframe with three columns (parameter estimates, std. errors, and p-values) 
param_estimates = function(model) {
  return(as_tibble(summary(model)$coefficients, rownames = "Variable") %>%
           select('Variable', 'Estimate', 'Std. Error') %>%
           mutate('p_value' = rep(NA, nrow(.))))
}

# For the given model, adds to the given dataframe three columns 
# (weighted parameter estimates, weighted std. errors, and weighted p-values) 
weighted_param_estimates = function(model, data_frame) {
  weighted_data_frame = as_tibble(summary(model)$coefficients, rownames = "Variable") %>%
    select('Variable', 'Estimate', 'Std. Error') %>%
    rename('Weighted Estimate' = 'Estimate', 'Weighted Std. Error' = 'Std. Error') %>%
    mutate('Weighted p_value' = rep(NA, nrow(.)))
  
  complete_data_frame = full_join(data_frame, weighted_data_frame, 'Variable')
  complete_data_frame = complete_data_frame[, c(1, 2, 5, 3, 6, 4, 7)]
}

# Use hypothesis tests on main effects and quadratic effects together
testing_main_and_quadratic = function(model, weighted_model, data_frame) {
  data_frame[nrow(data_frame)+1, 1] = "HWMDBMI + HWMDBMI^2"
  data_frame[nrow(data_frame), 2:8] = c(rep(NA, 4), anova(update(model, ~ .-HWMDBMI-I(HWMDBMI^2)), model, test = "LRT")$'Pr(>Chi)'[2], regTermTest(weighted_model, ~HWMDBMI+I(HWMDBMI^2), method = 'LRT')$p, NA)
  
  data_frame[nrow(data_frame)+1, 1] = "LAB_BCD + LAB_BCD^2"
  data_frame[nrow(data_frame), 2:8] = c(rep(NA, 4), anova(update(model, ~ .-LAB_BCD-I(LAB_BCD^2)), model, test = "LRT")$'Pr(>Chi)'[2], regTermTest(weighted_model, ~LAB_BCD+I(LAB_BCD^2), method = 'LRT')$p, NA)
  
  data_frame[nrow(data_frame)+1, 1] = "LAB_BHG + LAB_BHG^2"
  data_frame[nrow(data_frame), 2:8] = c(rep(NA, 4), anova(update(model, ~ .-LAB_BHG-I(LAB_BHG^2)), model, test = "LRT")$'Pr(>Chi)'[2], regTermTest(weighted_model, ~LAB_BHG+I(LAB_BHG^2), method = 'LRT')$p, NA)
  
  return(data_frame)
}


#----------------------------------------------------------------------------------------------------------------
# Use the package 'survey' to integrate weights into model
# Here, 'weights' represent the sample weights, and 'repweights' are the replicate weights 
glm_design = svrepdesign(study_data[, c(2:8)], repweights = select(study_data, starts_with("BS")), weights = study_data$WGT_FULL, type = 'bootstrap')
model_weighted = svyglm(as.factor(HIGHBP) ~ SMK_12+CLC_SEX+CLC_AGE+HWMDBMI+LAB_BCD+LAB_BHG, design = glm_design, family = quasibinomial(), data = study_data)

# Add some quadratic terms to logit model
model_logit_quad = update(model_logit, ~ .+I(HWMDBMI^2)+I(LAB_BCD^2)+I(LAB_BHG^2))

# Now add a few quadratic terms to weighted model
model_weighted_quad = update(model_weighted, ~ .+I(HWMDBMI^2)+I(LAB_BCD^2)+I(LAB_BHG^2))


#----------------------------------------------------------------------------------------------------------------
# Draw conclusions about individual parameters for unweighted and weighted model

# For logit model, create dataframe containing MLE estimates and s.e's for parameters
#logit_likelihood_tests = param_estimates(model_logit)

# Add likelihood ratio statistic p-values to the above dataframe
#logit_likelihood_tests[3:8, 'p_value'] = as_tibble(drop1(model_logit, test = "LRT"))[2:7, 'Pr(>Chi)']


# Since it seems likely that a weighted maximum likelihood method is being used to estimate parameters (see p. 268 of Fitting 
# Regression Models to Survey Data), we cannot just use the regular likelihood ratio test. Use regTermTest() from survey 
# package to perform hypothesis testing using the Rao-Scott working likelihood ratio test instead (see p.273 of same paper)

# For weighted logit model, create dataframe containing MLE estimates and s.e's for parameters
#logit_likelihood_tests = weighted_param_estimates(model_weighted, logit_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
#logit_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_weighted, 'SMK_12', method = 'LRT')$p
#for(i in c(4:8)) {
#  logit_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_weighted, logit_likelihood_tests$Variable[i], method = 'LRT')$p
#}

# Add Design Effects column
#logit_likelihood_tests = logit_likelihood_tests %>% 
#  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")
#logit_likelihood_tests


# Repeat with quadratic terms
# First unweighted logit model
logit_quad_likelihood_tests = param_estimates(model_logit_quad)
logit_quad_likelihood_tests[3:11, 'p_value'] = as_tibble(drop1(model_logit_quad, test = "LRT"), rownames = "Variable")[2:10, 'Pr(>Chi)']


# Now weighted model
logit_quad_likelihood_tests = weighted_param_estimates(model_weighted_quad, logit_quad_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
logit_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_weighted_quad, 'SMK_12', method = 'LRT')$p
logit_quad_likelihood_tests[4, 'Weighted p_value'] = regTermTest(model_weighted_quad, 'CLC_SEX', method = 'LRT')$p
for(i in c(5:11)) {
  logit_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_weighted_quad, logit_quad_likelihood_tests$Variable[i], method = 'LRT')$p
}

# Add Design Effects column
logit_quad_likelihood_tests = logit_quad_likelihood_tests %>% 
  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")

# Test main effects and higher-order terms simultaneously
logit_quad_likelihood_tests = testing_main_and_quadratic(model_logit_quad, model_weighted_quad, logit_quad_likelihood_tests)
logit_quad_likelihood_tests[c(2:14), c(1:3, 6:8)]


#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_SEX and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_SEX interaction terms
#model_gender = update(model_logit, ~ . +CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
#gender_likelihood_tests = param_estimates(model_gender)

# Add likelihood ratio statistic p-values to the above dataframe
#gender_likelihood_tests[c(3:8, 10:14), 'p_value'] = as_tibble(drop1(model_gender, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']


# Repeat for weighted model
# Build weighted model with CLC_SEX interaction terms
#model_gender_weighted = update(model_weighted, ~ .+CLC_SEX:SMK_12+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_SEX:LAB_BCD+CLC_SEX:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
#gender_likelihood_tests = weighted_param_estimates(model_gender_weighted, gender_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
#gender_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_gender_weighted, 'SMK_12', method = 'LRT')$p
#gender_likelihood_tests[10, 'Weighted p_value'] = regTermTest(model_gender_weighted, 'SMK_12:CLC_SEX', method = 'LRT')$p
#for(i in c(5:8, 11:14)) {
#  gender_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_gender_weighted, gender_likelihood_tests$Variable[i], method = 'LRT')$p
#}

# Add Design Effects column
#gender_likelihood_tests = gender_likelihood_tests %>% 
#  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")
#gender_likelihood_tests


# Repeat with quadratic terms
# First unweighted model with sex interaction terms
#model_gender_quad = update(model_logit_quad, ~ . +CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI)
#gender_quad_likelihood_tests = param_estimates(model_gender_quad)
#gender_quad_likelihood_tests[3:13, 'p_value'] = as_tibble(drop1(model_gender_quad, scope = ~ ., test = "LRT"), rownames = "Variable")[2:12, 'Pr(>Chi)']

# Now weighted model
#model_gender_weighted_quad = update(model_weighted_quad, ~ .+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI)
#gender_quad_likelihood_tests = weighted_param_estimates(model_gender_weighted_quad, gender_quad_likelihood_tests)

#gender_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_gender_weighted_quad, 'SMK_12', method = 'LRT')$p
#for(i in 4:13) {
#  gender_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_gender_weighted_quad, gender_quad_likelihood_tests$Variable[i], method = 'LRT')$p
#}

# Add Design Effects column
#gender_quad_likelihood_tests = gender_quad_likelihood_tests %>% 
#  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")

# Test main effects and higher-order terms simultaneously
#gender_quad_likelihood_tests = testing_main_and_quadratic(model_gender_quad, model_gender_weighted_quad, gender_quad_likelihood_tests)
#gender_quad_likelihood_tests[, c(1, 6:8)]


#----------------------------------------------------------------------------------------------------------------
# Now add Sex and Age interaction effects that involve the three significant risk factors: Age, Sex, and BMI

# First unweighted model 
model_interactions_quad = update(model_logit_quad, ~ . +CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_AGE:HWMDBMI)
interactions_quad_likelihood_tests = param_estimates(model_interactions_quad)
interactions_quad_likelihood_tests[3:14, 'p_value'] = as_tibble(drop1(model_interactions_quad, scope = ~ ., test = "LRT"), rownames = "Variable")[2:13, 'Pr(>Chi)']

# Now weighted model
model_interactions_weighted_quad = update(model_weighted_quad, ~ .+CLC_SEX:CLC_AGE+CLC_SEX:HWMDBMI+CLC_AGE:HWMDBMI)
interactions_quad_likelihood_tests = weighted_param_estimates(model_interactions_weighted_quad, interactions_quad_likelihood_tests)

interactions_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_interactions_weighted_quad, 'SMK_12', method = 'LRT')$p
interactions_quad_likelihood_tests[4, 'Weighted p_value'] = regTermTest(model_interactions_weighted_quad, 'CLC_SEX', method = 'LRT')$p
interactions_quad_likelihood_tests[12, 'Weighted p_value'] = regTermTest(model_interactions_weighted_quad, 'CLC_SEX:CLC_AGE', method = 'LRT')$p
interactions_quad_likelihood_tests[13, 'Weighted p_value'] = regTermTest(model_interactions_weighted_quad, 'CLC_SEX:HWMDBMI', method = 'LRT')$p
for(i in c(7:11, 14)) {
  interactions_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_interactions_weighted_quad, interactions_quad_likelihood_tests$Variable[i], method = 'LRT')$p
}

# Add Design Effects column
interactions_quad_likelihood_tests = interactions_quad_likelihood_tests %>% 
  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")

# Test main effects and higher-order terms simultaneously
interactions_quad_likelihood_tests = testing_main_and_quadratic(model_interactions_quad, model_interactions_weighted_quad, interactions_quad_likelihood_tests)
interactions_quad_likelihood_tests[, c(1:3, 6:8)]

#----------------------------------------------------------------------------------------------------------------
# Time for some log odds plots!

# Create a data frame that holds the explanatory variables and some of the estimated coefficients from both the weighted and unweighted model
log_odds = study_data[, 2:8] %>%
  mutate("SEX_coef" = model_interactions_quad$coefficients["CLC_SEX2"], "AGE_coef" = model_interactions_quad$coefficients["CLC_AGE"], "BMI_coef" = model_interactions_quad$coefficients["HWMDBMI"]) %>%
  mutate("SEX_AGE_coef" = model_interactions_quad$coefficients["CLC_SEX2:CLC_AGE"], "SEX_BMI_coef" = model_interactions_quad$coefficients["CLC_SEX2:HWMDBMI"], "AGE_BMI_coef" = model_interactions_quad$coefficients["CLC_AGE:HWMDBMI"]) %>%
  mutate("SEX_weight_coef" = model_interactions_weighted_quad$coefficients["CLC_SEX2"], "AGE_weight_coef" = model_interactions_weighted_quad$coefficients["CLC_AGE"], "BMI_weight_coef" = model_interactions_weighted_quad$coefficients["HWMDBMI"]) %>%
  mutate("SEX_AGE_weight_coef" = model_interactions_weighted_quad$coefficients["CLC_SEX2:CLC_AGE"], "SEX_BMI_weight_coef" = model_interactions_weighted_quad$coefficients["CLC_SEX2:HWMDBMI"], "AGE_BMI_weight_coef" = model_interactions_weighted_quad$coefficients["CLC_AGE:HWMDBMI"])
  

# For unweighted and weighted model, plot odds ratio for males and females, against Age at mean BMI
ggplot(log_odds, aes(CLC_AGE)) +
  geom_line(aes(y = exp(SEX_coef + CLC_AGE*SEX_AGE_coef + mean(HWMDBMI)*SEX_BMI_coef), color = 'Unweighted Model'), size = 1) +
  geom_line(aes(y = exp(SEX_weight_coef + CLC_AGE*SEX_AGE_weight_coef + mean(HWMDBMI)*SEX_BMI_weight_coef), color = 'Weighted Model'), size = 1) +
  xlab("Age") +
  ylab("Odds Ratio for Females to Males") + 
  ggtitle("Effect of Gender on HIGHBP Decreases with Age at Mean BMI") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  scale_colour_manual(name = "Model Type", values = c('Unweighted Model' = "red", 'Weighted Model' = "blue"))


Middle_Aged = 50 # Represents middle aged people (halfway between 40 and 60)
Older = 70 # Represents older people (halfway between 60 and 80)

# For unweighted and weighted model, plot odds ratio for middle aged and older people, against Sex at mean BMI
ggplot(log_odds, aes(CLC_SEX)) +
  geom_point(aes(y = exp((Middle_Aged-Older)*AGE_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_coef + (Middle_Aged-Older)*mean(HWMDBMI)*AGE_BMI_coef), color = 'Unweighted Model'), size = 3) +
  geom_point(aes(y = exp((Middle_Aged-Older)*AGE_weight_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_weight_coef + (Middle_Aged-Older)*mean(HWMDBMI)*AGE_BMI_weight_coef), color = 'Weighted Model'), size = 3) +
  xlab("Sex") +
  ylab("Odds Ratio for Middle Aged to Oldest Group") + 
  ggtitle("Effect of Age on HIGHBP Stronger for Females") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  scale_x_discrete(labels = c('Male', 'Female')) +
  scale_colour_manual(name = "Model Type", values = c('Unweighted Model' = "red", 'Weighted Model' = "blue"))

Overweight = 27.5  # Halfway between 25 and 30, which is typically considered to be the overweight group
Healthy = 22.5  # Halfway between 20 and 25, which is typically considered to be the healthy weight group 

# For unweighted and weighted model, plot odds ratio for Healthy and Overweight group, against Age for both male and female
ggplot(log_odds, aes(CLC_AGE)) +
  geom_line(aes(linetype = CLC_SEX, y = exp((Overweight-Healthy)*BMI_coef + (Overweight-Healthy)*(as.integer(study_data$CLC_SEX)-1)*SEX_BMI_coef + (Overweight-Healthy)*CLC_AGE*AGE_BMI_coef), color = 'Unweighted Model'), size = 1) +
  geom_line(aes(linetype = CLC_SEX, y = exp((Overweight-Healthy)*BMI_weight_coef + (Overweight-Healthy)*(as.integer(study_data$CLC_SEX)-1)*SEX_BMI_weight_coef + (Overweight-Healthy)*CLC_AGE*AGE_BMI_weight_coef), color = 'Weighted Model'), size = 1) +
  xlab("Age") +
  ylab("Odds Ratio for Overweight to Healthy Weight") + 
  ggtitle("Effect of Weight on HIGHBP Increases with Age at Mean BMI") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  scale_colour_manual(name = "Model Type", values=c('Unweighted Model' = "red", 'Weighted Model' = "blue")) +
  scale_linetype_discrete(name = "Sex", breaks = c(2, 1), labels = c('Female', 'Male'))

# For unweighted and weighted model, plot odds ratio for Healthy and Overweight group, against Sex at mean Age
ggplot(log_odds, aes(CLC_SEX)) +
  geom_point(aes(y = exp((Overweight-Healthy)*BMI_coef + (Overweight-Healthy)*(as.integer(study_data$CLC_SEX)-1)*SEX_BMI_coef + (Overweight-Healthy)*mean(CLC_AGE)*AGE_BMI_coef), color = 'Unweighted Model'), size = 3) +
  geom_point(aes(y = exp((Overweight-Healthy)*BMI_weight_coef + (Overweight-Healthy)*(as.integer(study_data$CLC_SEX)-1)*SEX_BMI_weight_coef + (Overweight-Healthy)*mean(CLC_AGE)*AGE_BMI_weight_coef), color = 'Weighted Model'), size = 3) +
  xlab("Sex") +
  ylab("Odds Ratio for Overweight to Healthy Weight") + 
  ggtitle("Effect of Weight on HIGHBP Weaker for Females") +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.x = element_text(size = 15), axis.title.y = element_text(size = 15)) +
  scale_x_discrete(labels = c('Male', 'Female')) +
  scale_colour_manual(name = "Model Type", values = c('Unweighted Model' = "red", 'Weighted Model' = "blue"))



  
  
  
  
  
  
#ggplot(log_odds, aes(HWMDBMI)) +
#  geom_line(aes(y = exp(SEX_coef + mean(CLC_AGE)*SEX_AGE_coef + HWMDBMI*SEX_BMI_coef), color = 'Unweighted Model'), size = 1) +
#  geom_line(aes(y = exp(SEX_weight_coef + mean(CLC_AGE)*SEX_AGE_weight_coef + HWMDBMI*SEX_BMI_weight_coef), color = 'Weighted Model'), size = 1) +
#  xlab("Body Mass Index") +
#  ylab("Odds Ratio for Females to Males") + 
#  ggtitle("Odds Ratio Decreases with BMI Faster in Weighted Model") +
#  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
#  scale_colour_manual(name = "Model Type", values = c('Unweighted Model' = "red", 'Weighted Model' = "blue"))
  
# For unweighted and weighted model, plot odds ratio for older and middle-aged people, against Sex at mean BMI
#ggplot(log_odds, aes(CLC_SEX)) +
#  geom_point(aes(y = exp((Middle_Aged-Older)*AGE_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_coef + (Middle_Aged-Older)*mean(HWMDBMI)*AGE_BMI_coef), color = 'Unweighted Model'), size = 3) +
#  geom_point(aes(y = exp((Middle_Aged-Older)*AGE_weight_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_weight_coef + (Middle_Aged-Older)*mean(HWMDBMI)*AGE_BMI_weight_coef), color = 'Weighted Model'), size = 3) +
#  xlab("Sex") +
#  ylab("Odds Ratio for Middle Aged Group to Oldest Group") + 
#  ggtitle("Odds Ratio Between Models More Different for Males") +
#  theme(plot.title = element_text(hjust = 0.5, size = 26)) +
#  scale_x_discrete(labels = c('Male', 'Female')) +
#  scale_colour_manual(name = "Model Type", values = c('Unweighted Model' = "red", 'Weighted Model' = "blue"))

# For unweighted and weighted model, plot odds ratio for older and middle-aged people, against BMI for males and females
#ggplot(log_odds, aes(HWMDBMI)) +
#  geom_line(aes(linetype = CLC_SEX, y = exp((Middle_Aged-Older)*AGE_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_coef + (Middle_Aged-Older)*HWMDBMI*AGE_BMI_coef), color = 'Unweighted Model'), size = 1) +
#  geom_line(aes(linetype = CLC_SEX, y = exp((Middle_Aged-Older)*AGE_weight_coef + (Middle_Aged-Older)*(as.integer(study_data$CLC_SEX)-1)*SEX_AGE_weight_coef + (Middle_Aged-Older)*HWMDBMI*AGE_BMI_weight_coef), color = 'Weighted Model'), size = 1) +
#  xlab("Body Mass Index") +
#  ylab("Odds Ratio for Middle Aged Group to Oldest Group") + 
#  ggtitle("Odds Ratio Decreases with BMI Faster in Weighted Model") +
#  theme(plot.title = element_text(hjust = 0.5, size = 30)) +
#  scale_colour_manual(name = "Model Type", values=c('Unweighted Model' = "red", 'Weighted Model' = "blue")) +
#  scale_linetype_discrete(name = "Sex", breaks = c(2, 1), labels = c('Female', 'Male'))
  

# For unweighted model, plot odds ratio for males and females, against age as bmi is fixed at 1st, 2nd, and 3rd quartiles
#ggplot(log_odds, aes(CLC_AGE)) +
#  geom_line(aes(y = exp(SEX_coef + CLC_AGE*SEX_AGE_coef + quantile(HWMDBMI)['25%']*SEX_BMI_coef)), color = 'green3', size = 1) +
#  geom_line(aes(y = exp(SEX_coef + CLC_AGE*SEX_AGE_coef + quantile(HWMDBMI)['50%']*SEX_BMI_coef)), color = 'blue', size = 1) +
#  geom_line(aes(y = exp(SEX_coef + CLC_AGE*SEX_AGE_coef + quantile(HWMDBMI)['75%']*SEX_BMI_coef)), color = 'red', size = 1) +
#  xlab("Age") +
#  ylab("Odds Ratio for Males to Females") + 
#  ggtitle("The Greater Risk of Hypertension for Males Decreases with Age") +
#  theme(plot.title = element_text(hjust = 0.5))
  
# Plot odds ratio for males and females, against bmi as age is fixed
#ggplot(log_odds, aes(HWMDBMI)) +
#  geom_line(aes(y = exp(SEX_coef + quantile(CLC_AGE)['25%']*SEX_AGE_coef + HWMDBMI*SEX_BMI_coef)), color = 'green3', size = 1) +
#  geom_line(aes(y = exp(SEX_coef + quantile(CLC_AGE)['50%']*SEX_AGE_coef + HWMDBMI*SEX_BMI_coef)), color = 'blue', size = 1) +
#  geom_line(aes(y = exp(SEX_coef + quantile(CLC_AGE)['75%']*SEX_AGE_coef + HWMDBMI*SEX_BMI_coef)), color = 'red', size = 1) +
#  xlab("BMI") +
#  ylab("Odds Ratio for Males to Females") + 
#  ggtitle("The Greater Risk of Hypertension for Males Decreases with BMI") +
#  theme(plot.title = element_text(hjust = 0.5))



#----------------------------------------------------------------------------------------------------------------
# Test whether interaction terms between CLC_AGE and other variables are non-zero, for both unweighted and weighted models

# Build logit model with CLC_AGE interaction terms
#model_age = update(model_logit, ~ . +CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
#age_likelihood_tests = param_estimates(model_age)

# Add likelihood ratio statistic p-values to the above dataframe
#age_likelihood_tests[c(3:8, 10:14), 'p_value'] = as_tibble(drop1(model_age, scope = ~ ., test = "LRT"))[2:12, 'Pr(>Chi)']

# Repeat for weighted model
# Build weighted model with CLC_AGE interaction terms
#model_age_weighted = update(model_weighted, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Create dataframe containing MLE estimates and s.e's for parameters
#age_likelihood_tests = weighted_param_estimates(model_age_weighted, age_likelihood_tests)

# Add Rao-Scott working likelihood ratio p-values to the above dataframe
#age_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_age_weighted, 'SMK_12', method = 'LRT')$p
#age_likelihood_tests[10, 'Weighted p_value'] = regTermTest(model_age_weighted, 'SMK_12:CLC_AGE', method = 'LRT')$p
#for(i in c(4, 6:8, 11:14)) {
#  age_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_age_weighted, age_likelihood_tests$Variable[i], method = 'LRT')$p
#}

# Add Design Effects column
#age_likelihood_tests = age_likelihood_tests %>% 
#  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")
#age_likelihood_tests


# Repeat with quadratic terms
# First unweighted model with age interaction terms
#model_age_quad = update(model_logit_quad, ~ . +CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
#age_quad_likelihood_tests = param_estimates(model_age_quad)
#age_quad_likelihood_tests[c(3:11, 13:17), 'p_value'] = as_tibble(drop1(model_age_quad, scope = ~ ., test = "LRT"), rownames = "Variable")[2:15, 'Pr(>Chi)']

# Now weighted model
#model_age_weighted_quad = update(model_weighted_quad, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
#age_quad_likelihood_tests = weighted_param_estimates(model_age_weighted_quad, age_quad_likelihood_tests)

#age_quad_likelihood_tests[3, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, 'SMK_12', method = 'LRT')$p
#age_quad_likelihood_tests[13, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, 'SMK_12:CLC_AGE', method = 'LRT')$p
#for(i in c(4, 6:11, 14:17)) {
#  age_quad_likelihood_tests[i, 'Weighted p_value'] = regTermTest(model_age_weighted_quad, age_quad_likelihood_tests$Variable[i], method = 'LRT')$p
#}

# Add Design Effects column
#age_quad_likelihood_tests = age_quad_likelihood_tests %>% 
#  mutate("Design Effect" = .$"Weighted Std. Error"/.$"Std. Error")

# Test main effects and higher-order terms simultaneously
#age_quad_likelihood_tests = testing_main_and_quadratic(model_age_quad, model_age_weighted_quad, age_quad_likelihood_tests)
#age_quad_likelihood_tests[, c(1:3, 6:8)]


#----------------------------------------------------------------------------------------------------------------
# Create a few more interesting models. All of these have all sex and age two way interactions
#model_agesex_quad = update(model_gender_quad, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
#model_agesex_weighted = update(model_gender_weighted, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)
#model_agesex_weighted_quad = update(model_gender_weighted_quad, ~ .+CLC_AGE:SMK_12+CLC_AGE:CLC_SEX+CLC_AGE:HWMDBMI+CLC_AGE:LAB_BCD+CLC_AGE:LAB_BHG)

# Check AIC and BIC for unweighted models
#AIC(model_logit, model_logit_quad, model_gender, model_gender_quad, model_age, model_age_quad, model_interaction_quad, model_all_interactions)
#BIC(model_logit, model_logit_quad, model_gender, model_gender_quad, model_age, model_age_quad, model_interaction_quad, model_all_interactions)

# Check AIC and BIC for weighted models
#AIC(model_weighted, model_weighted_quad, model_gender_weighted, model_gender_weighted_quad, model_age_weighted, model_age_weighted_quad, model_agesex_weighted, model_agesex_weighted_quad)
#BIC(model_weighted, model_weighted_quad, model_gender_weighted, model_gender_weighted_quad, model_age_weighted, model_age_weighted_quad, model_agesex_weighted, model_agesex_weighted_quad, maximal = model_agesex_weighted_quad)


#----------------------------------------------------------------------------------------------------------------
# Some Questions to Think About

# 1) In svrepdesign function, are the replication weights (our bootstrap weights) really of type 'bootstrap'? (see help page)
#    Need to better understand what Canadian Health Measures Survey (CHMS) means by 'bootstrap weights', and what author of survey package means by 'boostrap' replication weights

#    First talk about how variance is calculated (p. 240 and in section 2.3 of text)
#    From 'The Rao-Wu Rescaling Bootstrap: From theory to practice' (on github), seems likely that StatsCan uses the Rao-Wu rescaled bootstrap for the CHMS
#    as.svyrepdesign() seems to confirm that type = 'bootstrap' in svyrepdesign() can refer to Rao and Wu's rescaled boostrap
#    Section 2.3 mentions what scale in svrepdesign() should be for bootstrap. From experiments below, scale seems to equal 500-1 for bootstrap method
#    This agrees with 'Rescaled Bootstrap for stratified multistage sampling'

# 2) In svrepdesign function, not specifying repweights produces the warning: 'You must provide replication weights' 
#    Why is this?

#    I'm guessing that since replication weights are needed to estimate the variance of estimators, instead of just returning
#    point estimates without any kind of measure of error (which would be useless), the function just creates a warning

# 3) In svyglm function, letting family = binomial() produces the warning: 'In eval(family$initialize) : non-integer #successes in a binomial glm!'
#    Both the svyglm help file and the survey author's textbook (p.110) recommend using family = quasibinomial() to avoid these warnings, but no clear rationale is given

# 4) According to survey author's textbook (p.98), svyglm() does not seem to use maximum likelihood to estimate parameters,
#    and author concludes that we can't use likelihood ratio tests to compare nested models, but should use Wald statistic instead
#    Not sure if this is just for linear regression, or for logistic regression as well.

#    Actually, I think svyglm() uses pseudo-likelihood estimation (see p. 268 of Fitting Regression Models to Survey Data)
