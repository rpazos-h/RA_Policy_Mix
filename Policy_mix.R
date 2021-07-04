##--------------------------------------------------------------------##
##        Main: Fiscal and Monetary Policy Mix - OECD Countries       ## 
##                                                                    ##
## This script prepares the data set and estimates with a logit model ## 
##    the probability of the policy stance being pro-cyclical         ## 
##    or counter-cyclical, conditional on a vector of variables       ##  
##--------------------------------------------------------------------##

# Working directory

setwd('')


# Loading packages

libraries <- c('pglm','plm','tidyverse','readxl','lmtest',
               'StatMeasures', 'ExPanDaR', 'sampleSelection', 
               'sandwich', 'kableExtra', 'sjlabelled', 
               'texreg','reshape','httr')

new.packages <- libraries[!(libraries %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(libraries, require, character.only=TRUE)

#- A) Preparing dataset -----

# A.1) Load data

github_link <- "https://github.com/rpazos-h/RA_Policy_Mix/blob/main/Coordination_fiscal_monetary_database.xlsx?raw=true"

temp_file <- tempfile(fileext = ".xlsx")
req <- GET(github_link, 
           authenticate(Sys.getenv("GITHUB_PAT"), ""),
           write_disk(path = temp_file))
data <- readxl::read_excel(temp_file, sheet = "Dataset", skip = 3)

unlink(temp_file)


# A.2) Defining good/normal/bad times

growth_dummy <- data %>% 
  select(gdp_growth, country) %>%
  group_by(country) %>% 
  mutate(decile_gdp_growth = decile(vector = gdp_growth, decreasing=F)) %>% 
  ungroup

data <- cbind(data, "decile_gdp_growth" = growth_dummy[,3])
growth_dummy <- growth_dummy[complete.cases(growth_dummy),]

data <- data %>%
  mutate(good_times = ifelse(decile_gdp_growth>=7 & decile_gdp_growth<=8,1,0),
         vgood_times = ifelse(decile_gdp_growth>=9,1,0),
         bad_times = ifelse(decile_gdp_growth<=4 & decile_gdp_growth>=3,1,0),
         vbad_times = ifelse(decile_gdp_growth<=2,1,0),
         normal_times = ifelse(decile_gdp_growth>=4 & decile_gdp_growth<=5,1,0))


# A.3) Declaring dataset as panel data

data_ts <- pdata.frame(data, index=c("index","year"), drop.index = FALSE, row.names = TRUE)


# A.4) Creating lag and first differences variables

data_ts$output_gap_d <- diff(data_ts$output_gap)
data_ts$primary_balance_d <- diff(data_ts$primary_balance)
data_ts$int_rate_d <- diff(data_ts$int_rate)
data_ts$yields_10yr_lag <- lag(data_ts$yields_10yr)

data_ts$vix_lag <- lag(data_ts$vix)
data_ts$bank_crisis_lag <- lag(data_ts$bank_crisis)
data_ts$curr_crisis_lag <- lag(data_ts$curr_crisis)
data_ts$systemic_crisis_lag <- lag(data_ts$systemic_crisis)
data_ts$debt_crisis_lag <- lag(data_ts$debt_crisis)

data_ts$vgood_times_lag <- lag(data_ts$vgood_times)
data_ts$good_times_lag <- lag(data_ts$good_times)
data_ts$vbad_times_lag <- lag(data_ts$vbad_times)
data_ts$bad_times_lag <- lag(data_ts$bad_times)


# A.5) Adding new variables to data

data_1 <- data_ts %>% select(country, year, output_gap_d:debt_crisis_lag)
data <- merge(data, data_1, by=c("country","year"))


#- B) Table 2.1. Logit model for pro-cyclical fiscal policy----

# B.1) Models pc_fiscal

logit_1 <- glm(pc_fiscal ~  output_gap_d + debt  + yields_10yr +  er_reg,  data,
                family = binomial(link ='logit'))

logit_2 <- update(logit_1, . ~ . + inflation + mro_zerolb)

logit_3 <- update(logit_1, . ~ . + vix)

logit_4 <- update(logit_3, . ~ . + vix_lag)

logit_5 <- update(logit_4, . ~ . + fin_assist)

logit_6 <- update(logit_5, . ~ . + pc_monetary)

logit_7 <- update(logit_6, . ~ . + cc_monetary)


# B.2) Output tables

extract.glm <- function (model, include.nobs = TRUE, include.loglik = TRUE, ...) {
  s <- summary(model, ...)
  coefficient.names <- rownames(s$coefficients)
  coefficients <- s$coefficients[, 1]
  standard.errors <- s$coefficients[, 2]
  significance <- s$coefficients[, 4]
  loglik.value <- round(logLik(model),2)
  pseudoR <- round(1 - model$deviance / model$null.deviance,3)
  aic <- s$AIC
  bic <- s$BIC
  n <- nobs(model)
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.loglik == TRUE) {
    gof <- c(gof, loglik.value, pseudoR)
    gof.names <- c(gof.names, "Log-Likelihood", "Pseudo R sq")
    gof.decimal <- c(gof.decimal, TRUE, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Observations")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  tr <- createTexreg(coef.names = coefficient.names, coef = coefficients, 
                     se = standard.errors, pvalues = significance, gof.names = gof.names, 
                     gof = gof, gof.decimal = gof.decimal)
  return(tr)
}

setMethod("extract", signature = className("glm"), 
          definition = extract.glm)

list_models <- vector("list", 7)

for(i in 1:7) {
  Ps <- paste0("logit_",i)
  extract(get(Ps))
  list_models[[i]] <- get(Ps)
}

names(list_models) <- c('(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)')

labels = c('Intercept', 'og_diff', 'debt', 'yields_10yr', 'er_reg', 
           'Vix', 'Vix_lag', 'Fin assist', 'Inflation', 
           'ZLB', 'Pc_monetary', 'Cc_monetary')

screenreg(list_models,digits = 3, stars = c(0.01,0.05, 0.1), custom.coef.names = labels,
          reorder.coef = c(2:12,1))

htmlreg(list_models, file= "Output/Logit model/PC fiscal.doc", custom.coef.names = labels,reorder.coef = c(2:12,1),
        digits = 3, stars = c(0.01,0.05, 0.1),caption = "Table 2.1. Logit model for pro-cyclical fiscal policy", 
        caption.above = TRUE, margin = 4, 
        custom.header = list("Dependent Variable: Probability of pro-cyclical fiscal policy" = 1:7),
        inner.rules = 2, outer.rules = 3, padding=4, inline.css = FALSE, center = TRUE, doctype = TRUE, 
        html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)

# C) Table 2.2. Logit model for pro-cyclical monetary policy------

# C.1) Models pc_monetary

logit_1 <- glm(pc_monetary ~  output_gap_d + debt  + yields_10yr +  er_reg,  data,
               family = binomial(link ='logit'))

logit_2 <- update(logit_1, . ~ . + inflation + mro_zerolb)

logit_3 <- update(logit_1, . ~ . + vix + mro_zerolb)

logit_4 <- update(logit_3, . ~ . + vix_lag)

logit_5 <- update(logit_4, . ~ . + fin_assist)

logit_6 <- update(logit_3, . ~ . + pc_fiscal)

logit_7 <- update(logit_6, . ~ . + cc_fiscal)


# C.2) Output tables

list_models <- vector("list", 7)

for(i in 1:7) {
  Ps <- paste0("logit_",i)
  extract(get(Ps))
  list_models[[i]] <- get(Ps)
}

labels = c('Intercept', 'og_diff', 'debt', 'yields_10yr', 'er_reg', 
           'Inflation', 'ZLB', 'Vix', 'Vix_lag', 'Fin assist',
           'Pc_fiscal', 'Cc_fiscal')

names(list_models) <- c('(1)', '(2)', '(3)', '(4)', '(5)', '(6)', '(7)')

screenreg(list_models,digits = 3, stars = c(0.01,0.05, 0.1), custom.coef.names = labels,
          reorder.coef = c(2:5,8:10,6:7,11:12,1))

htmlreg(file= "Output/Logit model/PC monetary.doc",list_models, custom.coef.names = labels,reorder.coef = c(2:5,8:10,6:7,11:12,1),
        digits = 3, stars = c(0.01,0.05, 0.1), caption = "Table 2.2. Logit model for pro-cyclical monetary policy", 
        caption.above = TRUE, margin = 4, 
        custom.header = list("Dependent Variable: Probability of pro-cyclical monetary policy" = 1:7),
        inner.rules = 2, outer.rules = 3, padding=4, inline.css = FALSE, center = TRUE, doctype = TRUE, 
        html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)

# D) Table X.X. Logit model for counter-cyclical fiscal policy-----

# D.1) Models cc_fiscal

logit_1 <- glm(cc_fiscal ~  output_gap_d + debt  + yields_10yr +  er_reg + vix + vix_lag + fin_assist 
               + cc_monetary + pc_monetary,  data, family = binomial(link ='logit'))

# D.2) Output tables

labels = c('Intercept', 'og_diff', 'debt', 'yields_10yr', 'er_reg', 
           'Vix', 'Vix_lag', 'Fin assist', 'Pc_fiscal', 'Cc_fiscal')

extract(logit_1)
screenreg(logit_1, digits = 3, stars = c(0.01,0.05, 0.1), custom.coef.names = labels,
          reorder.coef = c(2:10,1))

htmlreg(file= "Output/Logit model/CC fiscal.doc",list('(1)' = logit_1), custom.coef.names = labels,reorder.coef = c(2:10,1),
        digits = 3, stars = c(0.01,0.05, 0.1), caption = "Table X.X Logit model for counter-cyclical fiscal policy", 
        caption.above = TRUE, margin = 4, 
        custom.header = list("Dependent Variable: cc_fiscal" = 1),
        inner.rules = 2, outer.rules = 3, padding=6, inline.css = FALSE, center = TRUE, doctype = TRUE, 
        html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)

# E) Table X.X. Logit model for counter-cyclical monetary policy-----

# E.1) Models cc_monetary

logit_1 <- glm(cc_monetary ~  output_gap_d + debt  + yields_10yr +  er_reg + vix + mro_zerolb + cc_fiscal 
               + pc_fiscal,  data, family = binomial(link ='logit'))

# E.2) Output tables

labels = c('Intercept', 'og_diff', 'debt', 'yields_10yr', 'er_reg', 
           'Vix', 'ZLB', 'Cc_fiscal', 'Pc_fiscal')

extract(logit_1)
screenreg(logit_1, digits = 3, stars = c(0.01,0.05, 0.1), custom.coef.names = labels,
          reorder.coef = c(2:9,1))

htmlreg(file= "Output/Logit model/CC monetary.doc",list('(1)' = logit_1), custom.coef.names = labels,reorder.coef = c(2:9,1),
        digits = 3, stars = c(0.01,0.05, 0.1), caption = "Table X.X Logit model for counter-cyclical monetary policy", 
        caption.above = TRUE, margin = 4, 
        custom.header = list("Dependent Variable: cc_monetary" = 1),
        inner.rules = 2, outer.rules = 3, padding=6, inline.css = FALSE, center = TRUE, doctype = TRUE, 
        html.tag = TRUE, head.tag = TRUE, body.tag = TRUE)
