##KATIE HOEMANN, Katholieke Universiteit Leuven; last updated 8 July 2021

# set directories
setwd("C:/Users/Katie/Documents/R/Category_Learning")
base_dir <- 'C:/Users/Katie/Documents/R/Category_Learning/'

# load required libraries (install if necessary)
library(readr)
library(brms)

# read in data
d <- read.csv("Category_Learning_IDA_text.csv")

# set up simple effects coding, following: https://stats.idre.ucla.edu/r/library/r-library-contrast-coding-systems-for-categorical-variables/#SIMPLE 
d$StudyS <- d$Study
d$StudyS <- as.factor(d$StudyS)
c <- contr.treatment(5)
coding.matrix <- matrix(rep(1/5,20),ncol=4)
simple.coding <- c-coding.matrix
contrasts(d$StudyS) <- simple.coding
d$ConditionS <- d$Condition
d$ConditionS <- as.factor(d$Condition)
contrasts(d$ConditionS) <- c(-0.5,+0.5)

# run Bayesian logistic mixed-effects model model - this is expensive so only do it once, save it, and then reload as necessary (comment lines out as appropriate)
# mlmB_IDA_full_S <- brm(Accuracy~StudyS*ConditionS*TrialCountC+(TrialCountC|PPID),
#                         data=d, family=bernoulli, prior=set_prior('normal(0, 3)'), iter=5000, chains=4, cores=4, save_pars=save_pars(all=TRUE), sample_prior=TRUE)
# save(mlmB_IDA_full_S,file=paste0(base_dir,'IDA_brms_full_S.rda'))
load(paste0(base_dir,'IDA_brms_full_S.rda'))

# get summary of fixed effects
mlmB_IDA_fixef <- fixef(mlmB_IDA_full_S)

# get odds ratios in addition to log odds for estimated betas, following: https://stats.idre.ucla.edu/r/dae/logit-regression/
mlmB_IDA_ORs <- exp(mlmB_IDA_fixef)

# test hypotheses, estimate Bayes factors, following: https://rstudio-pubs-static.s3.amazonaws.com/358672_09291d0b37ce43f08cf001cfd25c16c2.html#brmshypothesis
test_H1 <- hypothesis(mlmB_IDA_full_S,"TrialCountC>0")
H1 <- test_H1$hypothesis$Evid.Ratio # if Inf then use the number of samples (here 10000)

test_H2 <- hypothesis(mlmB_IDA_full_S,"ConditionS1=0")
H2 <- test_H2$hypothesis$Evid.Ratio 

## exploratory analysis A: studies 1 and 2 including factor for emotion category type (fear vs. novel)
d_s1s2 <- d[which(d$Study<3),]
d_s1s2$TypeS <- d_s1s2$Type
d_s1s2$TypeS <- as.factor(d_s1s2$Type)
contrasts(d_s1s2$TypeS) <- c(-0.5,+0.5)
d_s1s2$StudySa <- d_s1s2$Study
d_s1s2$StudySa <- as.factor(d_s1s2$StudySa)
contrasts(d_s1s2$StudySa) <- c(-0.5,+0.5)
# mlmB_IDA_s1s2_type_S <- brm(Accuracy~StudySa*ConditionS*TypeS*TrialCountC+(TrialCountC|PPID),
#                        data=d_s1s2, family=bernoulli, prior=set_prior('normal(0, 3)'), iter=5000, chains=4, cores=4, save_pars=save_pars(all=TRUE), sample_prior=TRUE)
# save(mlmB_IDA_s1s2_type_S,file=paste0(base_dir,'IDA_brms_s1s2_type_S.rda'))
load(paste0(base_dir,'IDA_brms_s1s2_type_S.rda'))
mlmB_IDA_s1s2_type_fixef <- fixef(mlmB_IDA_s1s2_type_S)
mlmB_IDA_s1s2_type_ORs <- exp(mlmB_IDA_s1s2_type_fixef)

## exploratory analysis B: studies 1-3 including factor for cue modality
d_s1s2s3 <- d[which(d$Study<4),]
d_s1s2s3$ModalityS <- d_s1s2s3$Modality
d_s1s2s3$ModalityS <- as.factor(d_s1s2s3$Modality)
c <- contr.treatment(3)
coding.matrix <- matrix(rep(1/3,6),ncol=2)
simple.coding <- c-coding.matrix
contrasts(d_s1s2s3$ModalityS) <- simple.coding
d_s1s2s3$StudySb <- d_s1s2s3$Study
d_s1s2s3$StudySb <- as.factor(d_s1s2s3$StudySb)
contrasts(d_s1s2s3$StudySb) <- simple.coding
# mlmB_IDA_s1s2s3_modality_S <- brm(Accuracy~StudySb*ConditionS*ModalityS*TrialCountC+(TrialCountC|PPID),
#                             data=d_s1s2s3, family=bernoulli, prior=set_prior('normal(0, 3)'), iter=5000, chains=4, cores=4, save_pars=save_pars(all=TRUE), sample_prior=TRUE)
# save(mlmB_IDA_s1s2s3_modality_S,file=paste0(base_dir,'IDA_brms_s1s2s3_modality_S.rda'))
load(paste0(base_dir,'IDA_brms_s1s2s3_modality_S.rda'))
mlmB_IDA_s1s2s3_modality_fixef <- fixef(mlmB_IDA_s1s2s3_modality_S)
mlmB_IDA_s1s2s3_modality_ORs <- exp(mlmB_IDA_s1s2s3_modality_fixef)