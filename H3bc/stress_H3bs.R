# hypothesis testing

library(lmerTest)
library(brms)
library(EMAtools)
library(sjstats)

load('Stress_aligned.RData')

# interactions
data.filtered$secondary <- scale(data.filtered$secondary)
data.filtered$pss <- scale(data.filtered$pss)
data.filtered$resilience <- scale(data.filtered$resilience)

# pss ~ secondary
options(width = 2000)

# prior -> normal distribution 
prior.coef <- brms::prior(cauchy(0.,1),class='b')

# model 0
pss.0 <- brms::brm(pss ~ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415)

# model 1
pss.1 <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model 2
pss.2 <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# bayes factors
pss.10 <- bayes_factor(pss.1, pss.0)
pss.20 <- bayes_factor(pss.2, pss.0)
pss.21 <- bayes_factor(pss.2,pss.1)

log (pss.21$bf[1])

# hypothesis testing
hypothesis(pss.2,'secondary>0')

# interaction?
data.filtered$gender<-as.factor(data.filtered$gender)
data.filtered$education<-as.factor(data.filtered$education)
data.filtered$occupation<-as.factor(data.filtered$occupation)

pss.int <- brms::brm(pss ~ secondary*gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1+secondary|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef, iter=5000)

bf.int <- bayes_factor(pss.int,pss.2)
# without interaction -> better (BF 0)

# only with gender
pss.int.gender <- brms::brm(pss ~ gender*secondary + work_location + age+
                              SSS_faml+ relationship_status+ education+
                              (1+secondary|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.gender <- bayes_factor(pss.int.gender,pss.2)
# without gender interaction -> better BF .00236)

# only with education
pss.int.edu <- brms::brm(pss ~ gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ 
                       (1+secondary|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.edu <- bayes_factor(pss.int.edu,pss.2)
# BF = 0

# only with occupation
pss.int.ocu <- brms::brm(pss ~gender +education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1+secondary|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef, iter= 5000)

bf.ocu <- bayes_factor(pss.int.ocu,pss.2)
# BF = 0

# gender / education
pss.int.gender.edu <- brms::brm(pss ~ gender*secondary + work_location + age+
                              SSS_faml+ relationship_status+ secondary*education+
                                occupation+
                              (1+secondary|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.gender.edu <- bayes_factor(pss.int.gender.edu,pss.2)
# BF = 0

# gender / occupation
pss.int.gender.occu <- brms::brm(pss ~ gender*secondary + work_location + age+
                                  SSS_faml+ relationship_status+ education+
                                   occupation*secondary+
                                  (1+secondary|residing_country),
                                data=data.filtered, family = gaussian(),
                                cores=4,chains=4, save_pars = save_pars(all = T),
                                sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.gender.occu <- bayes_factor(pss.int.gender.occu,pss.2)
# 0

# education / occupation
pss.int.edu.occu <- brms::brm(pss ~ gender + work_location + age+
                                   SSS_faml+ relationship_status+ education*secondary+
                                   occupation*secondary+
                                   (1+secondary|residing_country),
                                 data=data.filtered, family = gaussian(),
                                 cores=4,chains=4, save_pars = save_pars(all = T),
                                 sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.edu.occu <- bayes_factor(pss.int.edu.occu,pss.2)
# BF 0 

# double check with freq
freq.pss.2 <- lmer(pss ~ secondary+ work_location + age+
                       SSS_faml+ relationship_status+ education+
                       (1+secondary|residing_country),
                     data=data.filtered)
conft.pss.2<-confint(freq.pss.2)

# mumin
library(MuMIn)
MuMIn::r.squaredGLMM(freq.pss.2)

freq.pss.1 <- lmer(pss ~ secondary+ work_location + age+
                     SSS_faml+ relationship_status+ education+
                     (1|residing_country),
                   data=data.filtered)
EMAtools::lme.dscore(freq.pss.1,data.filtered,'lme4')

# ICC
performance::icc(freq.pss.2)

freq.pss.int <- lmer(pss ~ secondary*gender + secondary*education + work_location + age+
                     SSS_faml+ relationship_status+ secondary*occupation+
                     (1+secondary|residing_country),
                   data=data.filtered)
freq.pss.gender <- lmer(pss ~ secondary*gender + education + work_location + age+
                       SSS_faml+ relationship_status+ 
                       (1+secondary|residing_country),
                     data=data.filtered)

save.image(file='Vaccine_H3bs.RData')



##### resilience

# model 0
res.0 <- brms::brm(resilience ~ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415)

# model 1
res.1 <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model 2
res.2 <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# bayes factors
res.10 <- bayes_factor(res.1, res.0)
res.20 <- bayes_factor(res.2, res.0)
res.21 <- bayes_factor(res.2,res.1)

# 1 is better, so random intercept model will be used
hypothesis(res.1, 'secondary < 0')

# int
res.int <- brms::brm(resilience ~ secondary*gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)
res.int.1 <- bayes_factor(res.int,res.1)
#BF 0 -> no interaction!

# gender
res.int.gender <- brms::brm(resilience ~ secondary*gender + education + work_location + age+
                       SSS_faml+ relationship_status+ occupation+
                       (1|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)
res.gender.1 <- bayes_factor(res.int.gender,res.1)
# .00

# education
res.int.edu<- brms::brm(resilience ~ gender + secondary*education + work_location + age+
                              SSS_faml+ relationship_status+ occupation+
                              (1|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)
res.edu.1 <- bayes_factor(res.int.edu,res.1)
# 0

# occupation
res.int.occu<- brms::brm(resilience ~ gender + education + work_location + age+
                          SSS_faml+ relationship_status+ occupation*secondary+
                          (1|residing_country),
                        data=data.filtered, family = gaussian(),
                        cores=4,chains=4, save_pars = save_pars(all = T),
                        sample_prior ='yes', seed=1660415,prior=prior.coef)
res.occu.1 <- bayes_factor(res.int.occu,res.1)
# 0

# gender / edu
res.int.gender.edu <- brms::brm(resilience ~ secondary*gender + secondary*education + work_location + age+
                              SSS_faml+ relationship_status+ occupation+
                              (1|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)
res.gender.edu.1 <- bayes_factor(res.int.gender.edu,res.1)
# 0

# gender / occu
res.int.gender.occu <- brms::brm(resilience ~ secondary*gender + education + work_location + age+
                                  SSS_faml+ relationship_status+ secondary*occupation+
                                  (1|residing_country),
                                data=data.filtered, family = gaussian(),
                                cores=4,chains=4, save_pars = save_pars(all = T),
                                sample_prior ='yes', seed=1660415,prior=prior.coef)
res.gender.occu.1 <- bayes_factor(res.int.gender.occu,res.1)
# 0

# educ / occu
res.int.edu.occu <- brms::brm(resilience ~ gender + secondary*education + work_location + age+
                                   SSS_faml+ relationship_status+ secondary*occupation+
                                   (1|residing_country),
                                 data=data.filtered, family = gaussian(),
                                 cores=4,chains=4, save_pars = save_pars(all = T),
                                 sample_prior ='yes', seed=1660415,prior=prior.coef)
res.edu.occu.1 <- bayes_factor(res.int.edu.occu,res.1)

#0

# frequentist and icc
freq.res.1 <- lmer(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),data=data.filtered)
freq.res.int <- lmer(resilience ~ secondary*gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1|residing_country),data=data.filtered)

freq.res.00 <- lm(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status,data=data.filtered)

EMAtools::lme.dscore(freq.res.1,data=data.filtered,'lme4')
conf.res <- confint(freq.res.1)
r.squaredGLMM(freq.res.1)

# icc
performance::icc(freq.res.1)

save.image(file='Vaccine_H3bs.RData')
