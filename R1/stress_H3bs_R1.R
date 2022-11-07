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
                       SSS_faml+ relationship_status+ education+ gender+
                       (1+secondary|residing_country),
                     data=data.filtered)
conft.pss.2<-confint(freq.pss.2)

# mumin
library(MuMIn)
MuMIn::r.squaredGLMM(freq.pss.2)

freq.pss.1 <- lmer(pss ~ secondary+ work_location + age+
                     SSS_faml+ relationship_status+ education+gender+
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
res.21 <- bayes_factor(res.2,res.1) # -2.92144

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

# linearity
library(ggfortify)

test.pss <- lm(pss ~ secondary*identity+ gender + education + work_location + age+
                 SSS_faml+ relationship_status, data.filtered)
autoplot(test.pss, label=F)

test.res <- lm(resilience ~ secondary*identity+ gender + education + work_location + age+
                 SSS_faml+ relationship_status, data.filtered)
autoplot(test.res, label=F)


# each demographics. Main effect test from pss.2
hypothesis(pss.2,'genderMale<0')
hypothesis(pss.2,'SSS_faml<0')
hypothesis(res.1,'genderMale>0')
hypothesis(res.1,'SSS_faml>0')

# frequentist
summary(freq.pss.2)
summary(freq.res.1)

# bayesian test of secondary stressor ~ gender / ses
secondary.gender.0 <- brm(secondary~(1|residing_country),
                           data=data.filtered, family = gaussian(),
                           cores=4,chains=4, save_pars = save_pars(all = T),
                           sample_prior ='yes', seed=1660415)
secondary.gender.1<- brm(secondary~gender+(1|residing_country),
                          data=data.filtered, family = gaussian(),
                          cores=4,chains=4, save_pars = save_pars(all = T),
                          sample_prior ='yes', seed=1660415,prior=prior.coef)

secondary.ses.1 <- brm(secondary~SSS_faml+(1|residing_country),
                       data=data.filtered, family = gaussian(),
                       cores=4,chains=4, save_pars = save_pars(all = T),
                       sample_prior ='yes', seed=1660415,prior=prior.coef)
secondary.both.1 <-  brm(secondary~gender+SSS_faml+(1|residing_country),
                         data=data.filtered, family = gaussian(),
                         cores=4,chains=4, save_pars = save_pars(all = T),
                         sample_prior ='yes', seed=1660415,prior=prior.coef)

# comparison
bf.0.gender <- bayes_factor(secondary.gender.1,secondary.gender.0,log=T) #709.58980
bf.0.ses <- bayes_factor(secondary.ses.1,secondary.gender.0,log=T) # 21.46429
bf.0.both <- bayes_factor(secondary.both.1,secondary.gender.0,log=T) # 729.55055

# both best
hypothesis(secondary.both.1,'SSS_faml<0')
hypothesis(secondary.both.1,'genderMale<0')


#### no random
# PSS
pss.n <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status,
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
# comparison
pss.n0 <- bayes_factor(pss.n,pss.0,log=T                       ) # 1959.67492
pss.2n <- bayes_factor(pss.2, pss.n, log=T) # 217.40655
pss.1n <- bayes_factor(pss.1, pss.n, log=T) # 199.23953

# still better
# testing
hypothesis(pss.n,'secondary>0')
'
Hypothesis Tests for class b:
       Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
1 (secondary) > 0     0.36      0.01     0.34     0.37        Inf         1    *' # same

## RES
res.n <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status,
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
# comparison
res.n0 <- bayes_factor(res.n,res.0, log=T) # 1382.48734
res.1n <- bayes_factor(res.1,res.n, log=T) # 167.22297
res.2n <- bayes_factor(res.2,res.n, log=T) # 164.46898

# testing
hypothesis(res.n,'secondary<0')
'
Hypothesis Tests for class b:
       Hypothesis Estimate Est.Error CI.Lower CI.Upper Evid.Ratio Post.Prob Star
1 (secondary) < 0     -0.2      0.01    -0.22    -0.18        Inf         1    *' # same



#### sensivity check
prior.coef1 <- brms::prior(normal(0.,1000000),class='b')

pss.2.norm <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef1)

hypothesis(pss.2.norm,'secondary>0') # Inf
hypothesis(pss.2.norm,'genderMale < 0') # Inf
hypothesis(pss.2.norm,'SSS_faml<0') # Inf

res.2.norm <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef1)

hypothesis(res.2.norm,'secondary<0') # Inf
hypothesis(res.2.norm,'genderMale > 0') # Inf
hypothesis(res.2.norm,'SSS_faml>0') # Inf

# all same