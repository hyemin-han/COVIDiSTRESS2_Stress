# hypothesis testing

library(lmerTest)
library(brms)
library(EMAtools)

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

# only with gender
pss.int.gender <- brms::brm(pss ~ gender*secondary + work_location + age+
                              SSS_faml+ relationship_status+ education+
                              (1+secondary|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.gender <- bayes_factor(pss.int.gender,pss.2)

# only with education
pss.int.edu <- brms::brm(pss ~ gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ 
                       (1+secondary|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.edu <- bayes_factor(pss.int.edu,pss.2)

# only with occupation
pss.int.ocu <- brms::brm(pss ~gender +education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1+secondary|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef, iter= 5000)

bf.ocu <- bayes_factor(pss.int.ocu,pss.2)

# gender / education
pss.int.gender.edu <- brms::brm(pss ~ gender*secondary + work_location + age+
                              SSS_faml+ relationship_status+ secondary*education+
                              (1+secondary|residing_country),
                            data=data.filtered, family = gaussian(),
                            cores=4,chains=4, save_pars = save_pars(all = T),
                            sample_prior ='yes', seed=1660415,prior=prior.coef)

bf.gender.edu <- bayes_factor(pss.int.edu,pss.2)

# gender / occupation

# education / occupation


# double check with freq
freq.pss.2 <- lmer(pss ~ secondary+ work_location + age+
                       SSS_faml+ relationship_status+ education+
                       (1+secondary|residing_country),
                     data=data.filtered)
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

# int
res.int <- brms::brm(resilience ~ secondary*gender + secondary*education + work_location + age+
                       SSS_faml+ relationship_status+ secondary*occupation+
                       (1|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)

# gender
res.int.gender <- brms::brm(resilience ~ secondary*gender + education + work_location + age+
                       SSS_faml+ relationship_status+ 
                       (1|residing_country),
                     data=data.filtered, family = gaussian(),
                     cores=4,chains=4, save_pars = save_pars(all = T),
                     sample_prior ='yes', seed=1660415,prior=prior.coef)
# bayes factors
res.int.all <- bayes_factor(res.int,res.1)
bf.res.gender <- bayes_factor(res.int.gender,res.1)

save.image(file='Vaccine_H3bs.RData')
