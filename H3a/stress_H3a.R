# hypothesis testing

library(lmerTest)
library(brms)
library(EMAtools)
library(sjstats)

# H1

# load dataset
load('../Stress_aligned.RData')

# H1 pss / resilience ~ primary

# standardize
data.filtered$primary_stressor_avg <- scale(data.filtered$primary_stressor_avg)
data.filtered$pss <- scale(data.filtered$pss)
data.filtered$resilience <- scale(data.filtered$resilience)

# pss
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
pss.1 <- brms::brm(pss ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model 2
pss.2 <- brms::brm(pss ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+primary_stressor_avg|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model comparison
pss.10 <- bayes_factor(pss.1,pss.0)
pss.20 <- bayes_factor(pss.2,pss.0)
pss.21 <- bayes_factor(pss.2,pss.1)
log(pss.21$bf[1])

# hypothesis testing
hypothesis(pss.2,'primary_stressor_avg > 0')

# effect size
lmer.pss.1 <- lmer (pss ~ primary_stressor_avg + gender + education + work_location + age+
                      SSS_faml+ relationship_status+
                      (1|residing_country),
                    data=data.filtered)
EMAtools::lme.dscore(lmer.pss.1,data.filtered,'lme4')

# ICC calculation
lmer.pss.2 <- lmer (pss ~ primary_stressor_avg + gender + education + work_location + age+
                      SSS_faml+ relationship_status+
                      (1+primary_stressor_avg|residing_country),
                    data=data.filtered)
performance::icc(lmer.pss.2)

library(nlme)
pss.conf <- confint(lmer.pss.2)

######
# resilience

# model 0
res.0 <- brms::brm(resilience ~ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415)

# model 1
res.1 <- brms::brm(resilience ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model 2
res.2 <- brms::brm(resilience ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+primary_stressor_avg|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# model comparison
res.10 <- bayes_factor(res.1,res.0)
log(res.10$bf[1])
res.20 <- bayes_factor(res.2,res.0)
log(res.20$bf[1])
res.21 <- bayes_factor(res.2,res.1)
res.21

# model 2 slightly better than res.1 -> inconclusive
# hypothesis testing (both)

hypothesis(res.1,'primary_stressor_avg < 0')
hypothesis(res.2,'primary_stressor_avg < 0')

# effect size
lmer.res.1 <- lmer (resilience ~ primary_stressor_avg + gender + education + work_location + age+
                      SSS_faml+ relationship_status+
                      (1|residing_country),
                    data=data.filtered)
EMAtools::lme.dscore(lmer.res.1,data.filtered,'lme4')

# icc test
lmer.res.2 <- lmer (resilience ~ primary_stressor_avg + gender + education + work_location + age+
                      SSS_faml+ relationship_status+
                      (1+primary_stressor_avg|residing_country),
                    data=data.filtered)
performance::icc(lmer.res.2)

res.conf <- confint(lmer.res.2)

# result save
save.image(file='stress_H3a.RData')

# R2s
library(MuMIn)

# pss
MuMIn::r.squaredGLMM(lmer.pss.2)
MuMIn::r.squaredGLMM(lmer.res.2)

# vif
library(car)
vif(lmer.pss.2)
vif(lmer.res.2)
