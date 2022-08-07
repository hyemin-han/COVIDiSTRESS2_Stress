# hypothesis testing

library(lmerTest)
library(brms)
library(EMAtools)

load('Stress_aligned.RData')

# interactions
data.filtered$secondary <- scale(data.filtered$secondary)
data.filtered$pss <- scale(data.filtered$pss)
data.filtered$resilience <- scale(data.filtered$resilience)
data.filtered$identity <- scale(data.filtered$identity)
data.filtered$sps<-scale(data.filtered$sps)

# pss ~ secondary
options(width = 2000)

# prior -> normal distribution 
prior.coef <- brms::prior(cauchy(0.,1),class='b')

### PSS : vs. model 2
# model 2
pss.2 <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# by identity
pss.identity <- brms::brm(pss ~ secondary*identity+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary+identity|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# only main effect
pss.main.identity <- brms::brm(pss ~ secondary+identity+ gender + education + work_location + age+
                                 SSS_faml+ relationship_status+
                                 (1+secondary+identity|residing_country),
                               data=data.filtered, family = gaussian(),
                               cores=4,chains=4, save_pars = save_pars(all = T),
                               sample_prior ='yes', seed=1660415,prior=prior.coef)


# bfs
bf.int <- bayes_factor(pss.identity,pss.2)
bf.main <- bayes_factor(pss.main.identity,pss.2)

log(bf.int$bf[1])
log(bf.main$bf[1])


# main effect test
hypothesis (pss.identity, 'secondary > 0') # Inf
hypothesis (pss.identity, 'identity < 0') # Inf
hypothesis (pss.identity, 'secondary:identity < 0')

# frequentist
freq.pss.int <- lmer(pss ~ secondary*identity+ gender + education + work_location + age+
                       SSS_faml+ relationship_status+
                       (1+secondary+identity|residing_country),data.filtered)
conf.pss <- confint(freq.pss.int)

# R2
library(MuMIn)

r.squaredGLMM(freq.pss.int)

freq.pss.int.1 <- lmer(pss ~ secondary*identity+ gender + education + work_location + age+
                       SSS_faml+ relationship_status+
                       (1|residing_country),data.filtered)
lme.dscore(freq.pss.int.1,data.filtered,'lme4')

### H3d

## H3da
# outcomes ~ id

# pss
pss.id.0 <- brms::brm(pss ~ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)

pss.id.1 <- brms::brm(pss ~ identity+ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)

pss.id.2 <- brms::brm(pss ~ identity+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+identity|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
# bfs
bf.pss.id.10 <- bayes_factor(pss.id.1,pss.id.0)
bf.pss.id.20 <- bayes_factor(pss.id.2,pss.id.0)
bf.pss.id.21 <- bayes_factor(pss.id.2,pss.id.1)

# 2 best
hypothesis(pss.id.2,'identity<0')

performance::icc(freq.pss.int)


# only with the main effect of stressor


# res
res.id.0 <- brms::brm(resilience ~  gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)
res.id.1 <- brms::brm(resilience ~ identity+ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)

res.id.2 <- brms::brm(resilience ~ identity+ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1+identity|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)


res.2 <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1+secondary|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)
res.int <- brms::brm(resilience ~ secondary*identity+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary+identity|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

res.1 <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
res.int.1<- brms::brm(resilience ~ secondary*identity+ gender + education + work_location + age+
                        SSS_faml+ relationship_status+
                        (1|residing_country),
                      data=data.filtered, family = gaussian(),
                      cores=4,chains=4, save_pars = save_pars(all = T),
                      sample_prior ='yes', seed=1660415,prior=prior.coef)

# bfs
bf.res.int <- bayes_factor(res.int,res.2,log=T)
bf.res.int1 <- bayes_factor(res.int.1,res.1,log=T)

# hypothesis
hypothesis(res.int.1, 'secondary < 0')
hypothesis(res.int.1, 'identity > 0')
hypothesis(res.int.1, 'secondary:identity > 0')

# bfs
#bf.res.id.10 <- bayes_factor(res.id.1,res.id.0)
#bf.res.id.20 <- bayes_factor(res.id.2,res.id.0)
#bf.res.id.21 <- bayes_factor(res.id.2,res.id.1)

# 1 best
#hypothesis(res.id.1,'identity>0')


# frequentist
freq.res.int <- lmer(resilience ~ secondary*identity+ gender + education + work_location + age+
  SSS_faml+ relationship_status+
  (1+secondary+identity|residing_country),data.filtered)
freq.res.int1 <- lmer(resilience ~ secondary*identity+ gender + education + work_location + age+
                       SSS_faml+ relationship_status+
                       (1|residing_country),data.filtered)
lme.dscore(freq.res.int1,data.filtered,'lme4')
conf.res <- confint(freq.res.int1)

r.squaredGLMM(freq.res.int1)
performance::icc(freq.res.int1)

## H3db

# mediation
model.mediator.1 <- bf(sps ~ identity + gender + education + work_location + age+
                         SSS_faml+ relationship_status+(1|residing_country))
model.pss.1 <- bf(pss ~ sps+identity + gender + education + work_location + age+
                    SSS_faml+ relationship_status+(1|residing_country))
med.pss.1 = brm(
  model.mediator.1 + model.pss.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# with random slopes
model.mediator.2 <- bf(sps ~ identity + gender + education + work_location + age+
                         SSS_faml+ relationship_status+(1+identity|residing_country))
model.pss.2 <- bf(pss ~ sps+identity + gender + education + work_location + age+
                    SSS_faml+ relationship_status+(1+identity+sps|residing_country))
med.pss.2 = brm(
  model.mediator.2 + model.pss.2 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# bf
bf.med.pss.21 <- bayes_factor(med.pss.2,med.pss.1)
# 2 better
# mediation
library(bayestestR)
bayestestR::mediation(med.pss.2)


### res

model.res.1 <- bf(resilience ~ sps+identity + gender + education + work_location + age+
                    SSS_faml+ relationship_status+(1|residing_country))
med.res.1 = brm(
  model.mediator.1 + model.res.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# with random slopes

model.res.2 <- bf(resilience ~ sps+identity + gender + education + work_location + age+
                    SSS_faml+ relationship_status+(1+identity+sps|residing_country))
med.res.2 = brm(
  model.mediator.2 + model.res.2 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# bf
bf.med.res.21 <- bayes_factor(med.res.2,med.res.1)
# 1 better
# mediation
library(bayestestR)
bayestestR::mediation(med.res.1)

save.image('Vaccine_H3_groupid.RData')

# VIF check
library(car)
pss.vif <- car::vif(freq.pss.int)
res.vif <- car::vif(freq.res.int1)
library(effectsize)
effectsize::interpret_vif(pss.vif[,1])
effectsize::interpret_vif(res.vif[,1])
