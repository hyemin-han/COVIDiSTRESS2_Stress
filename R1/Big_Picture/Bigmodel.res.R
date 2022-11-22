# all in one
library(psych)
library(BayesFactor)
library(doMC)
library(lmerTest)
library(brms)
library(bayestestR)

# outcome ~ primary + secondary + groupid + social support + demo + (interaction)
# social support ~ groupid + demo
# group id * secondary only

load('../../Stress_aligned.RData')
data.filtered$pss <- data.filtered$resilience

# resilience

# standardize
data.filtered$primary_stressor_avg <- scale(data.filtered$primary_stressor_avg)
data.filtered$secondary <- scale(data.filtered$secondary)
data.filtered$pss <- scale(data.filtered$pss)
data.filtered$sps <- scale(data.filtered$sps)
data.filtered$resilience <- scale(data.filtered$resilience)
data.filtered$identity <- scale(data.filtered$identity)
data.filtered$age <- scale(data.filtered$age)

# pss ~ secondary
#options(width = 2000)

# prior -> normal distribution 
prior.coef <- brms::prior(cauchy(0.,1),class='b')

# model definition
model.mediator.1 <- bf(sps ~ identity+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1|residing_country))
model.pss.1 <- bf(pss ~ primary_stressor_avg + primary_stressor_avg*sps +
                    secondary*identity + secondary*sps +
                    gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1|residing_country))
pss.1 <- brms::brm(model.mediator.1 + model.pss.1 + set_rescor(F),
                        data=data.filtered, family = gaussian(),
                        cores=4,chains=4, save_pars = save_pars(all = T),
                        sample_prior ='yes', seed=1660415,prior=prior.coef)

# slope
model.mediator.2 <- bf(sps ~ identity+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1+identity|residing_country))
model.pss.2 <- bf(pss ~ primary_stressor_avg + primary_stressor_avg*sps +
                    secondary*identity + secondary*sps +
                    gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1+primary_stressor_avg+secondary+
                       identity+sps|residing_country))
pss.2 <- brms::brm(model.mediator.2 + model.pss.2 + set_rescor(F),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# bf
bf.pss.2.1 <- bayes_factor(pss.2,pss.1,log=T)
#Estimated log Bayes factor in favor of pss.2 over pss.1: -12.81481

# two other models: null, fixed only

model.mediator.0 <- bf(sps ~  gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1|residing_country))
model.pss.0 <- bf(pss ~ 
                    gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1|residing_country))
pss.0 <- brms::brm( model.pss.0+model.mediator.0+set_rescor(F) ,
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# fixed effects only
model.mediator.f <- bf(sps ~ identity+ gender + education + work_location + age+
                         SSS_faml+ relationship_status
                         )
model.pss.f <- bf(pss ~ primary_stressor_avg + primary_stressor_avg*sps +
                    secondary*identity + secondary*sps +
                    gender + education + work_location + age+
                    SSS_faml+ relationship_status)
pss.f <- brms::brm(model.mediator.f + model.pss.f + set_rescor(F),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)

# bfs
bf.pss.1.0 <- bayes_factor(pss.1,pss.0,log=T) # 4009.51726

bf.pss.2.0 <- bayes_factor(pss.2,pss.0,log=T)
bf.pss.1.2 <- bayes_factor(pss.1,pss.2,log=T) # 12.21
bf.pss.1.f<- bayes_factor(pss.1,pss.f,log=T) # 220.50337


# hypothesis testing
# primary
hypothesis(pss.1, 'pss_primary_stressor_avg <0') #-0.09      0.01    -0.11    -0.07        Inf         1    *
hypothesis(pss.1, 'pss_primary_stressor_avg:sps =0') # -0.02      0.01    -0.04     0.01      44.34      0.98    
# secondary
hypothesis(pss.1, 'pss_secondary <0') #-0.15      0.01    -0.17    -0.13        Inf         1    *
hypothesis(pss.1, 'pss_sps:secondary=0 ') #-0.03      0.01    -0.05    -0.01       4.75      0.83    * #NOTE 
# group id 
hypothesis(pss.1 ,'pss_identity > 0') #0.08      0.01     0.06      0.1        Inf         1    *

hypothesis(pss.1,'pss_secondary:identity=0') #0.02      0.01    -0.01     0.04       49.9      0.98   
# sps
hypothesis(pss.1,'pss_sps>0') #0.24      0.01     0.22     0.26        Inf         1    *
hypothesis(pss.1,'sps_identity>0')

mediation(pss.2)
'
  Treatment: identity
  Mediator : sps
  Response : pss

Effect                 | Estimate |        95% ETI
--------------------------------------------------
Direct Effect (ADE)    |    0.087 | [0.056, 0.123]
Indirect Effect (ACME) |    0.048 | [0.037, 0.060]
Mediator Effect        |    0.223 | [0.182, 0.260]
Total Effect           |    0.135 | [0.102, 0.173]

Proportion mediated: 35.30% [24.66%, 45.94%]
'

# so far.
save.image('big.res.RData')

