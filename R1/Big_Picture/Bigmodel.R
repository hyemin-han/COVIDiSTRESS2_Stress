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
bf.pss.2.0 <- bayes_factor(pss.2,pss.0,log=T)
bf.pss.2.1 <- bayes_factor(pss.2,pss.1,log=T)

#Estimated log Bayes factor in favor of pss.2 over pss.1: 10.43362

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

# fixe only
model.mediator.f <- bf(sps ~  identity+gender + education + work_location + age+
                         SSS_faml+ relationship_status)
model.pss.f <- bf(pss ~ primary_stressor_avg + primary_stressor_avg*sps +
                    secondary*identity + secondary*sps +
                    gender + education + work_location + age+
                    SSS_faml+ relationship_status)
pss.f <- brms::brm( model.pss.f+model.mediator.f+set_rescor(F) ,
                    data=data.filtered, family = gaussian(),
                    cores=4,chains=4, save_pars = save_pars(all = T),
                    sample_prior ='yes', seed=1660415,prior=prior.coef)

# bfs
bf.pss.1.0 <- bayes_factor(pss.1,pss.0,log=T)
bf.pss.2.0 <- bayes_factor(pss.2,pss.0,log=T)
bf.pss.2.f <- bayes_factor(pss.2,pss.f,log=T) #486.67946

# hypothesis testing
# primary
hypothesis(pss.2, 'pss_primary_stressor_avg >0') #0.14      0.02     0.11     0.17        Inf         1    *
hypothesis(pss.2, 'pss_primary_stressor_avg:sps =0') #0.01      0.01    -0.01     0.03      40.98      0.99    
# secondary
hypothesis(pss.2, 'pss_secondary >0') #0.26      0.01     0.24     0.29        Inf         1    *
hypothesis(pss.2, 'pss_sps:secondary=0 ') #0.01      0.01    -0.01     0.03     79.51      0.99 
# group id
hypothesis(pss.2 ,'pss_identity < 0') #-0.03      0.01    -0.06    -0.01     221.22         1    *
#hypothesis(pss.2,'pss_primary_stressor_avg:identity=0') #0      0.01    -0.02     0.02      129.3      0.99  
hypothesis(pss.2,'pss_secondary:identity=0') #-0.01      0.01    -0.03     0.01     68.19      0.99   
# sps
hypothesis(pss.2,'pss_sps<0') #-0.26      0.02    -0.28    -0.23        Inf         1    *
hypothesis(pss.2,'sps_identity>0') #0.21      0.01     0.19     0.24        Inf         1    *

mediation(pss.2)
'
 Treatment: identity
  Mediator : sps
  Response : pss

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |   -0.034 | [-0.063, -0.008]
Indirect Effect (ACME) |   -0.054 | [-0.065, -0.044]
Mediator Effect        |   -0.256 | [-0.290, -0.221]
Total Effect           |   -0.088 | [-0.117, -0.059]

Proportion mediated: 60.87% [40.13%, 81.61%]
'

# so far.
save.image('big.pss.RData')

