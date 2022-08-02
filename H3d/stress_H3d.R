# hypothesis testing

library(lmerTest)
library(brms)
library(EMAtools)
library(bayestestR)

load('../Stress_aligned.RData')

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

#####
# 1. PSS

# null model
model.mediator.0 <- bf(sps ~ 
                     (1|residing_country))
model.pss.0 <- bf(pss ~ 
                         (1|residing_country)
                       )

med.pss.0 = brm(
  model.mediator.0 + model.pss.0 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415, iter=4000
)


# random intercept
model.mediator.1 <- bf(sps ~ identity+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1|residing_country))
model.pss.1 <- bf(pss ~ identity+ sps+gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1|residing_country))

med.pss.1 = brm(
  model.mediator.1 + model.pss.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# random slope
model.mediator.2 <- bf(sps ~ identity+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1+identity|residing_country))
model.pss.2 <- bf(pss ~ identity+ sps+gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1+identity+sps|residing_country))

med.pss.2 = brm(
  model.mediator.2 + model.pss.2 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# bayes factor
bf.pss.10 <- bayes_factor(med.pss.1,med.pss.0, log=T)
bf.pss.20 <- bayes_factor(med.pss.2,med.pss.0,log = T)
bf.pss.21 <- bf.pss.20$bf - bf.pss.10$bf
# random slope model is the best
# mediation analysis
bayestestR::mediation(med.pss.2)
'
# Causal Mediation Analysis for Stan Model

  Treatment: identity
  Mediator : sps
  Response : pss

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |   -0.050 | [-0.083, -0.020]
Indirect Effect (ACME) |   -0.061 | [-0.072, -0.051]
Mediator Effect        |   -0.311 | [-0.346, -0.275]
Total Effect           |   -0.111 | [-0.147, -0.079]

Proportion mediated: 54.88% [38.39%, 71.36%]
'
# each path bayes test
hypothesis(med.pss.2, 'sps_identity>0') # Inf
hypothesis(med.pss.2, 'pss_identity<0') # 1999
hypothesis(med.pss.2, 'pss_sps<0') # Inf

save.image('stress_H3d.RData')



####
# 2. resilience


# null model

model.res.pss.0 <- bf(resilience ~ 
                    (1|residing_country)
)

med.res.0 = brm(
  model.mediator.0 + model.res.pss.0 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415, iter=2000
)


# random intercept

model.res.1 <- bf(resilience ~ identity+ sps+gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1|residing_country))

med.res.1 = brm(
  model.mediator.1 + model.res.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# random slope

model.res.2 <- bf(resilience ~ identity+ sps+gender + education + work_location + age+
                    SSS_faml+ relationship_status+
                    (1+identity+sps|residing_country))

med.res.2 = brm(
  model.mediator.2 + model.res.2 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef
)

# bayes factor
bf.res.10 <- bayes_factor(med.res.1,med.res.0, log=T)
bf.res.20 <- bayes_factor(med.res.2,med.res.0,log = T, maxiter=5000)
bf.res.21 <- bf.res.20$bf - bf.res.10$bf
# .7744457. Inconclusive. However, random slope model is marginally better.


# mediation analysis
bayestestR::mediation(med.res.2)
'
# Causal Mediation Analysis for Stan Model

  Treatment: identity
  Mediator : sps
  Response : resilience

Effect                 | Estimate |        95% ETI
--------------------------------------------------
Direct Effect (ADE)    |    0.096 | [0.066, 0.129]
Indirect Effect (ACME) |    0.049 | [0.039, 0.061]
Mediator Effect        |    0.252 | [0.210, 0.292]
Total Effect           |    0.146 | [0.114, 0.180]

Proportion mediated: 33.80% [24.66%, 42.93%]
'
# just in case
bayestestR::mediation(med.res.1)
'
# Causal Mediation Analysis for Stan Model

  Treatment: identity
  Mediator : sps
  Response : resilience

Effect                 | Estimate |        95% ETI
--------------------------------------------------
Direct Effect (ADE)    |    0.082 | [0.059, 0.105]
Indirect Effect (ACME) |    0.050 | [0.043, 0.058]
Mediator Effect        |    0.256 | [0.233, 0.279]
Total Effect           |    0.132 | [0.108, 0.156]

Proportion mediated: 38.02% [29.81%, 46.24%]
' # similar

# path bayes factor
hypothesis(med.res.2, 'sps_identity>0') # Inf
hypothesis(med.res.2, 'resilience_identity>0') # Inf
hypothesis(med.res.2, 'resilience_sps>0') # Inf

save.image('stress_H3d.RData')
