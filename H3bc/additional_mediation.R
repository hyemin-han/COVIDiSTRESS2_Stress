# additional mediation
# gender/SES -> secondary -> pss/brms

library(lmerTest)
library(brms)
library(EMAtools)
library(sjstats)

load('../Stress_aligned.RData')

# interactions
data.filtered$secondary <- scale(data.filtered$secondary)
data.filtered$pss <- scale(data.filtered$pss)
data.filtered$resilience <- scale(data.filtered$resilience)

# being male
# exclude the last
data.filtered <- data.filtered[data.filtered$gender != 'Other/Would rather not say',]
data.filtered$man <- data.filtered$gender=='Male'

# pss 
options(width = 2000)

# prior -> normal distribution 
prior.coef <- brms::prior(cauchy(0.,1),class='b')

# null model
model.mediator.0 <- bf(secondary ~ 
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
model.mediator.1 <- bf(secondary ~man+
                         SSS_faml+ 
                         (1|residing_country))
model.pss.1 <- bf(pss ~ secondary+ man+
                    SSS_faml+ 
                    (1|residing_country))
med.pss.1 = brm(
  model.mediator.1 + model.pss.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef, iter=4000
)


# random slope
model.mediator.2 <- bf(secondary ~man+
                         SSS_faml+ 
                         (1+man+SSS_faml|residing_country))
model.pss.2 <- bf(pss ~ secondary+ man+
                    SSS_faml+ 
                    (1+secondary+man+SSS_faml|residing_country))
med.pss.2 = brm(
  model.mediator.1 + model.pss.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef, iter=4000
)

# bayes factor
bf.pss.10 <- bayes_factor(med.pss.1,med.pss.0, log=T) # 1961.45553
bf.pss.20 <- bayes_factor(med.pss.2,med.pss.0,log = T) #1961.47915
bf.pss.21 <- bf.pss.20$bf - bf.pss.10$bf #0.02362304
# no substantial difference, but 2 is very slightly better.

bayestestR::mediation(med.pss.2)
'
# Causal Mediation Analysis for Stan Model

  Treatment: manTRUE
  Mediator : secondary
  Response : pss

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |   -0.196 | [-0.234, -0.157]
Indirect Effect (ACME) |   -0.047 | [-0.063, -0.032]
Mediator Effect        |    0.355 | [ 0.335,  0.374]
Total Effect           |   -0.243 | [-0.284, -0.202]

Proportion mediated: 19.34% [13.25%, 25.43%]
'
bayestestR::mediation(med.pss.2,treatment='SSS_faml')
'
# Causal Mediation Analysis for Stan Model

  Treatment: SSS_faml
  Mediator : secondary
  Response : pss

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |   -0.101 | [-0.113, -0.089]
Indirect Effect (ACME) |   -0.062 | [-0.067, -0.056]
Mediator Effect        |    0.355 | [ 0.335,  0.374]
Total Effect           |   -0.163 | [-0.175, -0.151]

Proportion mediated: 37.94% [34.23%, 41.65%]'

# path significance

# male -> secondary
hypothesis(med.pss.2,'secondary_manTRUE<0')
# ses -> secondary
hypothesis(med.pss.2,'secondary_SSS_faml<0')
#  male ->  pss
hypothesis(med.pss.2,'pss_manTRUE<0')
# ses -> pss
hypothesis(med.pss.2,'pss_SSS_faml<0')

# all significant
# partial mediation

save.image('additional_mediation_demo.RData')


### BRS

# null model

model.res.0 <- bf(resilience ~ 
                    (1|residing_country)
)
med.res.0 = brm(
  model.mediator.0 + model.res.0 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415, iter=4000
)

# random intercept

model.res.1 <- bf(resilience ~ secondary+ man+
                    SSS_faml+ 
                    (1|residing_country))
med.res.1 = brm(
  model.mediator.1 + model.res.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef, iter=4000
)


# random slope

model.res.2 <- bf(resilience ~ secondary+ man+
                    SSS_faml+ 
                    (1+secondary+man+SSS_faml|residing_country))
med.res.2 = brm(
  model.mediator.1 + model.res.1 + set_rescor(F),
  data=data.filtered,
  family = gaussian(),
  cores=4,chains=4, save_pars = save_pars(all = T),
  sample_prior ='yes', seed=1660415,prior=prior.coef, iter=4000
)

# bayes factor
bf.res.10 <- bayes_factor(med.res.1,med.res.0, log=T) # 1192.11018
bf.res.20 <- bayes_factor(med.res.2,med.res.0,log = T) #1192.15034
bf.res.21 <- bf.res.20$bf - bf.res.10$bf #0.04016309

# model 2 very slightly better


bayestestR::mediation(med.res.2,iv = 'manTRUE')
'
# Causal Mediation Analysis for Stan Model

  Treatment: manTRUE
  Mediator : secondary
  Response : resilience

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |    0.200 | [ 0.151,  0.248]
Indirect Effect (ACME) |    0.030 | [ 0.020,  0.042]
Mediator Effect        |   -0.226 | [-0.249, -0.202]
Total Effect           |    0.230 | [ 0.181,  0.280]

Proportion mediated: 13.20% [8.04%, 18.36%]
'
bayestestR::mediation(med.res.2,treatment = 'SSS_faml')
'
# Causal Mediation Analysis for Stan Model

  Treatment: SSS_faml
  Mediator : secondary
  Response : resilience

Effect                 | Estimate |          95% ETI
----------------------------------------------------
Direct Effect (ADE)    |    0.097 | [ 0.084,  0.111]
Indirect Effect (ACME) |    0.038 | [ 0.034,  0.044]
Mediator Effect        |   -0.226 | [-0.249, -0.202]
Total Effect           |    0.136 | [ 0.122,  0.149]

Proportion mediated: 28.27% [23.98%, 32.57%]

'

# path significance

# male -> secondary
hypothesis(med.res.2,'secondary_manTRUE<0')
# ses -> secondary
hypothesis(med.res.2,'secondary_SSS_faml<0')
#  male ->  brs
hypothesis(med.res.2,'resilience_manTRUE>0')
# ses -> brs
hypothesis(med.res.2,'resilience_SSS_faml>0')

# all significant
# partial mediation

save.image('additional_mediation_demo.RData')
