# hypothesis testing

library(lmerTest)
library(brms)

# H3a

# load dataset
load('Stress_aligned.RData')

# do regression
reg.h3a.pss <- lm(data.aligned$pss ~ data.aligned$primary_stressor_avg) 
reg.h3a.rs <- lm(data.aligned$resilience ~ data.aligned$primary_stressor_avg)

# do mlm
lreg.h3a.pss.0 <- lmer(pss~(1|residing_country), data=data.aligned)
lreg.h3a.pss.1 <- lmer(pss~primary_stressor_avg+(1|residing_country), data=data.aligned)
lreg.h3a.pss.2 <- lmer(pss~primary_stressor_avg+(1+primary_stressor_avg|
                                                   residing_country), data=data.aligned)
AIC(lreg.h3a.pss.0,lreg.h3a.pss.1,lreg.h3a.pss.2)
BIC(lreg.h3a.pss.0,lreg.h3a.pss.1,lreg.h3a.pss.2)

# brms
prior.coef <- brms::prior(cauchy(0.,1),class='b')
breg.h3a.pss.2 <- brm(pss~primary_stressor_avg+(1+primary_stressor_avg|
                                                   residing_country), data=data.aligned,
                      family = gaussian(),
                      cores=4,chains=4, save_all_pars = TRUE,
                      sample_prior ='yes',prior=prior.coef, seed=1660415)



### mediation pilot
model.3da.pss.2 <- bf(sps~identity+(1+identity|
                                      residing_country))
model.3da.pss.m.2 <- bf(pss~identity+sps+(1+identity+sps|
                                            residing_country))


h3da.pss.test <- brm(model.3da.pss.m.2 + model.3da.pss.2 + 
                   set_rescor(F),data=data.aligned,family = gaussian(),
                 cores=4,chains=4, save_all_pars = TRUE,
                 sample_prior ='yes',prior=prior.coef, seed=1660415)
bayestestR::mediation(h3da.pss.test)


model.3da.rs.m.2 <- bf(resilience~identity+sps+(1+identity+sps|
                                            residing_country))


h3da.rs.test <- brm(model.3da.rs.m.2 + model.3da.pss.2 + 
                   set_rescor(F),data=data.aligned,family = gaussian(),
                 cores=4,chains=4, save_all_pars = TRUE,
                 sample_prior ='yes',prior=prior.coef, seed=1660415)
bayestestR::mediation(h3da.rs.test)
