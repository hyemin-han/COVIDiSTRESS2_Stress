# additional

# stressors * social support


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
data.filtered$secondary <- scale(data.filtered$secondary)
data.filtered$sps<-scale(data.filtered$sps)

# pss
# prior -> normal distribution 
prior.coef <- brms::prior(cauchy(0.,1),class='b')

### 1
# PSS ~ primary * sps
# compare with the best model
pss.2 <- brms::brm(pss ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+primary_stressor_avg|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
pss.2.sps <- brms::brm(pss ~ primary_stressor_avg*sps + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+primary_stressor_avg+sps|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
# bf comparison
bf.pss.2.sps <- bayes_factor(pss.2.sps,pss.2,log=T) #498.93396
# test
hypothesis(pss.2.sps, 'primary_stressor_avg>0') # Inf
hypothesis(pss.2.sps, 'sps <0') # Inf
hypothesis(pss.2.sps,'primary_stressor_avg:sps=0')
# only main effects mattered

### 2
# RES ~ primary*sps
res.2 <- brms::brm(resilience ~ primary_stressor_avg + gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+primary_stressor_avg|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
res.2.sps <- brms::brm(resilience ~ primary_stressor_avg*sps + gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1+primary_stressor_avg+sps|residing_country),
                       data=data.filtered, family = gaussian(),
                       cores=4,chains=4, save_pars = save_pars(all = T),
                       sample_prior ='yes', seed=1660415,prior=prior.coef)
# bf comparison
bf.res.2.sps <- bayes_factor(res.2.sps,res.2,log=T) # 268.21327
# better
hypothesis (res.2.sps, 'primary_stressor_avg < 0') # Inf
hypothesis (res.2.sps, 'sps >0') # Inf
hypothesis(res.2.sps,' primary_stressor_avg:sps=0')




### 3
# PSS ~ sec * sps

pss.sec.2 <- brms::brm(pss ~ secondary+ gender + education + work_location + age+
                     SSS_faml+ relationship_status+
                     (1+secondary|residing_country),
                   data=data.filtered, family = gaussian(),
                   cores=4,chains=4, save_pars = save_pars(all = T),
                   sample_prior ='yes', seed=1660415,prior=prior.coef)
pss.sec.2.sps <- brms::brm(pss ~ secondary*sps+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1+secondary+sps|residing_country),
                       data=data.filtered, family = gaussian(),
                       cores=4,chains=4, save_pars = save_pars(all = T),
                       sample_prior ='yes', seed=1660415,prior=prior.coef)
# comparison
bf.pss.sec.2.sps <- bayes_factor(pss.sec.2.sps,pss.sec.2,log=T) # 330.23330
# testing
hypothesis(pss.sec.2.sps, 'secondary > 0') # Inf
hypothesis(pss.sec.2.sps, 'sps <0') # Inf
hypothesis(pss.sec.2.sps, 'secondary:sps=0') # 39.14. including 0
# only main effects

### 4
# RES ~ Sec * sps

res.sec.2 <- brms::brm(resilience ~ secondary+ gender + education + work_location + age+
                         SSS_faml+ relationship_status+
                         (1+secondary|residing_country),
                       data=data.filtered, family = gaussian(),
                       cores=4,chains=4, save_pars = save_pars(all = T),
                       sample_prior ='yes', seed=1660415,prior=prior.coef)
res.sec.2.sps <- brms::brm(resilience ~ secondary*sps+ gender + education + work_location + age+
                             SSS_faml+ relationship_status+
                             (1+secondary+sps|residing_country),
                           data=data.filtered, family = gaussian(),
                           cores=4,chains=4, save_pars = save_pars(all = T),
                           sample_prior ='yes', seed=1660415,prior=prior.coef)

# comparison
bf.res.sec.2.sps <- bayes_factor(res.sec.2.sps,res.sec.2,log=T) # 205.32347

# testing
hypothesis(res.sec.2.sps, 'secondary < 0') # Inf
hypothesis(res.sec.2.sps, 'sps >0') # Inf
hypothesis(res.sec.2.sps, 'secondary:sps=0') # .48. slight minus, but not sig.
1/.48 # 2.083333
