library(psych)
library(lavaan)
library(sirt)
library(MASS)



aligned.factor.scores <- function(lambda,nu,y){
  #calculate inverse matrix
  lambda1 <- ginv((lambda))
  #create matrix for nu
  ns <- nrow(y)
  nus <- matrix(nu,nrow=ns, ncol=length(nu), byrow=T)
  # y - nu
  y_nu <- y - nu
  # f = inv(lambda)*(y-nu)
  F <- lambda1 %*% t(as.matrix(y_nu))
}

#extract compliance_related variable names
# new model: without compliance_3 and 7

items.PSUP <- c('secondary_stressors__1','secondary_stressors__2',
                'secondary_stressors__3','secondary_stressors__4')
cfa.model.PSUP <-'RES =~ secondary_stressors__1 + secondary_stressors__2+secondary_stressors__3+
  secondary_stressors__4'


#Load the cleaned cvs file

#load data
data<-read.csv('../../Final_COVIDiSTRESS_Vol2_cleaned.csv')



#extract languages with n>=100
n.langs <- table(data$UserLanguage)
list.langs <-labels(n.langs)[[1]]
langs.include <-list.langs[n.langs>=100]
n.include <- n.langs[n.langs>=100]

#extract data
for (i in 1:length(langs.include)){
  if (i == 1){
    data.mi <-data[data$UserLanguage == langs.include[i],]
 }else{
   current <- data[data$UserLanguage == langs.include[i],]
   data.mi <- rbind(data.mi,current)
 }
}

#set and examine fitmeasures
fits <- c('rmsea.scaled','srmr','cfi.scaled','tli.scaled')


#####
#resilience

#general CFA: resilience
cfa.whole.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                      'UserLanguage')
fitMeasures(cfa.whole.sps)[fits]
# 0.1042919    0.0309894    0.9678559    0.9035678  mediocre

# then metric
cfa.metric.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                      'UserLanguage',group.equal ='loadings')
fitMeasures(cfa.metric.sps)[fits]
#   0.06707948   0.04305488   0.96746798   0.96010672 Accet

# then scalar
cfa.scalar.sps <-cfa(model=cfa.model.PSUP, data=data.mi,estimator='WLSMV', group=
                       'UserLanguage',group.equal =c('loadings','intercepts'))
fitMeasures(cfa.scalar.sps)[fits]
#    0.11497732   0.08032175   0.84791306   0.88279539 NG

# extract parameters
par.sps <- invariance_alignment_cfa_config(dat = data.mi[,items.PSUP], 
                                            group = data.mi$UserLanguage,
                                           estimator='WLSMV')
# do alignment
mod1.sps <- invariance.alignment(lambda = par.sps$lambda, nu =
                                   par.sps$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
# test performance
mod1.sps$es.invariance['R2',]
#  loadings intercepts 
#   0.9774234  0.9809317    
# all ≥ 75%. well addressed.

# monte carlo. IRT test (optional)
#--- find parameter constraints for prespecified tolerance
cmod1.sps <- sirt::invariance_alignment_constraints(model=mod1.sps, nu_parm_tol=.4,
                                                lambda_parm_tol=.4 )

par.cle <- par.sps

# extract params directly from CFA model

params <- parameterEstimates(cfa.whole.sps)
par.cle$lambda <-  matrix( params[ params$op=="=~", "est"], nrow=28,  byrow=T)
colnames(par.cle$lambda) <- items.PSUP
par.cle$nu <- matrix( params[ params$op=="~1"  & params$se !=0, "est" ], nrow=28,  byrow=T)
colnames(par.cle$nu) <- items.PSUP

# simulation function
simulation_CLE <- function(times,n,data,n.include,seed=1){
  # get data
  data.mi <- data
  items.PSUP <- c('secondary_stressors__1','secondary_stressors__2',
                  'secondary_stressors__3','secondary_stressors__4')

  compliance <- items.PSUP
  
  # list for return
  cor.mean <- rep(0,times)
  cor.var <- rep(0,times)
  R2.loading <- rep(0,times)
  R2.intercept <- rep(0,times)
  
  # simulation replication
  set.seed(seed)
  # repeat
  j <- 1
  while(1){

    G <- 28 # groups
    I <- 4 # items

    
    # lambda, nu, and error_var
    err_var.cle <- matrix(1, nrow=G,ncol=I)
    # simulate data
    # enter group mu and sigma
    data.mi$COMP <- rowMeans(data.mi[,compliance])
    mu<-scale(aggregate(x=data.mi$COMP,
              by = list(data.mi$UserLanguage),
              FUN=mean, na.rm=T)[,2])[,1]
    sigma <- (aggregate(x=data.mi$COMP,
                       by = list(data.mi$UserLanguage),
                       FUN=sd, na.rm=T)[,2])
    N <- rep(n,G)
    dat <- invariance_alignment_simulate(
      par.cle$nu,par.cle$lambda,err_var.cle,mu,sigma,N
    )
    #par.simul <- invariance_alignment_cfa_config(dat = dat[,items.PSUP], 
    #                                             group = dat$group,
    #                                             estimator='WLSMV'
    #                                             )
    
    # directly extract from cfa

    cfa.model.PSUP <-'PSUP =~ secondary_stressors__1 + secondary_stressors__2+
 secondary_stressors__3 +secondary_stressors__4'

    cfa.test <- cfa(cfa.model.PSUP,dat,group='group',estimator='WLSMV')
    params <- parameterEstimates(cfa.test)
    lambda <-  matrix( params[ params$op=="=~", "est"], nrow=28,  byrow=T)
    nu <- matrix( params[ params$op=="~1"  & params$se !=0, "est" ], nrow=28,  byrow=T)
    err_var <- matrix( params[ params$op=="~~"  & params$rhs=='PSUP', "est" ], nrow=28,  byrow=T)
    
    #flag <- sum(par.simul$err_var < 0)
    flag <- sum(err_var < 0)
    if (flag > 0){
      next
    }
    


    #mod1.simul <- invariance.alignment(lambda = par.simul$lambda, nu =
    #                                     par.simul$nu, align.scale = c(0.2, 0.4), align.pow = c(0.25, 0.25))
    
    tryCatch(
    expr={
      mod1.simul <- invariance.alignment(lambda = lambda, 
                                         nu =  nu, align.scale = c(0.2, 0.4), 
                                         align.pow = c(0.25, 0.25))}
    ,
    error = function(e){
      flag <- 1
    },
    warning = function(wan){
      flag <- 1
    })
   
    if (exists('mod1.simul')==FALSE){
      flag <- 1
    }
    if (flag > 0){
      next
    }
    
     # true vs aligned scores
    
   
    
    cfa.model.simul <- cfa.model.PSUP
    cfa.simul <- cfa(cfa.model.simul,dat,estimator='WLSMV',group='group',
                     group.equal=c('loadings','intercepts'),meanstructure=T)
    
    # get group mean
    params.simul <- parameterEstimates(cfa.simul)
    alpha.simul <- params.simul[(params.simul$op=='~1')&(params.simul$lhs=='PSUP'),'est']
 
    # group mean correlation (Muthen 2018)
    correlation <- corr.test(alpha.simul,mod1.simul$pars$alpha0,method='spearman')$r
    
    
    
    psi.simul <- params.simul[(params.simul$op=='~~')&(params.simul$lhs=='PSUP'),'est']
    correlation.psi <- corr.test(psi.simul,mod1.simul$pars$psi0,method='spearman')$r
    
    cor.mean[j] <- correlation
    cor.var[j] <- correlation.psi
    
    R2.loading[j] <- mod1.simul$es.invariance['R2',1]
    R2.intercept[j] <- mod1.simul$es.invariance['R2',2]
    
    # minus error detection
    if (correlation < 0 ){
      flag <- 1
      
    }
    if (correlation.psi < 0){
      flag <- 1
    }
    if (mod1.simul$es.invariance['R2',1] < 0){
      flag <- 1
    }
    if (mod1.simul$es.invariance['R2',2] < 0){
      flag <- 1
    }
    
    
    if (flag ==0){
      break
    }
  }
  # make matrix
  # with error flag

  to.return <-cbind(cor.mean,cor.var,R2.loading,R2.intercept)
  to.return <- data.matrix(to.return)
  
  message  (sprintf('%d/%d Done',j,times))
  
  return(to.return)
  
}



# use five cores

library(foreach)
library(parallel)
library(doParallel)

cores <- 5
times <- 500

cl <- parallel::makeCluster(cores,type='FORK')
doParallel::registerDoParallel(cl)

# n = 100

start_100 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,100,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_100<-Sys.time()
elapsed_100 <- end_100 - start_100
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_100 <- now[[1]]
  }else{
    simulate_100 <- rbind(simulate_100,now[[i]])
  }
}
# save n = 100
write.csv(data.frame(simulate_100),file='simulate_100.csv',row.names = F)


# n=200
# multicore processing for n=200
start_200 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,200,data.mi,28,i)
#  message(sprintf('%d',i))
}
end_200<-Sys.time()
elapsed_200 <- end_200 - start_200
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_200 <- now[[1]]
  }else{
    simulate_200 <- rbind(simulate_200,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_200),file='simulate_200.csv',row.names = F)

# n= 500
start_500 <-Sys.time()
now <- foreach (i = seq(1,times)) %dopar%{
  # get result
  simulation_CLE(1,500,data.mi,28,i)
  #  message(sprintf('%d',i))
}
end_500<-Sys.time()
elapsed_500 <- end_500 - start_500
# merge result
for (i in 1:times){
  if (i == 1){
    simulate_500 <- now[[1]]
  }else{
    simulate_500 <- rbind(simulate_500,now[[i]])
  }
}
# save n = 200
write.csv(data.frame(simulate_500),file='simulate_500.csv',row.names = F)
# stop cluster
parallel::stopCluster(cl)