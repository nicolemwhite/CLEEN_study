#99_functions.R


cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
g.theme =   theme_bw() + theme(legend.text = element_text(size=12),axis.text.y = element_text(size=12),axis.text.x = element_text(size=12),legend.position = 'top',legend.direction = 'horizontal',legend.title = element_blank(),panel.grid.minor = element_blank())

roundz = function(x, digits){
  dformat = paste('%.', digits, 'f', sep='')
  x = sprintf(dformat, round(x, digits))
  return(x)
}

# function to change 0/1 to no/yes
change_factor_lvls <- function(x){
  return(factor(x,levels=0:1,labels=c('No','Yes')))
}
#GLMM fitting and summary function
inv.logit<-function(x){exp(x)/(1+exp(x))}
fit_glmm_hai <- function(indat=dat_primary,y='total_hai_all',mod.form='~(1|ward)+data_collection_period+phase',run.boot=F,n.boot=NULL,glm.opt=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))){
  
  #create model formula
  ad.form<- paste0('cbind(',y,',n_surveys-',y,')',mod.form)
  mod.form <- formula(ad.form)
  
  #fit binomial glmm
  mod.fit <- glmer(mod.form,data=indat,family = binomial(link = "logit"),control=glm.opt)
  aic.mod <- roundz(AIC(mod.fit),0)
  
  #summaries
  ftab.mod = summary(mod.fit)$coefficients %>% as.data.frame() %>% rownames_to_column(var='Parameter') %>% select(Parameter,Estimate,`Std. Error`,`Pr(>|z|)`) %>% mutate_at('Estimate',~exp(.)) %>%
    rename('p-value'=`Pr(>|z|)`) %>% mutate_at('p-value',~ifelse(.<0.0001,'<0.0001',roundz(.,5)))
  glmer_ci = exp(confint.merMod(mod.fit,method="Wald")) %>% as.data.frame() %>% rownames_to_column(var="Parameter")
  ftab.mod = ftab.mod %>% left_join(glmer_ci,by="Parameter")
  ftab.mod = ftab.mod %>% 
    mutate('Estimate (95% CI)' = paste0(roundz(Estimate,2),' (',roundz(`2.5 %`,2),' to ',roundz(`97.5 %`,2),')'),
           '% change in odds (95% CI)' = paste0(roundz(100*(Estimate-1),1),' (',roundz(100*(`2.5 %`-1),1),' to ',roundz(100*(`97.5 %`-1),1),')')) %>%
    select(Parameter,`Estimate (95% CI)`,`% change in odds (95% CI)`,`p-value`) %>% 
    mutate_at("Parameter", ~ case_when(
      .=='(Intercept)' ~ 'Intercept',
      grepl('phase',.) ~ 'Intervention exposure')) %>% filter(!is.na(Parameter)) %>% rename('OR (95% CI)'=`Estimate (95% CI)`)
  
  
  rtab.mod = lme4::ranef(mod.fit) %>% data.frame() %>% rename('ward'=grp) %>% mutate(conf.low = condval - 1.96*condsd,conf.high=condval + 1.96*condsd)
  #marginal effects - run outside of function
  time_effect = ggpredict(mod.fit,terms="data_collection_period [all]") %>% data.frame()
  #phase_effect = ggemmeans(mod.fit,terms=c("phase [0,1]"))  %>% data.frame()

  #phase_effect = NULL
  mod.name <- paste(y)
  
  if (run.boot==T){
    phase_var <- rownames(summary(mod.fit)$coefficients)[grepl('phase',rownames(summary(mod.fit)$coefficients))]
    nd = data.frame(phase=0:1,data_collection_period='1') %>% rename_at('phase',~phase_var)
    boot_output = bootMer(mod.fit,FUN=function(.) c(inv.logit(predict(object=.,newdata=nd,re.form=~0)),summary(.)$coefficients[phase_var,c("Estimate","z value")]),nsim=n.boot)
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,boot_output=boot_output$t,time_effect=time_effect)
  }
  if (run.boot==F){
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,time_effect=time_effect)
    
  }
  
  return(out)
  
}

fit_glm_hai <- function(indat=dat_primary,y='total_hai_all',mod.form='~data_collection_period+phase',run.boot=F,n.boot=NULL){
  
  #create model formula
  ad.form<- paste0('cbind(',y,',n_surveys-',y,')',mod.form)
  mod.form <- formula(ad.form)
  
  #fit binomial glmm
  mod.fit <- glm(mod.form,data=indat,family = binomial(link = "logit"))
  aic.mod <- roundz(AIC(mod.fit),0)
  
  #summaries
  ftab.mod = summary(mod.fit)$coefficients %>% as.data.frame() %>% rownames_to_column(var='Parameter') %>% select(Parameter,Estimate,`Std. Error`,`Pr(>|z|)`) %>% mutate_at('Estimate',~exp(.)) %>%
    rename('p-value'=`Pr(>|z|)`) %>% mutate_at('p-value',~ifelse(.<0.001,'<0.001',roundz(.,3)))
  glmer_ci = exp(confint(mod.fit,method="Wald")) %>% as.data.frame() %>% rownames_to_column(var="Parameter")
  ftab.mod = ftab.mod %>% left_join(glmer_ci,by="Parameter")
  ftab.mod = ftab.mod %>% 
    mutate('Estimate (95% CI)' = paste0(roundz(Estimate,2),' (',roundz(`2.5 %`,2),' to ',roundz(`97.5 %`,2),')'),
           '% change in odds (95% CI)' = paste0(roundz(100*(Estimate-1),1),' (',roundz(100*(`2.5 %`-1),1),' to ',roundz(100*(`97.5 %`-1),1),')')) %>%
    select(Parameter,`Estimate (95% CI)`,`% change in odds (95% CI)`,`p-value`) %>% 
    mutate_at("Parameter", ~ case_when(
      .=='(Intercept)' ~ 'Intercept',
      grepl('phase',.) ~ 'Intervention exposure')) %>% filter(!is.na(Parameter)) %>% rename('OR (95% CI)'=`Estimate (95% CI)`)
  
  
  #marginal effects - run outside of function
  time_effect = ggpredict(mod.fit,terms="data_collection_period [all]") %>% data.frame()
  #phase_effect = ggemmeans(mod.fit,terms=c("phase [0,1]"))  %>% data.frame()
  
  #phase_effect = NULL
  mod.name <- paste(y)
  
  #if (run.boot==T){
  #  phase_var <- rownames(summary(mod.fit)$coefficients)[grepl('phase',rownames(summary(mod.fit)$coefficients))]
  #  nd = data.frame(phase=0:1,data_collection_period='1') %>% rename_at('phase',~phase_var)
  #  boot_output = bootMer(mod.fit,FUN=function(.) c(inv.logit(predict(object=.,newdata=nd,re.form=~0)),summary(.)$coefficients[phase_var,c("Estimate","z value")]),nsim=n.boot)
  #  out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,boot_output=boot_output$t,time_effect=time_effect)
  #}
  #if (run.boot==F){
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,time_effect=time_effect)
    
  #}
  
  return(out)
  
}

glmm_finv <- function(x,fam=NA) {case_when(fam=='poisson' ~ exp(x),fam=='binomial' ~ exp(x)/(1+exp(x)),TRUE ~ x)}

summarise_glmm <- function(mod.fit,est.label = 'OR (95% CI)'){
  glm.fam = ifelse(!is.null(summary(mod.fit)$family),summary(mod.fit)$family,NA)
  
  #summaries
  ftab.mod = summary(mod.fit)$coefficients %>% as.data.frame() %>% rownames_to_column(var='Parameter') %>% select(Parameter,Estimate,`Std. Error`,any_of(c('Pr(>|z|)','Pr(>|t|)'))) %>% mutate_at('Estimate',~glmm_finv(.,fam=glm.fam)) %>%
    rename('p-value'=starts_with('Pr')) %>% mutate_at('p-value',~ifelse(.<0.0001,'<0.0001',roundz(.,5)))
  glmer_ci = glmm_finv(confint.merMod(mod.fit,method="Wald"),fam=glm.fam) %>% as.data.frame() %>% rownames_to_column(var="Parameter")
  ftab.mod = ftab.mod %>% left_join(glmer_ci,by="Parameter")
  ftab.mod = ftab.mod %>% 
    mutate('Estimate (95% CI)' = paste0(roundz(Estimate,2),' (',roundz(V1,2),' to ',roundz(V2,2),')'),
           '% change (95% CI)' = paste0(roundz(100*(Estimate-1),1),' (',roundz(100*(V1-1),1),' to ',roundz(100*(V2-1),1),')')) %>%
    select(Parameter,`Estimate (95% CI)`,`% change (95% CI)`,`p-value`) %>% 
    mutate_at("Parameter", ~ case_when(
      .=='(Intercept)' ~ 'Intercept',
      grepl('phase',.) ~ 'Intervention exposure')) %>% filter(!is.na(Parameter)) %>% rename(!!est.label:=`Estimate (95% CI)`)
  
  return(ftab.mod)
}


fit_glmm_uv <- function(indat=dat_secondary,y='cbind(n_clean,n_audited - n_clean)',mod.form='~(1|ward)+data_collection_period+phase',run.boot=F,n.boot=NULL,glm.opt=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))){
  
  #create model formula
  ad.form<- paste0(y,mod.form)
  mod.form <- formula(ad.form)
  
  #fit binomial glmm
  mod.fit <- glmer(mod.form,data=indat,family = binomial(link = "logit"),control=glm.opt)
  aic.mod <- roundz(AIC(mod.fit),0)
  
  #summaries
  ftab.mod = summary(mod.fit)$coefficients %>% as.data.frame() %>% rownames_to_column(var='Parameter') %>% select(Parameter,Estimate,`Std. Error`,`Pr(>|z|)`) %>% mutate_at('Estimate',~exp(.)) %>%
    rename('p-value'=`Pr(>|z|)`) %>% mutate_at('p-value',~ifelse(.<0.0001,'<0.0001',roundz(.,5)))
  glmer_ci = exp(confint.merMod(mod.fit,method="Wald")) %>% as.data.frame() %>% rownames_to_column(var="Parameter")
  ftab.mod = ftab.mod %>% left_join(glmer_ci,by="Parameter")
  ftab.mod = ftab.mod %>% 
    mutate('Estimate (95% CI)' = paste0(roundz(Estimate,2),' (',roundz(`2.5 %`,2),' to ',roundz(`97.5 %`,2),')'),
           '% change in odds (95% CI)' = paste0(roundz(100*(Estimate-1),1),' (',roundz(100*(`2.5 %`-1),1),' to ',roundz(100*(`97.5 %`-1),1),')')) %>%
    select(Parameter,`Estimate (95% CI)`,`% change in odds (95% CI)`,`p-value`) %>% 
    mutate_at("Parameter", ~ case_when(
      .=='(Intercept)' ~ 'Intercept',
      .=='phase' ~ 'Intervention exposure',
      .=='time_since_intervention' ~ 'Time since intervention exposure (+2 weeks)')) %>% filter(!is.na(Parameter)) %>% rename('OR (95% CI)'=`Estimate (95% CI)`)
  
  
  rtab.mod = lme4::ranef(mod.fit) %>% data.frame() %>% rename('ward'=grp) %>% mutate(conf.low = condval - 1.96*condsd,conf.high=condval + 1.96*condsd)
  time_effect = ggpredict(mod.fit,'data_collection_period') %>% data.frame()
  
  #phase_covars = intersect(rownames(summary(mod.fit)$coefficients),c('phase','time_since_intervention'))
  #phase_effect = ggeffect(mod.fit,terms=paste(phase_covars, '[all]')) #may need to debug
  mod.name <- paste(y)
  
  if (run.boot==T){
    phase_var <- rownames(summary(mod.fit)$coefficients)[grepl('phase|time_since_intervention',rownames(summary(mod.fit)$coefficients))]
    #newdat = data.frame(time_label = c(seq(-12,12,2),-0.001),data_collection_period='1') %>% mutate(phase=ifelse(time_label<0,0,1),time_since_intervention = pmax(0,time_label))
    newdat = data.frame(time_label = sort(c(seq(-15,15,1),0)),data_collection_period='1',phase=rep(0:1,each=16)) %>% mutate(time_since_intervention=pmax(0,time_label))
    boot_output = bootMer(mod.fit,FUN=function(.) c(inv.logit(predict(object=.,newdata=newdat,re.form=~0)),summary(.)$coefficients[phase_var,c("Estimate","z value")]),nsim=n.boot)
    colnames(boot_output$t) <- c(paste('time_label = ',newdat$time_label),paste(phase_var,rep(c('Estimate','z'),each=length(phase_var)),sep='_'))
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,boot_output=boot_output$t,newdat=newdat,time_effect=time_effect)
  }
  if (run.boot==F){
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,time_effect=time_effect)
    
  }
  
  return(out)
  
}

#poisson log linear models by HAI type
fit_glmm_hai_log <- function(indat=dat_primary,y='total_hai_all',mod.form='~(1|ward)+data_collection_period+phase+offset(log(n_surveys))',run.boot=F,n.boot=NULL,glm.opt=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))){
  
  #create model formula
  ad.form<- paste0(y,mod.form)
  mod.form <- formula(ad.form)
  
  #fit binomial glmm
  mod.fit <- glmer(mod.form,data=indat,family = poisson(link = "log"),control=glm.opt)
  aic.mod <- roundz(AIC(mod.fit),0)
  
  #summaries
  ftab.mod = summary(mod.fit)$coefficients %>% as.data.frame() %>% rownames_to_column(var='Parameter') %>% select(Parameter,Estimate,`Std. Error`,`Pr(>|z|)`) %>% mutate_at('Estimate',~exp(.)) %>%
    rename('p-value'=`Pr(>|z|)`) %>% mutate_at('p-value',~ifelse(.<0.001,'<0.001',roundz(.,3)))
  glmer_ci = exp(confint.merMod(mod.fit,method="Wald")) %>% as.data.frame() %>% rownames_to_column(var="Parameter")
  ftab.mod = ftab.mod %>% left_join(glmer_ci,by="Parameter")
  ftab.mod = ftab.mod %>% 
    mutate('Estimate (95% CI)' = paste0(roundz(Estimate,2),' (',roundz(`2.5 %`,2),' to ',roundz(`97.5 %`,2),')'),
           '% change (95% CI)' = paste0(roundz(100*(Estimate-1),1),' (',roundz(100*(`2.5 %`-1),1),' to ',roundz(100*(`97.5 %`-1),1),')')) %>%
    select(Parameter,`Estimate (95% CI)`,`% change (95% CI)`,`p-value`) %>% 
    mutate_at("Parameter", ~ case_when(
      .=='(Intercept)' ~ 'Intercept',
      grepl('phase',.) ~ 'Intervention exposure')) %>% filter(!is.na(Parameter)) %>% rename('RR (95% CI)'=`Estimate (95% CI)`)
  
  
  rtab.mod = lme4::ranef(mod.fit) %>% data.frame() %>% rename('ward'=grp) %>% mutate(conf.low = condval - 1.96*condsd,conf.high=condval + 1.96*condsd)
  #marginal effects - run outside of function
  time_effect = ggpredict(mod.fit,terms="data_collection_period [all]") %>% data.frame()
  #phase_effect = ggemmeans(mod.fit,terms=c("phase [0,1]"))  %>% data.frame()
  
  #phase_effect = NULL
  mod.name <- paste(y)
  
  if (run.boot==T){
    phase_var <- rownames(summary(mod.fit)$coefficients)[grepl('phase',rownames(summary(mod.fit)$coefficients))]
    nd = data.frame(phase=0:1,data_collection_period='1') %>% rename_at('phase',~phase_var)
    boot_output = bootMer(mod.fit,FUN=function(.) c(exp(predict(object=.,newdata=nd,re.form=~0)),summary(.)$coefficients[phase_var,c("Estimate","z value")]),nsim=n.boot)
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,boot_output=boot_output$t,time_effect=time_effect)
  }
  if (run.boot==F){
    out = list(mod.name = mod.name,mod.fit=mod.fit,mod.form = mod.form,aic.mod=aic.mod,ftab.mod = ftab.mod,rtab.mod=rtab.mod,time_effect=time_effect)
    
  }
  
  return(out)
  
}


#bootstrapped p-values for glmm
#steps:
##i. Fit the null model to the data
##ii. Simulate R datasets from the null model
##iii. Fit null and alternative models to each simulated dataset
##iv. Store test stats
### Bootstrapped control and intervention estimates added
require(progress)
bootstrap_glmm_pvalue = function(indat=dat_primary,y='total_hai_all',n='n_surveys',mod.form='~(1|ward)+data_collection_period',phase_var='phase',R=500){
  pb <- progress_bar$new(format = "[:bar] :current/:total (:percent)", total = R)  #create model formulas for null and alt models
  ad.form<- paste0('cbind(',y,',',n,'-',y,')',mod.form)
  mod.form_0 <- formula(ad.form)
  mod.form_1 <- formula(paste(ad.form,paste(phase_var,collapse='+'),sep='+'))
  
  #fit binomial glmm
  mod_0 = glmer(mod.form_0,data=indat,family=binomial('logit'),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  mod_1 = glmer(mod.form_1,data=indat,family=binomial('logit'),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000)))
  
  #observed test statistics
  ll_obs = summary(mod_1)$AICtab['logLik'] - summary(mod_0)$AICtab['logLik']
  t_obs = summary(mod_1)$coefficients[phase_var,'z value']
  
  #random effect SD
  sd_0 = attr(summary(mod_0)$varcor$ward,'stddev')
  sd_1 = attr(summary(mod_1)$varcor$ward,'stddev')  
  
  #newdat = data.frame(time_label = c(seq(-12,12,2),-0.001),data_collection_period='1') %>% mutate(phase=ifelse(time_label<0,0,1),time_since_intervention = pmax(0,time_label))
  newdat = data.frame(time_label = sort(c(seq(-12,12,2),0)),data_collection_period='1',phase=rep(0:1,each=7)) %>% mutate(time_since_intervention=pmax(0,time_label))
  
  
  newdat = indat %>% select(ward,data_collection_period,all_of(phase_var),all_of(n)) %>% rename('n'=n)
  
  #simulate R datasets for the null model, fit null and alt models to each simulated dataset
  ll_r = rep(0,R)
  t_r = array(0,c(length(phase_var),R))
  est_r = NULL
  conv_r = NULL
  if(n=='n_surveys'){nd = data.frame(phase=0:1,data_collection_period='1') %>% rename_at('phase',~phase_var)}
  if(n=='n_audited'){nd = data.frame(time_label = sort(c(seq(-12,12,2),0)),data_collection_period='1',phase=rep(0:1,each=7)) %>% mutate(time_since_intervention=pmax(0,time_label))}
  K = length(unique(newdat$ward))
  for (r in 1:R){
    u_0 = rnorm(K,0,sd_0)
    u_1 = rnorm(K,0,sd_1)
    simdat = newdat %>% mutate(phat_0 = inv.logit(predict(mod_0,type='link',re.form=~0) + u_0[ward]),
                               phat_1 = inv.logit(predict(mod_1,type='link',re.form=~0) + u_1[ward])) %>% 
      rowwise() %>% mutate(yhat_0 = rbinom(1,n,phat_0),yhat_1 = rbinom(1,n,phat_1))
    mod_0_r = glmer(paste0('cbind(yhat_0,n-yhat_0)',mod.form),data=simdat,family=binomial('logit'),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) #null
    mod_1_r = glmer(paste0('cbind(yhat_0,n-yhat_0)',mod.form,'+',paste(phase_var,collapse='+')),data=simdat,family=binomial('logit'),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) #alt
    mod_2_r = glmer(paste0('cbind(yhat_1,n-yhat_1)',mod.form,'+',paste(phase_var,collapse='+')),data=simdat,family=binomial('logit'),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) #fixed effect(s)
    
    conv_r = rbind(conv_r,c(mod_0=ifelse(is.null(summary(mod_0_r)$optinfo$conv$lme4$code),0,summary(mod_0_r)$optinfo$conv$lme4$code),
                            mod_1=ifelse(is.null(summary(mod_1_r)$optinfo$conv$lme4$code),0,summary(mod_1_r)$optinfo$conv$lme4$code),
                            mod_2=ifelse(is.null(summary(mod_0_r)$optinfo$conv$lme4$code),0,summary(mod_0_r)$optinfo$conv$lme4$code))) #convergence codes
    
    
    ll_0_r = summary(mod_0_r)$AICtab['logLik']
    ll_1_r = summary(mod_1_r)$AICtab['logLik']
    ll_r[r] = ll_1_r - ll_0_r
    t_r[,r] = summary(mod_1_r)$coefficients[phase_var,'z value']
    est_r = rbind(est_r,c(inv.logit(predict(object=mod_2_r,newdata=nd,re.form=~0)),summary(mod_2_r)$coefficients[phase_var,"Estimate"]))
    
    pb$tick()
    Sys.sleep(1 / 100)
  }
  
  if(n=='n_surveys'){colnames(est_r) = c('Control','Intervention','Estimate')}
  if(n=='n_audited'){time_labs = paste(nd$phase,nd$time_label,sep=':'); colnames(est_r)<-c(time_labs,paste('Estimate',phase_var,sep=':'))}
  
  #p_boot_ll = (1+sum(ll_r>=ll_obs))/(R+1)
  #p_boot_t = (1+sum(abs(t_r)>=abs(t_obs)))/(R+1)
  conv_boot = apply(conv_r,1,function(x) max(abs(x)))
  
  return(list(t_obs = t_obs, ll_obs = ll_obs,est_boot = est_r,stat_boot = data.frame(ll=ll_r,zval=t_r),conv_boot=conv_boot))
  
}





#baseline characteristics functions


make_cat_tab_row = function(indata=dat,invar='gender',choose.level='Female',grouped=F,group_var=NULL,dps=1,label=NULL){
  if (choose.level=='all'){
    choose.level = indata %>% distinct(get(invar)) %>% pull()
    fdat = indata %>% filter(!is.na(get(invar))) %>% 
      group_by(level=get(invar)) %>% 
      summarise(n=n(),.groups='drop') %>% 
      mutate(denom=sum(n),perc=roundz(100*n/sum(n),dps)) %>% 
      mutate_if(is.numeric,function(x) format(x,big.mark=',',trim=T)) %>%
      mutate(stat = paste0(n,'/',denom,' (',perc,'%)')) %>% 
      filter(level %in% choose.level) %>% select(level,stat) %>% 
      add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1)  %>% add_column('grp'='All patients',.after='stat')
    
    if(grouped==T){
      ad = indata %>% filter(!is.na(get(invar)),!is.na(get(group_var))) %>% rename('grp'=(group_var)) %>%
        group_by(level=get(invar),grp) %>% 
        summarise(n=n(),.groups='drop') %>% group_by(grp) %>%
        mutate(denom=sum(n),perc=roundz(100*n/sum(n),dps)) %>% 
        mutate_if(is.numeric,function(x) format(x,big.mark=',',trim=T)) %>%
        mutate(stat = paste0(n,'/',denom,' (',perc,'%)')) %>% 
        filter(level %in% choose.level) %>% select(level,stat,grp) %>%
        add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1) 
      fdat = bind_rows(fdat,ad) %>% spread(grp,stat,fill='0 (0.0%)') %>% select(VAR,`All patients`,everything())
    }
  }
  else if (choose.level!='all'){
    fdat = indata %>% filter(!is.na(get(invar))) %>% 
      summarise(n=n(),r=sum(get(invar) %in% choose.level),perc=ifelse(n>0,roundz(100*r/n,dps),roundz(0,dps))) %>% 
      mutate_if(is.numeric,function(x) format(x,big.mark=',',trim=T)) %>%
      mutate(stat = paste0(r,'/',n,' (',perc,'%)')) %>% select(stat) %>% 
      add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1) %>% add_column('grp'='All patients',.after='stat')
    
    if(grouped==T){
      ad = indata %>% filter(!is.na(get(invar)),!is.na(get(group_var))) %>% rename('grp'=(group_var)) %>%
        group_by(grp) %>% 
        summarise(n=n(),r=sum(get(invar)==choose.level),perc=ifelse(n>0,roundz(100*r/n,dps),roundz(0,dps))) %>% 
        mutate_if(is.numeric,function(x) format(x,big.mark=',',trim=T)) %>%
        mutate(stat = paste0(r,'/',n,' (',perc,'%)')) %>% 
        select(stat,grp) %>% 
        add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1)  
      fdat = bind_rows(fdat,ad) %>% spread(grp,stat,fill='0 (0%)') %>% select(VAR,`All patients`,everything())
    }
    
    
  }
  
  return(fdat)
}

make_cont_tab_row = function(indata=dat,invar='age',grouped=F,group_var=NULL,dps=1,label=NULL){
  fdat = indata %>% filter(!is.na(get(invar))) %>%
    summarise(n=n(),med=median(get(invar)),q1=quantile(get(invar),.25),q3=quantile(get(invar),.75),avg=mean(get(invar)),stdev=sd(get(invar)),.groups='drop') %>% 
    mutate_if(is.numeric,function(x) format(round(x,dps),big.mark=',',trim=T)) %>%
    mutate(stat = paste0(med,' (',q1,' to ',q3,')\n',avg,' (',stdev,')')) %>% select(stat) %>% add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1) 
  if(grouped==T){
    ad = indata %>% filter(!is.na(get(invar)),!is.na(get(group_var))) %>% rename('grp'=(group_var)) %>%
      group_by(grp) %>%
      summarise(n=n(),med=median(get(invar)),q1=quantile(get(invar),.25),q3=quantile(get(invar),.75),avg=mean(get(invar)),stdev=sd(get(invar)),.groups='drop') %>% 
      mutate_if(is.numeric,function(x) format(round(x,dps),big.mark=',',trim=T)) %>%
      mutate(stat = paste0(med,' (',q1,' to ',q3,')\n',avg,' (',stdev,')')) %>% select(stat,grp) %>% 
      add_column(VAR=ifelse(!is.null(label),paste0(label),paste0(invar)),.before=1)  %>% spread(grp,stat,fill='--')
    fdat = full_join(fdat,ad,by='VAR') %>% rename('All patients'=stat)
  }
  
  return(fdat)
}

