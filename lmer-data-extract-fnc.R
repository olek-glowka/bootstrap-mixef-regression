lmer.data.extract = function(lmer.mod, name=deparse(substitute(lmer.mod))){
  
  #extract predictor names & create data frame, attach other cols to new data frame
  mod.data = data.frame(summary(lmer.mod)$coefficients[,1])
  names(mod.data) = "estimate"
  mod.data$std.error = as.numeric(summary(lmer.mod)$coefficients[,2]) #std errors
  mod.data$df = as.numeric(summary(lmer.mod)$coefficients[,3]) #degrees of freedom
  mod.data$t.val = as.numeric(summary(lmer.mod)$coefficients[,4]) #t-values
  mod.data$p.val = as.numeric(summary(lmer.mod)$coefficients[,5]) #p-values 
  
  #extract AIC, BIC, logLik, deviance df.resid
  mod.data$AIC = as.numeric(summary(lmer.mod)$AIC[1])
  mod.data$BIC = as.numeric(summary(lmer.mod)$AICtab[2][1])
  mod.data$logLik = as.numeric(summary(lmer.mod)$AICtab[3][1])
  mod.data$deviance = as.numeric(summary(lmer.mod)$AICtab[4][1])
  mod.data$df.resid = as.numeric(summary(lmer.mod)$AICtab[5][1])
  
  #add number of datapoints
  mod.data$N = as.numeric(summary(lmer.mod)$devcomp$dims[1])
  
  #add model name
  mod.data$model = name
  
  return(mod.data)
  
}

lmer.ranef.data.extract = function(lmer.mod, name=deparse(substitute(lmer.mod))){

  #extract random effect variance, standard error, correlations between slope and intercept
  mod.data.ranef = as.data.frame(VarCorr(lmer.mod))

  mod.data.ranef$n.subj = as.numeric(summary(lmer.mod)$ngrps[1]) #number of subjects
  mod.data.ranef$n.item = as.numeric(summary(lmer.mod)$ngrps[2]) #number of items

  #add number of datapoints
  mod.data.ranef$N = as.numeric(summary(lmer.mod)$devcomp$dims[1])
  
  #add model name
  mod.data.ranef$model = name
  
  return(mod.data.ranef)

}

lmer.optim.data.extract = function(lmer.mod, name=deparse(substitute(lmer.mod))){
  
  optim.data = as.data.frame(lmer.mod@optinfo$derivs$gradient)
  names(optim.data) = "gradient"
  optim.data$optimizer = lmer.mod@optinfo$optimizer
  optim.data$conv.opt = lmer.mod@optinfo$conv$opt
  optim.data$code  = lmer.mod@optinfo$conv$lme4$code
  optim.data$fval = lmer.mod@optinfo$feval
  optim.data$val = lmer.mod@optinfo$val
  optim.data$model = name
  
  return(optim.data)    
  
}

