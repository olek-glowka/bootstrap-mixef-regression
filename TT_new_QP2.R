#packages
require(lme4)
require(lmerTest)
require(languageR)
require(rms)
require(MASS)
require(plyr)
require(ggplot2)
require(grid)
require(gridExtra)
require(lattice)
require(effects)
require(MuMIn)
require(graphics)
require(glmnet)

#setup
#working directory
setwd("C:/Users/aglowka/Desktop/QP2/new_QP2/stats")

#ancillary functions (needs to be in same folder as working directory)
source("lmer-diag-fnc.R") #diagnostic plots for lmer mods
source("lmer-data-extract-fnc.R") #selected lmer results extractor

#load processed dset
TT.data = read.csv("TT.proc.data.csv", header=TRUE)

#############################################3

#iterative model fitting

# I. Holistic model
holistic.m.no.slopes = lmer(RT.log ~ FreqABCD.log.std +
                                freq.sub.PC1 +
                                freq.sub.PC2 +
                                freq.sub.PC3 +
                                freq.sub.PC4 +
                                PrevTrialPC1.log.std +
                                Length.log.std +
                                word.type.ratio.std +
                                PhraseABCD +
                                (1|Subject) +
                                (1|Item),
                              data = TT.data,
                              REML= FALSE)

holistic.m.subj.slope = lmer(RT.log ~ FreqABCD.log.std +
                              freq.sub.PC1 +
                              freq.sub.PC2 +
                              freq.sub.PC3 +
                              freq.sub.PC4 +
                              PrevTrialPC1.log.std +
                              Length.log.std +
                              word.type.ratio.std +
                              PhraseABCD +
                              (FreqABCD.log.std|Subject) +
                              (1|Item),
                            data = TT.data,
                            REML= FALSE)

anova(holistic.m.no.slopes, holistic.m.subj.slope)
#no improvement

holistic.m.item.slope = lmer(RT.log ~ FreqABCD.log.std +
                               freq.sub.PC1 +
                               freq.sub.PC2 +
                               freq.sub.PC3 +
                               freq.sub.PC4 +
                               PrevTrialPC1.log.std +
                               Length.log.std +
                               word.type.ratio.std +
                               PhraseABCD +
                               (1|Subject) +
                               (FreqABCD.log.std|Item),
                             data = TT.data,
                             REML= FALSE)

anova(holistic.m.no.slopes, holistic.m.item.slope)
#model failed to converge, keeping no slopes model

#compare best model to null model
holistic.m.null = lmer(RT.log ~ 
                              freq.sub.PC1 +
                              freq.sub.PC2 +
                              freq.sub.PC3 +
                              freq.sub.PC4 +
                              PrevTrialPC1.log.std +
                              Length.log.std +
                              word.type.ratio.std +
                              PhraseABCD +
                              (1|Subject) +
                              (1|Item),
                            data = TT.data,
                            REML=FALSE)

anova(holistic.m.null, holistic.m.no.slopes)
#significantly better than null at p < 0.001***

holistic.best.m = holistic.m.no.slopes

##########################

#model diagnostics

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(holistic.best.m), 2:9)$cnumber
#2.50 low collinearity

lmer.diag.mplot(holistic.best.m)

#Fixed effects
dotplot(fixef(holistic.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#Random effects
dotplot(ranef(holistic.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#######################################

#regression data extraction
write.csv(lmer.data.extract(holistic.best.m), "TT.holistic.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(holistic.best.m), "TT.holistic.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(holistic.best.m), "TT.holistic.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(holistic.best.m@optinfo$derivs$Hessian), "TT.holistic.m.hessian.orig.csv", row.names=TRUE)
write(unlist(holistic.best.m@optinfo$conv$lme4$messages), "TT.holistic.m.warnings.orig.txt")

#####################################

# II. Incremental model
incr.m.no.slopes = lmer(RT.log ~ LogitABCD.neg.log.std +
                                 MIABCD.neg.log.std +
                              PrevTrialPC1.log.std +
                              Length.log.std +
                              word.type.ratio.std +
                              PhraseABCD +
                              (1|Subject) +
                              (1|Item),
                            data = TT.data,
                            REML= FALSE)

incr.m.logit.subj.slope.MI.no.slopes = lmer(RT.log ~ LogitABCD.neg.log.std +
                                                     MIABCD.neg.log.std +
                                                  PrevTrialPC1.log.std +
                                                  Length.log.std +
                                                  word.type.ratio.std +
                                                  PhraseABCD +
                                                  (LogitABCD.neg.log.std|Subject) +
                                                  (1|Item),
                                                data = TT.data,
                                                REML= FALSE)

anova(incr.m.no.slopes, incr.m.logit.subj.slope.MI.no.slopes)
#no improvement

incr.m.logit.no.slopes.MI.subj.slope = lmer(RT.log ~ LogitABCD.neg.log.std +
                                              MIABCD.neg.log.std +
                                              PrevTrialPC1.log.std +
                                              Length.log.std +
                                              word.type.ratio.std +
                                              PhraseABCD +
                                              (MIABCD.neg.log.std|Subject) +
                                              (1|Item),
                                            data = TT.data,
                                            REML= FALSE)

anova(incr.m.no.slopes, incr.m.logit.no.slopes.MI.subj.slope)
#no improvement

incr.m.logit.item.slope.MI.no.slopes = lmer(RT.log ~ LogitABCD.neg.log.std +
                                             MIABCD.neg.log.std +
                                             PrevTrialPC1.log.std +
                                             Length.log.std +
                                             word.type.ratio.std +
                                             PhraseABCD +
                                             (1|Subject) +
                                             (LogitABCD.neg.log.std|Item),
                                           data = TT.data,
                                           REML= FALSE)
#model failed to converge

incr.m.logit.no.slopes.MI.item.slope = lmer(RT.log ~ LogitABCD.neg.log.std +
                                              MIABCD.neg.log.std +
                                              PrevTrialPC1.log.std +
                                              Length.log.std +
                                              word.type.ratio.std +
                                              PhraseABCD +
                                              (1|Subject) +
                                              (MIABCD.neg.log.std|Item),
                                            data = TT.data,
                                            REML= FALSE)

anova(incr.m.no.slopes, incr.m.logit.no.slopes.MI.item.slope)
# no improvement, keeping no slopes model

#developing stepwise null models to compare to full model

#starting with just logit no slopes
incr.m.logit.only.no.slopes = lmer(RT.log ~ LogitABCD.neg.log.std +
                                              PrevTrialPC1.log.std +
                                              Length.log.std +
                                              word.type.ratio.std +
                                              PhraseABCD +
                                              (1|Subject) +
                                              (1|Item),
                                            data = TT.data,
                                            REML= FALSE)

#logit only with subject slope
incr.m.logit.only.subj.slope = lmer(RT.log ~ LogitABCD.neg.log.std +
                                     PrevTrialPC1.log.std +
                                     Length.log.std +
                                     word.type.ratio.std +
                                     PhraseABCD +
                                     (LogitABCD.neg.log.std|Subject) +
                                     (1|Item),
                                   data = TT.data,
                                   REML= FALSE)

anova(incr.m.logit.only.no.slopes, incr.m.logit.only.subj.slope)
# bo improvement keeping no slopes model for logit only

#logit with item slope
incr.m.logit.only.item.slope = lmer(RT.log ~ LogitABCD.neg.log.std +
                                     PrevTrialPC1.log.std +
                                     Length.log.std +
                                     word.type.ratio.std +
                                     PhraseABCD +
                                     (1|Subject) +
                                     (LogitABCD.neg.log.std|Item),
                                   data = TT.data,
                                   REML= FALSE)

#model failed to converge, keeping no slopes logit model for comparison with null model

incr.m.null = lmer(RT.log ~
                     PrevTrialPC1.log.std +
                     Length.log.std +
                     word.type.ratio.std +
                     PhraseABCD +
                     (1|Subject) +
                     (1|Item),
                   data = TT.data,
                   REML= FALSE)

anova(incr.m.null, incr.m.logit.only.no.slopes)
#no significant improvement in fit, IMPORTANT: using null model for model with added MI

#trying MI only, selecting best MI model to compare to null
incr.m.MI.only.no.slopes = lmer(RT.log ~ MIABCD.neg.log.std +
                                  PrevTrialPC1.log.std +
                                  Length.log.std +
                                  word.type.ratio.std +
                                  PhraseABCD +
                                  (1|Subject) +
                                  (1|Item),
                                data = TT.data,
                                REML= FALSE)

# MI only with subject slope
incr.m.MI.only.subj.slope = lmer(RT.log ~ MIABCD.neg.log.std +
                                  PrevTrialPC1.log.std +
                                  Length.log.std +
                                  word.type.ratio.std +
                                  PhraseABCD +
                                  (MIABCD.neg.log.std|Subject) +
                                  (1|Item),
                                data = TT.data,
                                REML= FALSE)

anova(incr.m.MI.only.no.slopes, incr.m.MI.only.subj.slope)
#no improvement keeping no slopes model

# MI only with item slope
incr.m.MI.only.item.slope = lmer(RT.log ~ MIABCD.neg.log.std +
                                   PrevTrialPC1.log.std +
                                   Length.log.std +
                                   word.type.ratio.std +
                                   PhraseABCD +
                                   (1|Subject) +
                                   (MIABCD.neg.log.std|Item),
                                 data = TT.data,
                                 REML= FALSE)

anova(incr.m.MI.only.no.slopes, incr.m.MI.only.item.slope)
#no improvement, choosing no slopes model as intermediate model to compare with null

anova(incr.m.null, incr.m.MI.only.no.slopes)
#significant improvement in fit with p<0.001***

#comparing best intermediate model to no slopes incremental model
anova(incr.m.MI.only.no.slopes, incr.m.no.slopes)
#no improvement either way

#comparing no slopes incremental model to null model
anova(incr.m.null, incr.m.no.slopes)
#significant improvement with p<0.001***

############################################

incr.best.m = incr.m.no.slopes

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(incr.best.m), 2:7)$cnumber
#2.73 low collinearity

lmer.diag.mplot(incr.best.m)

#Fixed effects
dotplot(fixef(incr.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#Random effects
dotplot(ranef(incr.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

###############################################

#regression data extraction
write.csv(lmer.data.extract(incr.best.m), "TT.incr.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(incr.best.m), "TT.incr.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(incr.best.m), "TT.incr.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(incr.best.m@optinfo$derivs$Hessian), "TT.incr.m.hessian.orig.csv", row.names=TRUE)
write(unlist(incr.best.m@optinfo$conv$lme4$messages), "incr.m.warnings.orig.txt")

###############################################

# III. Dual model

dual.m.no.slopes = lmer(RT.log ~     FreqABCD.log.std +
                                     MIABCD.neg.log.std +
                                     LogitABCD.neg.log.std +
                                   freq.sub.PC1 +
                                   freq.sub.PC2 +
                                   freq.sub.PC3 +
                                   freq.sub.PC4 +
                                   PrevTrialPC1.log.std +
                                   Length.log.std +
                                   word.type.ratio.std +
                                   PhraseABCD +
                                   (1|Subject) +
                                   (1|Item),
                                 data = TT.data,
                                 REML= FALSE)


### freq random effects 

dual.m.freq.subj.slope = lmer(RT.log ~     FreqABCD.log.std +
                                           MIABCD.neg.log.std +
                                           LogitABCD.neg.log.std +
                                            freq.sub.PC1 +
                                            freq.sub.PC2 +
                                            freq.sub.PC3 +
                                            freq.sub.PC4 +
                                            PrevTrialPC1.log.std +
                                            Length.log.std +
                                            word.type.ratio.std +
                                            PhraseABCD +
                                            (FreqABCD.log.std|Subject) +
                                            (1|Item),
                                          data = TT.data,
                                          REML= FALSE)

anova(dual.m.no.slopes, dual.m.freq.subj.slope)
#no improvement

dual.m.freq.item.slope = lmer(RT.log ~     FreqABCD.log.std +
                                           MIABCD.neg.log.std +
                                           LogitABCD.neg.log.std +
                                            freq.sub.PC1 +
                                            freq.sub.PC2 +
                                            freq.sub.PC3 +
                                            freq.sub.PC4 +
                                            PrevTrialPC1.log.std +
                                            Length.log.std +
                                            word.type.ratio.std +
                                            PhraseABCD +
                                            (1|Subject) +
                                            (FreqABCD.log.std|Item),
                                          data = TT.data,
                                          REML= FALSE)

anova(dual.m.no.slopes, dual.m.freq.item.slope)

#### logit random effects

dual.m.logit.subj.slope = lmer(RT.log ~     FreqABCD.log.std +
                                MIABCD.neg.log.std +
                                LogitABCD.neg.log.std +
                                 freq.sub.PC1 +
                                 freq.sub.PC2 +
                                 freq.sub.PC3 +
                                 freq.sub.PC4 +
                                PrevTrialPC1.log.std +
                                Length.log.std +
                                word.type.ratio.std +
                                PhraseABCD +
                                (LogitABCD.neg.log.std|Subject) +
                                (1|Item),
                              data = TT.data,
                              REML= FALSE)

anova(dual.m.no.slopes, dual.m.logit.subj.slope)
#no improvment

dual.m.logit.item.slope = lmer(RT.log ~     FreqABCD.log.std +
                                 MIABCD.neg.log.std +
                                 LogitABCD.neg.log.std +
                                 freq.sub.PC1 +
                                 freq.sub.PC2 +
                                 freq.sub.PC3 +
                                 freq.sub.PC4 +
                                 PrevTrialPC1.log.std +
                                 Length.log.std +
                                 word.type.ratio.std +
                                 PhraseABCD +
                                 (1|Subject) +
                                 (LogitABCD.neg.log.std|Item),
                               data = TT.data,
                               REML= FALSE)

#model failed to converge, keeping null model

### Mutual Information random effects
dual.m.MI.subj.slope = lmer(RT.log ~     FreqABCD.log.std +
                                 MIABCD.neg.log.std +
                                 LogitABCD.neg.log.std +
                                  freq.sub.PC1 +
                                  freq.sub.PC2 +
                                  freq.sub.PC3 +
                                  freq.sub.PC4 +
                                 PrevTrialPC1.log.std +
                                 Length.log.std +
                                 word.type.ratio.std +
                                 PhraseABCD +
                                 (MIABCD.neg.log.std|Subject) +
                                 (1|Item),
                               data = TT.data,
                               REML= FALSE)

anova(dual.m.no.slopes, dual.m.MI.subj.slope)
#no improvement

dual.m.MI.item.slope = lmer(RT.log ~     FreqABCD.log.std +
                              MIABCD.neg.log.std +
                              LogitABCD.neg.log.std +
                              freq.sub.PC1 +
                              freq.sub.PC2 +
                              freq.sub.PC3 +
                              freq.sub.PC4 +
                              PrevTrialPC1.log.std +
                              Length.log.std +
                              word.type.ratio.std +
                              PhraseABCD +
                              (MIABCD.neg.log.std|Subject) +
                              (1|Item),
                            data = TT.data,
                            REML= FALSE)

anova(dual.m.no.slopes, dual.m.MI.item.slope)
#no improvement

#No random slopes

#rename best model
dual.best.m = dual.m.no.slopes

###########################

#dual model diagnostics

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(dual.best.m), 2:12)$cnumber
#7.82 medium collinearity

lmer.diag.mplot(dual.best.m)

#Fixed effects
dotplot(fixef(dual.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#Random effects
dotplot(ranef(dual.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

###########################

#regression data extraction
write.csv(lmer.data.extract(dual.best.m), "TT.dual.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(dual.best.m), "TT.dual.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(dual.best.m), "TT.dual.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(dual.best.m@optinfo$derivs$Hessian), "TT.dual.m.hessian.orig.csv", row.names=TRUE)
write(unlist(dual.best.m@optinfo$conv$lme4$messages), "TT.dual.m.warnings.orig.txt")

###########################

# *** Regression bagging ****
setwd("D:/Q2_bootstraps/main_analysis/tt_main")
source("lmer-data-extract-boot-fnc.R") #selected lmer results extractor

# divert messages stream to file, so you can log warnings and track non-converged models

options(warn=1)
wngs=file("warnings_log.txt",open="w+", blocking=TRUE)
sink(wngs,type="message")

for(iter in 7808:10000){
  
  setwd("D:/Q2_bootstraps/main_analysis/tt_main")
  
  data = read.csv(paste("tt_bootstrap_sample_", iter,".csv", sep=""), header=TRUE)
  
  data$LogitABCD.neg.log.std =  -data$LogitABCD.log.std
  data$MIABCD.neg.log.std = -data$MIABCD.log.std
  
  #log model number to file so that any potential warnings appear underneath
  message(paste("holistic model #", iter, sep=""))
  
  holistic.mod = lmer(RT.log ~ FreqABCD.log.std +
                                freq.sub.PC1 +
                                freq.sub.PC2 +
                                freq.sub.PC3 +
                                freq.sub.PC4 +
                                PrevTrialPC1.log.std +
                                Length.log.std +
                                word.type.ratio.std +
                                PhraseABCD +
                                (1|Subject) +
                                (1|Item),
                              data = data,
                              REML= FALSE)
  
  message(paste("probabilistic mod #", iter, sep=""))
  
  prob.mod = lmer(RT.log ~ LogitABCD.neg.log.std +
                           MIABCD.neg.log.std +
                            PrevTrialPC1.log.std +
                            Length.log.std +
                            word.type.ratio.std +
                            PhraseABCD +
                            (1|Subject) +
                            (1|Item),
                          data = data,
                          REML= FALSE)
  
  message(paste("dual mod #", iter, sep=""))
  
  dual.mod = lmer(RT.log ~ FreqABCD.log.std +
                           MIABCD.neg.log.std +
                           LogitABCD.neg.log.std +
                              freq.sub.PC1 +
                              freq.sub.PC2 +
                              freq.sub.PC3 +
                              freq.sub.PC4 +
                              PrevTrialPC1.log.std +
                              Length.log.std +
                              word.type.ratio.std +
                              PhraseABCD +
                              (1|Subject) +
                              (1|Item),
                            data = data,
                            REML= FALSE)
  
  #write model results to file
  setwd("D:/Q2_bootstraps/main_analysis/tt_main/reg_boot_holistic")
  write.csv(lmer.data.extract.boot(holistic.mod, iter), paste("holistic.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(holistic.mod, iter), paste("holistic.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(holistic.mod, iter), paste("holistic.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(holistic.mod@optinfo$derivs$Hessian), paste("holistic.m.hessian_", iter, ".csv", sep=""))
  write(unlist(holistic.mod@optinfo$conv$lme4$messages), paste("holistic.m.warnings_", iter, ".txt", sep=""))
  
  setwd("D:/Q2_bootstraps/main_analysis/tt_main/reg_boot_prob")
  write.csv(lmer.data.extract.boot(prob.mod, iter), paste("prob.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(prob.mod, iter), paste("prob.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(prob.mod, iter), paste("prob.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(prob.mod@optinfo$derivs$Hessian), paste("prob.m.hessian_", iter, ".csv", sep=""))
  write(unlist(prob.mod@optinfo$conv$lme4$messages), paste("prob.m.warnings_", iter, ".txt", sep=""))
  
  setwd("D:/Q2_bootstraps/main_analysis/tt_main/reg_boot_dual")
  write.csv(lmer.data.extract.boot(dual.mod, iter), paste("dual.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(dual.mod, iter), paste("dual.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(dual.mod, iter), paste("dual.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(dual.mod@optinfo$derivs$Hessian), paste("dual.m.hessian_", iter, ".csv", sep=""))
  write(unlist(dual.mod@optinfo$conv$lme4$messages), paste("dual.m.warnings_", iter, ".txt", sep=""))
  
  
  cat("Bagging iteration", iter, "completed!\n")
  
}

#close log file & restore warnings stream to console
closeAllConnections()
