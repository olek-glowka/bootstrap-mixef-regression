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

#read in data from experiments 1 & 2
AS.data.raw = read.csv("as-meta.csv", header=TRUE)

#exclude error trials
AS.data.corr = AS.data.raw[AS.data.raw$corr==1,]

# How much data was lost after removing error trials?
paste("Out of ", nrow(AS.data.raw) ," points, ", nrow(AS.data.raw)-nrow(AS.data.corr) ," points were excluded after removing error trials, resulting in the loss of ", round((1 - (nrow(AS.data.corr)/nrow(AS.data.raw)))*100, digits=1), "% of the data.", sep = "")

#########################################################

#data transformation and further preprocessing

#make block order a factor
AS.data.corr$block = as.factor(AS.data.corr$block)

#convert RTs from raw to log and scale for trimming
AS.data.corr$RT.log = log(AS.data.corr$time)
AS.data.corr$RT.log.std = scale(AS.data.corr$RT.log)

#exclude outlier RTs and create a new trimmed dataset
AS.data = AS.data.corr[abs(AS.data.corr$RT.log.std) <= 2,]

# How much data was lost after trimming by z-scored RTs?
paste("Out of ", nrow(AS.data.corr) ," points, ", nrow(AS.data.corr)-nrow(AS.data) ," points were removed after z-scoring the depdendent variable (RT.log) by condition and excluding points 2SDs above the mean, resulting in the loss of ", round((1 - (nrow(AS.data)/nrow(AS.data.corr)))*100, digits=1), "% of the data from correct trials only. This constituents a further loss of ", round(((nrow(AS.data.corr) - nrow(AS.data))/nrow(AS.data.raw))*100, digits=1), "% of the entire data.", sep = "")

#log frequency counts
AS.data$FreqA.log = log(AS.data$FreqA)
AS.data$FreqB.log = log(AS.data$FreqB)
AS.data$FreqC.log = log(AS.data$FreqC)
AS.data$FreqD.log = log(AS.data$FreqD)
AS.data$FreqAB.log = log(AS.data$FreqAB)
AS.data$FreqBC.log = log(AS.data$FreqBC)
AS.data$FreqCD.log = log(AS.data$FreqCD)
AS.data$FreqABC.log = log(AS.data$FreqABC)
AS.data$FreqBCD.log = log(AS.data$FreqBCD)
AS.data$FreqABCD.log = log(AS.data$FreqABCD)

#compute logits (a proxy for conditional probabilites)
AS.data$LogitABCD.log = AS.data$FreqABCD.log / (AS.data$FreqABC.log - AS.data$FreqABCD.log)

#compute Mutual Information
AS.data$MIABCD.log = AS.data$FreqABCD.log / (AS.data$FreqA.log * AS.data$FreqB.log * AS.data$FreqC.log * AS.data$FreqD.log)

#create unscaled neg conditional probability and neg mutual information for correlation matrix
AS.data$LogitABCD.neg.log =  -AS.data$LogitABCD.log
AS.data$MIABCD.neg.log = -AS.data$MIABCD.log

#scale continuous variables
#scale number of letters
AS.data$nletter.std = scale(AS.data$nletter)

#scale frequency variables
AS.data$FreqA.log.std = scale(AS.data$FreqA.log)
AS.data$FreqB.log.std = scale(AS.data$FreqB.log)
AS.data$FreqC.log.std = scale(AS.data$FreqC.log)
AS.data$FreqD.log.std = scale(AS.data$FreqD.log)
AS.data$FreqAB.log.std = scale(AS.data$FreqAB.log)
AS.data$FreqBC.log.std = scale(AS.data$FreqBC.log)
AS.data$FreqCD.log.std = scale(AS.data$FreqCD.log)
AS.data$FreqABC.log.std = scale(AS.data$FreqABC.log)
AS.data$FreqBCD.log.std = scale(AS.data$FreqBCD.log)
AS.data$FreqABCD.log.std = scale(AS.data$FreqABCD.log)

#scale Logit
AS.data$LogitABCD.log.std = scale(AS.data$LogitABCD.log)

#scale mutual Information
AS.data$MIABCD.log.std = scale(AS.data$MIABCD.log)

#create negative conditional probability and negative mutual information
AS.data$LogitABCD.neg.log.std =  -AS.data$LogitABCD.log.std
AS.data$MIABCD.neg.log.std = -AS.data$MIABCD.log.std

#save pre-processed data
write.csv(AS.data, "AS.data.preproc.csv")

#########################################################
#DV correlation matrix
#Correlation tests to examine degree of independence of logged but unstandardized critical predictors

#freq vs. logit
cor.test(AS.data$FreqABCD.log, AS.data$LogitABCD.neg.log, method="kendall")

#freq vs. MI
cor.test(AS.data$FreqABCD.log, AS.data$MIABCD.neg.log, method="kendall")

#logit vs. MI
cor.test(AS.data$LogitABCD.neg.log, AS.data$MIABCD.neg.log, method="kendall")

#create a small dataset with just the tree critical variables to plot pairwise correlations

pairs.data = data.frame(AS.data$FreqABCD.log)
names(pairs.data) = "FreqABCD.log"
pairs.data$LogitABCD.neg.log = AS.data$LogitABCD.neg.log
pairs.data$MIABCD.neg.log = AS.data$MIABCD.neg.log

panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- as.numeric(cor.test(x, y, method= "kendall")[4][1])
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 2 #0.3/strwidth(txt) 
  
  test <- cor.test(x,y, method="kendall") 
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(.5, .5, txt, cex = cex) 
  text(.7, .55, Signif, cex=cex, col=1) 
}

pairs(pairs.data, lower.panel=panel.smooth, upper.panel=panel.cor)

#########################################################

#Perform PCA on Frequency subchunk controls to reduce the number of predictors in frequency models
freq.pca.dat = data.frame(AS.data$FreqA.log.std)
names(freq.pca.dat) = "FreqA.log.std"
freq.pca.dat$FreqB.log.std = AS.data$FreqB.log.std
freq.pca.dat$FreqC.log.std = AS.data$FreqC.log.std
freq.pca.dat$FreqD.log.std = AS.data$FreqD.log.std
freq.pca.dat$FreqAB.log.std = AS.data$FreqAB.log.std
freq.pca.dat$FreqBC.log.std = AS.data$FreqBC.log.std
freq.pca.dat$FreqCD.log.std = AS.data$FreqCD.log.std
freq.pca.dat$FreqABC.log.std = AS.data$FreqABC.log.std
freq.pca.dat$FreqBCD.log.std = AS.data$FreqBCD.log.std

#perform PCA
pr.out = prcomp(freq.pca.dat, scale=FALSE) #no need to scale again

#the variance explained by each PC component is obtained by squaring its standard deviation
pr.var = pr.out$sdev^2

#to compute the proportion of variance explained by each PC, we didive the variance explained by any given PC by the total variance explained by all PCs.
pve = pr.var/sum(pr.var)

#plot pve and cumulative pve
plot(pve, xlab = "Principal Compnent", ylab="Proportion of Variance Explained", ylim=c(0,1), type='b')
plot(cumsum(pve), xlab = "Principal Compnent", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1), type='b')

#alternative visualization
par(mar = rep(2, 4))
plot(pr.out)

sum(pve[1:4])

#add the scores for the first four Principal components to the trimmed data frame
AS.data$freq.sub.PC1 = pr.out$x[,1]
AS.data$freq.sub.PC2 = pr.out$x[,2]
AS.data$freq.sub.PC3 = pr.out$x[,3]
AS.data$freq.sub.PC4 = pr.out$x[,4]

###########

# I. Holistic Model

#no freq slopes model
holistic.no.slopes.m = lmer(RT.log ~ FreqABCD.log.std +
                                     AS.data$freq.sub.PC1 +
                                     AS.data$freq.sub.PC2 +
                                     AS.data$freq.sub.PC3 +
                                     AS.data$freq.sub.PC4 +
                                     block + 
                                     nletter.std +
                                     (1|subj) +
                                     (1|item),
                                 data = AS.data,
                                 REML = FALSE)
#freq subj slope
holistic.subj.m = lmer(RT.log ~ FreqABCD.log.std +
                                   AS.data$freq.sub.PC1 +
                                   AS.data$freq.sub.PC2 +
                                   AS.data$freq.sub.PC3 +
                                   AS.data$freq.sub.PC4 +
                                   block + 
                                   nletter.std +
                                   (FreqABCD.log.std|subj) +
                                   (1|item),
                                 data = AS.data,
                                 REML = FALSE)

anova(holistic.no.slopes.m, holistic.subj.m)
#no improvement in fit, keeping no freq slopes model

#freq item slope
holistic.item.m = lmer(RT.log ~ FreqABCD.log.std +
                              AS.data$freq.sub.PC1 +
                              AS.data$freq.sub.PC2 +
                              AS.data$freq.sub.PC3 +
                              AS.data$freq.sub.PC4 +
                              block + 
                              nletter.std +
                              (1|subj) +
                              (FreqABCD.log.std|item),
                            data = AS.data,
                            REML = FALSE)

anova(holistic.no.slopes.m, holistic.item.m)
#improvement in fit, adopting freq item mod

#freq subj & item slope
holistic.subj.item.m = lmer(RT.log ~ FreqABCD.log.std +
                              AS.data$freq.sub.PC1 +
                              AS.data$freq.sub.PC2 +
                              AS.data$freq.sub.PC3 +
                              AS.data$freq.sub.PC4 +
                              block + 
                              nletter.std +
                              (FreqABCD.log.std|subj) +
                              (FreqABCD.log.std|item),
                            data = AS.data,
                            REML = FALSE)

anova(holistic.item.m, holistic.subj.item.m)
#no improvement in fit improvement, keeping item slope model as best model

#null model without critical predictor for comparison
#random slopes not included in null model
holistic.null.m = lmer(RT.log ~ AS.data$freq.sub.PC1 +
                                AS.data$freq.sub.PC2 +
                                AS.data$freq.sub.PC3 +
                                AS.data$freq.sub.PC4 +
                                block + 
                                nletter.std +
                                (1|subj) +
                                (1|item),
                                data = AS.data,
                                REML = FALSE)
#better than null model?
anova(holistic.null.m, holistic.item.m)
#better than null!

holistic.best.m = holistic.item.m

##########################

#model diagnostics

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(holistic.best.m), 2:8)$cnumber
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
write.csv(lmer.data.extract(holistic.best.m), "AS.holistic.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(holistic.best.m), "AS.holistic.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(holistic.best.m), "AS.holistic.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(holistic.best.m@optinfo$derivs$Hessian), "AS.holistic.m.hessian.orig.csv", row.names=TRUE)
write(unlist(holistic.best.m@optinfo$conv$lme4$messages), "AS.holistic.m.warnings.orig.txt")

#####################################

# II. Incremental Model

#no slopes model
incr.logit.no.slopes.MI.no.slopes.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                                    MIABCD.neg.log.std + 
                                                    block + 
                                                    nletter.std +
                                                    (1|subj) +
                                                    (1|item),
                                                data = AS.data,
                                                REML = FALSE)
#adding logit subj slope
incr.logit.subj.MI.no.slopes.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                              MIABCD.neg.log.std + 
                                              block + 
                                              nletter.std +
                                              (LogitABCD.neg.log.std|subj) +
                                              (1|item),
                                        data = AS.data,
                                        REML = FALSE)

anova(incr.logit.no.slopes.MI.no.slopes.m, incr.logit.subj.MI.no.slopes.m)
#no improvement in fit, keeping no slopes model

#adding logit item slope
incr.logit.item.MI.no.slopes.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                        MIABCD.neg.log.std + 
                                        block + 
                                        nletter.std +
                                        (1|subj) +
                                        (LogitABCD.neg.log.std|item),
                                      data = AS.data,
                                      REML = FALSE) 

anova(incr.logit.no.slopes.MI.no.slopes.m, incr.logit.item.MI.no.slopes.m)
#improved fit, adopting logit item slope mod

#adding logit subj to item slope
incr.logit.subj.item.MI.no.slopes.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                        MIABCD.neg.log.std + 
                                        block + 
                                        nletter.std +
                                        (LogitABCD.neg.log.std|subj) +
                                        (LogitABCD.neg.log.std|item),
                                      data = AS.data,
                                      REML = FALSE)

anova(incr.logit.item.MI.no.slopes.m, incr.logit.subj.item.MI.no.slopes.m)
# no improvement, keeping logit item slope mod

#adding MI subj slope
incr.logit.item.MI.subj.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                        MIABCD.neg.log.std + 
                                        block + 
                                        nletter.std +
                                        (MIABCD.neg.log.std|subj) +
                                        (LogitABCD.neg.log.std|item),
                                      data = AS.data,
                                      REML = FALSE)

anova(incr.logit.item.MI.no.slopes.m, incr.logit.item.MI.subj.m)
# no improvement, keeping logit item MI no slopes mod

#adding MI item slope
incr.logit.item.MI.item.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                   MIABCD.neg.log.std + 
                                   block + 
                                   nletter.std +
                                   (1|subj) +
                                   (LogitABCD.neg.log.std|item) +
                                   (MIABCD.neg.log.std|item),
                                 data = AS.data,
                                 REML = FALSE)

anova(incr.logit.item.MI.no.slopes.m, incr.logit.item.MI.item.m)
#no improvement, keeping logit item MI no slopes

#adding MI item and removing logit item slope
incr.logit.no.slopes.MI.item.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                                   MIABCD.neg.log.std + 
                                   block + 
                                   nletter.std +
                                   (1|subj) +
                                   (MIABCD.neg.log.std|item),
                                 data = AS.data,
                                 REML = FALSE)

#If both logit item and MI item are better than no slopes, which model should I choose? 
# The model with larger decrease in AIC.

#how much better is MI item better than no slopes?
anova(incr.logit.no.slopes.MI.no.slopes.m, incr.logit.no.slopes.MI.item.m)

#how much better is logit item better than no slopes?
anova(incr.logit.no.slopes.MI.no.slopes.m, incr.logit.item.MI.no.slopes.m)

#MI item mod showed a larger decrease in AIC relative to the logit item model
#so adopting MI item model as best model

########################

#null model comparison

#null mod
incr.null.m = lmer(RT.log ~ block + 
                            nletter.std +
                            (1|subj) +
                            (1|item),
                            data = AS.data,
                            REML = FALSE)

#Best to compare null model with no slopes or item slope or subject slope or both?

#logit only with no slope for comparison
incr.logit.only.no.slope.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                           block + 
                           nletter.std +
                           (1|subj) +
                           (1|item),
                         data = AS.data,
                         REML = FALSE)

#logit only with subj slope
incr.logit.only.subj.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                           block + 
                           nletter.std +
                           (LogitABCD.neg.log.std|subj) +
                           (1|item),
                         data = AS.data,
                         REML = FALSE)

anova(incr.logit.only.no.slope.m, incr.logit.only.subj.m)
#no improvement compared to logit with no slopes

#logit with item slope
incr.logit.only.item.m = lmer(RT.log ~ LogitABCD.neg.log.std +
                           block + 
                           nletter.std +
                           (1|subj) +
                           (LogitABCD.neg.log.std|item),
                         data = AS.data,
                         REML = FALSE)

anova(incr.logit.only.no.slope.m, incr.logit.only.item.m)
#improvement compared to logit with no slopes, so using this model compare to null

anova(incr.null.m, incr.logit.only.item.m)
#significant improvement in fit, IMPORTANT: logit item slope becomes the null model for model with added MI

#now add MI
#Is it better to add MI with no slopes or MI with subj and/or item slope?
#That model selection has been done above, we saw that model with no logit slopes and MI item slope is best.
#So use that best random effect structure model to compare to null.
#The logic here is to compare the most optimal model with just logit 
#to most optimal model when both logit and MI are present, this is case the logit no slopes, MI item model 
#This justifies the apparent random effect "leaps" and allows us to obvite the complexity explosion by
#performing optimal random effect and fixed effect model selection separately.

anova(incr.logit.only.item.m,incr.logit.no.slopes.MI.item.m)
# significant improvement in fit!

# rename best model for ease of comparison
incr.best.m = incr.logit.no.slopes.MI.item.m 

###########################

#incremental model diagnostics

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(incr.best.m), 2:5)$cnumber
#2.73 low collinearity

lmer.diag.mplot(incr.best.m)
#majorly impacted by leverage, consider removing 8 points with leverage < 0.06

#Fixed effects
dotplot(fixef(incr.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#Random effects
dotplot(ranef(incr.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#########################

#regression data extraction
write.csv(lmer.data.extract(incr.best.m), "AS.incr.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(incr.best.m), "AS.incr.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(incr.best.m), "AS.incr.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(incr.best.m@optinfo$derivs$Hessian), "AS.incr.m.hessian.orig.csv", row.names=TRUE)
write(unlist(incr.best.m@optinfo$conv$lme4$messages), "AS.incr.m.warnings.orig.txt")

#########################

#If you need to refit the model after removing 7 high leverage points, use his code
AS.data$incr.low.lev = TRUE

#for (i in 1:nrow(AS.data)){
#  
#  if (as.numeric(hatvalues(incr.best.m)[i]) > 0.07){
#    
#    AS.data$incr.low.lev[i] = FALSE
#    
#  } 
#}

#sum(AS.data$incr.low.lev == FALSE) # 8 points with high leverage removed

######################################

# III. Dual-stream model

# There are three critical variables (X1, X2, X3). Each of them could potentially have a random slope.
# How to choose the best model without sucking out critical variable variance?
# Procedure:
# Step 1. Find best random effect structure independently for X1, X2 and X3

#*** freq slopes ***

#no freq slopes model
dual.no.slopes.m = lmer(RT.log ~ FreqABCD.log.std +
                              LogitABCD.neg.log.std +
                              MIABCD.neg.log.std +
                              AS.data$freq.sub.PC1 +
                              AS.data$freq.sub.PC2 +
                              AS.data$freq.sub.PC3 +
                              AS.data$freq.sub.PC4 +
                              block + 
                              nletter.std +
                              (1|subj) +
                              (1|item),
                            data = AS.data,
                            REML = FALSE)

#freq subj slope model
dual.freq.subj.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                          LogitABCD.neg.log.std +
                          MIABCD.neg.log.std +
                          AS.data$freq.sub.PC1 +
                          AS.data$freq.sub.PC2 +
                          AS.data$freq.sub.PC3 +
                          AS.data$freq.sub.PC4 +
                          block + 
                          nletter.std +
                          (FreqABCD.log.std|subj) +
                          (1|item),
                        data = AS.data,
                        REML = FALSE)

anova(dual.no.slopes.m, dual.freq.subj.slope.m)
#no improvement compared to no slopes mod

#freq item slope model
dual.freq.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                LogitABCD.neg.log.std +
                                MIABCD.neg.log.std +
                                AS.data$freq.sub.PC1 +
                                AS.data$freq.sub.PC2 +
                                AS.data$freq.sub.PC3 +
                                AS.data$freq.sub.PC4 +
                                block + 
                                nletter.std +
                                (1|subj) +
                                (FreqABCD.log.std|item),
                              data = AS.data,
                              REML = FALSE)

anova(dual.no.slopes.m, dual.freq.item.slope.m)
#improvement, adopting dual.freq.item.slope.m

#freq subj slope model
dual.freq.subj.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                LogitABCD.neg.log.std +
                                MIABCD.neg.log.std +
                                AS.data$freq.sub.PC1 +
                                AS.data$freq.sub.PC2 +
                                AS.data$freq.sub.PC3 +
                                AS.data$freq.sub.PC4 +
                                block + 
                                nletter.std +
                                (FreqABCD.log.std|subj) +
                                (FreqABCD.log.std|item),
                              data = AS.data,
                              REML = FALSE)

anova(dual.freq.item.slope.m, dual.freq.subj.item.slope.m)
#no improvement, keeping dual.freq.item.slope.m

#*** logit slopes ***


dual.logit.subj.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                          LogitABCD.neg.log.std +
                          MIABCD.neg.log.std +
                          AS.data$freq.sub.PC1 +
                          AS.data$freq.sub.PC2 +
                          AS.data$freq.sub.PC3 +
                          AS.data$freq.sub.PC4 +
                          block + 
                          nletter.std +
                          (LogitABCD.neg.log.std|subj) +
                          (1|item),
                        data = AS.data,
                        REML = FALSE)

anova(dual.no.slopes.m, dual.logit.subj.slope.m)
#no improvwment, keeping null mod

dual.logit.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                 LogitABCD.neg.log.std +
                                 MIABCD.neg.log.std +
                                 AS.data$freq.sub.PC1 +
                                 AS.data$freq.sub.PC2 +
                                 AS.data$freq.sub.PC3 +
                                 AS.data$freq.sub.PC4 +
                                 block + 
                                 nletter.std +
                                 (1|subj) +
                                 (LogitABCD.neg.log.std|item),
                               data = AS.data,
                               REML = FALSE)

anova(dual.no.slopes.m, dual.logit.item.slope.m)
#significant improvement, adopting logit item slope as best

#re-adding subject slope
dual.logit.subj.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                 LogitABCD.neg.log.std +
                                 MIABCD.neg.log.std +
                                 AS.data$freq.sub.PC1 +
                                 AS.data$freq.sub.PC2 +
                                 AS.data$freq.sub.PC3 +
                                 AS.data$freq.sub.PC4 +
                                 block + 
                                 nletter.std +
                                 (LogitABCD.neg.log.std|subj) +
                                 (LogitABCD.neg.log.std|item),
                               data = AS.data,
                               REML = FALSE)

anova(dual.logit.item.slope.m, dual.logit.subj.item.slope.m)
#no improvement

#*** MI slopes ***

dual.MI.subj.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                 LogitABCD.neg.log.std +
                                 MIABCD.neg.log.std +
                                 AS.data$freq.sub.PC1 +
                                 AS.data$freq.sub.PC2 +
                                 AS.data$freq.sub.PC3 +
                                 AS.data$freq.sub.PC4 +
                                 block + 
                                 nletter.std +
                                 (MIABCD.neg.log.std|subj) +
                                 (1|item),
                               data = AS.data,
                               REML = FALSE)

anova(dual.no.slopes.m, dual.MI.subj.slope.m)
#no improvwment, keeping null mod

dual.MI.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                 LogitABCD.neg.log.std +
                                 MIABCD.neg.log.std +
                                 AS.data$freq.sub.PC1 +
                                 AS.data$freq.sub.PC2 +
                                 AS.data$freq.sub.PC3 +
                                 AS.data$freq.sub.PC4 +
                                 block + 
                                 nletter.std +
                                 (1|subj) +
                                 (MIABCD.neg.log.std|item),
                               data = AS.data,
                               REML = FALSE)

anova(dual.no.slopes.m, dual.MI.item.slope.m)
#significant improvement, adopting MI item slope as best

#re-adding subject slope
dual.MI.subj.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                      LogitABCD.neg.log.std +
                                      MIABCD.neg.log.std +
                                      AS.data$freq.sub.PC1 +
                                      AS.data$freq.sub.PC2 +
                                      AS.data$freq.sub.PC3 +
                                      AS.data$freq.sub.PC4 +
                                      block + 
                                      nletter.std +
                                      (MIABCD.neg.log.std|subj) +
                                      (MIABCD.neg.log.std|item),
                                    data = AS.data,
                                    REML = FALSE)

anova(dual.MI.item.slope.m, dual.MI.subj.item.slope.m)
#no improvement

#for all three models, item slope only was best
#we choose the item slope for the variable whose model has the lowest AIC

as.numeric(summary(dual.freq.item.slope.m)$AIC[1]) #freq
as.numeric(summary(dual.logit.item.slope.m)$AIC[1]) #logit
as.numeric(summary(dual.MI.item.slope.m)$AIC[1]) #MI

#freq is best

#repeating the procedure assuming freq mod as null model for further random effect structure selection
#model with logit item slope nearly unidentifiable because of large eigenvalue ratio

dual.freq.MI.item.slope.m = lmer(RT.log ~ FreqABCD.log.std +
                                LogitABCD.neg.log.std +
                                MIABCD.neg.log.std +
                                AS.data$freq.sub.PC1 +
                                AS.data$freq.sub.PC2 +
                                AS.data$freq.sub.PC3 +
                                AS.data$freq.sub.PC4 +
                                block + 
                                nletter.std +
                                (1|subj) +
                                (FreqABCD.log.std|item) +
                                (MIABCD.log.std|item),
                              data = AS.data,
                              REML = FALSE)

anova(dual.freq.item.slope.m, dual.freq.MI.item.slope.m)
#no improvment, keeping freq item slope only

#rename best model
dual.best.m = dual.freq.item.slope.m

###########################

#dual model diagnostics

#checking for collinearity
#exclude intercept
collin.fnc(model.matrix(dual.best.m), 2:10)$cnumber
#13.05 medium collinearity

lmer.diag.mplot(dual.best.m)
#some leverage issues, similar to incremental model

#Fixed effects
dotplot(fixef(dual.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))

#Random effects
dotplot(ranef(dual.best.m,condVar=TRUE),
        lattice.options=list(layout=c(1,2)))


###################

#regression data extraction
write.csv(lmer.data.extract(dual.best.m), "AS.dual.m.orig.csv", row.names=TRUE) #fixed effects
write.csv(lmer.ranef.data.extract(dual.best.m), "AS.dual.m.ranef.orig.csv", row.names=FALSE) #random effects
write.csv(lmer.optim.data.extract(dual.best.m), "AS.dual.m.optim.orig.csv", row.names=TRUE)
write.csv(as.data.frame(dual.best.m@optinfo$derivs$Hessian), "AS.dual.m.hessian.orig.csv", row.names=TRUE)
write(unlist(dual.best.m@optinfo$conv$lme4$messages), "AS.dual.m.warnings.orig.txt")

###############################################################################################################
###############################################################################################################

# *** Regression bagging ****
setwd("D:/Q2_bootstraps/main_analysis/as_main")
source("lmer-data-extract-boot-fnc.R") #selected lmer results extractor

# divert messages stream to file, so you can log warnings and track non-converged models

options(warn=1)
wngs=file("warnings_log.txt",open="w+",blocking=TRUE)
sink(wngs,type="message")

for(iter in 1:3){

  setwd("D:/Q2_bootstraps/main_analysis/as_main")
  
  data = read.csv(paste("AS_bootstrap_sample_", iter,".csv", sep=""), header=TRUE)

  data$LogitABCD.neg.log.std =  -data$LogitABCD.log.std
  data$MIABCD.neg.log.std = -data$MIABCD.log.std
  
  #log model number to file so that any potential warnings appear underneath
  message(paste("holistic model #", iter, sep=""))
  
  holistic.mod = lmer(RT.log ~ FreqABCD.log.std +
                           freq.sub.PC1 +
                           freq.sub.PC2 +
                           freq.sub.PC3 +
                           freq.sub.PC4 +
                           block + 
                           nletter.std +
                           (1|subj) +
                           (FreqABCD.log.std|item),
                         data = data,
                         REML = FALSE)
  
  message(paste("probabilistic mod #", iter, sep=""))
  
  prob.mod = lmer(RT.log ~ LogitABCD.neg.log.std +
                           MIABCD.neg.log.std + 
                           block + 
                           nletter.std +
                           (1|subj) +
                           (MIABCD.neg.log.std|item),
                         data = data,
                         REML = FALSE)
  
  message(paste("dual mod #", iter, sep=""))
  
  dual.mod = lmer(RT.log ~ FreqABCD.log.std +
                           LogitABCD.neg.log.std +
                           MIABCD.neg.log.std +
                           freq.sub.PC1 +
                           freq.sub.PC2 +
                           freq.sub.PC3 +
                           freq.sub.PC4 +
                           block + 
                           nletter.std +
                           (1|subj) +
                           (FreqABCD.log.std|item),
                         data = data,
                         REML = FALSE)
  
  #write model results to file
  setwd("D:/Q2_bootstraps/main_analysis/as_main/reg_boot_holistic")
  write.csv(lmer.data.extract.boot(holistic.mod, iter), paste("holistic.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(holistic.mod, iter), paste("holistic.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(holistic.mod, iter), paste("holistic.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(holistic.mod@optinfo$derivs$Hessian), paste("holistic.m.hessian_", iter, ".csv", sep=""))
  write(unlist(holistic.mod@optinfo$conv$lme4$messages), paste("holistic.m.warnings_", iter, ".txt", sep=""))
  
  setwd("D:/Q2_bootstraps/main_analysis/as_main/reg_boot_prob")
  write.csv(lmer.data.extract.boot(prob.mod, iter), paste("prob.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(prob.mod, iter), paste("prob.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(prob.mod, iter), paste("prob.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(prob.mod@optinfo$derivs$Hessian), paste("prob.m.hessian_", iter, ".csv", sep=""))
  write(unlist(prob.mod@optinfo$conv$lme4$messages), paste("prob.m.warnings_", iter, ".txt", sep=""))

  setwd("D:/Q2_bootstraps/main_analysis/as_main/reg_boot_dual")
  write.csv(lmer.data.extract.boot(dual.mod, iter), paste("dual.m.fixef_", iter, ".csv", sep=""), row.names=TRUE) #fixed effects
  write.csv(lmer.ranef.data.extract.boot(dual.mod, iter), paste("dual.m.ranef_", iter, ".csv", sep=""), row.names=FALSE) #random effects
  write.csv(lmer.optim.data.extract.boot(dual.mod, iter), paste("dual.m.optim_", iter, ".csv", sep=""), row.names=TRUE)
  write.csv(as.data.frame(dual.mod@optinfo$derivs$Hessian), paste("dual.m.hessian_", iter, ".csv", sep=""))
  write(unlist(dual.mod@optinfo$conv$lme4$messages), paste("dual.m.warnings_", iter, ".txt", sep=""))
  
  
  cat("Bagging iteration", iter, "completed!\n")
  
}

#close log file & restore warnings stream to console
closeAllConnections()

##################################


