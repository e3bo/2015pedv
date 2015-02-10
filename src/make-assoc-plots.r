library(car)
library(DescTools)
library(vcd)
load('mgAge.RData')
set.seed(3234)

#' Test for conditional independence
#+
ct <- xtabs(cases~State+age, data=mgAge)
test <- rowSums(ct) != 0
ct <- as.table(ct[test,])
colnames(ct) <- c('Grower/Finisher', 'Nursery', 'Sow/Boar', 'Suckling')
names(dimnames(ct))[[2]] <- 'Age Class'
coindep_test(ct, indepfun=function(x) sum(x^2))
assoc(ct, shade=TRUE, labeling_args=list(rot_labels=c(0,0,90,0),
                          offset_varnames=c(0,0,0,0.5)))
mosaic(ct, shade=TRUE, labeling_args=list(rot_labels=c(0,0,90,0),
                          offset_varnames=c(0,0,0,0.5)))
pdf('assoc-plots.pdf', width=6.8, height=7.8)
assoc(ct, shade=TRUE, labeling_args=list(rot_labels=c(0,0,90,0),
                          offset_varnames=c(0,0,0,0.5)))
dev.off()
Desc(ct)
Assocs(ct)

#' See if model probablities of suckling cases are predictive of observed.
#'
#' First, let's get the data into the right format.
#+
mgAge$isSuckling <- mgAge$age == 'Suckling'
cs <- xtabs(cases~State+isSuckling, data=mgAge)
test <- rowSums(cs) != 0
cs <- cs[test,]
cs <- cs[,c(2,1)]
wt <- xtabs(samplingWeight~State + isSuckling, data=mgAge)
wt <- wt/rowSums(wt)
wt <- wt[test,]

#' Do empirical and theoretical probabilities seem to be correlated?
#+
pEmp <- cs[,1]/rowSums(cs)
pThe <- wt[, 'TRUE']
plot(x=pThe,y=pEmp, type='n')
text(x=pThe,y=pEmp, labels=rownames(wt))

#' Not particularly. Let's address this more carefully with a logistic
#' regression.
#+
m <- glm(cs~1 + pThe, family=binomial)
summary(m)
par(mfrow=c(2,2))
plot(m)
outlierTest(m)

#' There's no effect, but the diagnostics indicate that NC, OK, and IA
#' are having large effects on the fit, and IA and NC are
#' outliers. Let's refit without IA and NC.
#+
test2 <- rownames(cs) %nin% c('IA', 'NC')
cs2 <- cs[test2,]
pThe2 <- pThe[test2]
m2 <- glm(cs2~1 + pThe2, family=binomial)
summary(m2)
par(mfrow=c(2,2))
plot(m2)
outlierTest(m2)

#' Still no effect, and now KS, OK, and MN have the highest Cook's
#' distance. Let's refit without them.
#+
test3 <- rownames(cs) %nin% c('IA', 'NC', 'KS', 'OK', 'MN')
cs3 <- cs[test3,]
pThe3 <- pThe[test3]
m3 <- glm(cs3~1 + pThe3, family=binomial)
summary(m3)
par(mfrow=c(2,2))
plot(m3)
outlierTest(m3)

#' IL and CO have somewhat large Cook's distances. Let's refit without
#' them.
#+
test4 <- rownames(cs) %nin% c('IA', 'NC', 'KS', 'OK', 'MN', 'CO', 'IL')
cs4 <- cs[test4,]
pThe4 <- pThe[test4]
m4 <- glm(cs4~1 + pThe4, family=binomial)
summary(m4)
par(mfrow=c(2,2))
plot(m4)
outlierTest(m4)

#' No significan effect, but we've taken out most of the states with
#' lot's of cases. Perhaps they could be fit to their own model.
#'
#+
test5 <- rownames(cs) %in% c('IA', 'NC', 'KS', 'OK', 'MN', 'CO', 'IL')
cs5 <- cs[test5,]
pThe5 <- pThe[test5]
m5 <- glm(cs5~1 + pThe5, family=binomial)
summary(m5)
par(mfrow=c(2,2))
plot(m5)
outlierTest(m5)

#' The diagnostics don't look particularly good, but at least this
#' shows which directions the major states deviate from predictions
#' in.
#'
#' Let's get an estimate of the average difference between predicted
#' and observed logits.
#+
offset <- qlogis(pThe)
mm <- glm(cs~1, offset=offset, family=binomial)
summary(mm)
exp(coef(mm))
exp(confint(mm))
par(mfrow=c(2,2))
plot(mm)
cbind(pEmp, pFit=predict(mm, type='response'), pThe)

#' The logit is about 4 times higher than expected on average. More so
#' in IA and NC. In OK, the empirical and theoretical probabilities
#' are close, but that gives OK a large negative residual in the
#' fitted model.
#'
#' Now let's refit without the IA, NC, and OK.
#+
test2m <- rownames(cs) %nin% c('IA', 'NC', 'OK')
csm <- cs[test2m,]
offsetm <- offset[test2m]
mm2 <- glm(csm~1, offset=offsetm, family=binomial)
summary(mm2)
exp(coef(mm))
exp(confint(mm2))
par(mfrow=c(2,2))
plot(mm2)

#' No major changes, but KS and IL look influential now. Let's check
#' their influence.
#'
#+
test3m <- rownames(cs) %nin% c('IA', 'NC', 'OK', 'KS', 'IL')
csm <- cs[test3m,]
offsetm <- offset[test3m]
mm3 <- glm(csm~1, offset=offsetm, family=binomial)
summary(mm3)
exp(coef(mm3))
exp(confint(mm3))
par(mfrow=c(2,2))
plot(mm3)

#' The difference is starting to look smaller, and MN now has the
#' largest Cook's distance. Let check it's influence.
#+
test4m <- rownames(cs) %nin% c('IA', 'NC', 'OK', 'KS', 'IL', 'MN')
csm <- cs[test4m,]
offsetm <- offset[test4m]
mm4 <- glm(csm~1, offset=offsetm, family=binomial)
summary(mm4)
exp(coef(mm4))
exp(confint(mm4))
par(mfrow=c(2,2))
plot(mm4)

#' We've moved down in size a little again, but there are no
#' observations with large Cook's distances now, and there's still a
#' clear positive difference between the means.
