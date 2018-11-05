library("lme4")
library("multcomp")
library("contrast")
library("lsmeans")

sink("R_Exp2-Exp3.out")

SCAData <- read.table("E2-patient_data.txt")
SCAData$Group = "SCA"
CTLData <- read.table("E2-oldercontrol_data.txt")
CTLData$Group = "CTL"
GENData <- read.table("E3-patient_data.txt")
GENData$Group = "GEN"

YNGData <- read.table("E2-youngcontrol_data.txt")
YNGData$Group = "YNG"


dat <- rbind(SCAData, CTLData, YNGData)
dat$Bout <- factor(dat$Bout)
dat$Group <- factor(dat$Group)
dat$EarlyLate <- factor(dat$EarlyLate)
dat$Subj <- factor(dat$Subj)

EarlyData <- subset(dat, EarlyLate == 'Early')
LateData <- subset(dat, EarlyLate == 'Late')
WOData <- subset(dat, EarlyLate == 'Washout')

SCAEarlyData <- subset(EarlyData,Group == "SCA")
SCALateData <- subset(LateData,Group == "SCA")
SCAWOData <- subset(WOData,Group == "SCA")

CTLEarlyData <- subset(EarlyData,Group == "CTL")
CTLLateData <- subset(LateData,Group == "CTL")
CTLWOData <- subset(WOData,Group == "CTL")

YNGEarlyData <- subset(EarlyData,Group == "YNG")
YNGLateData <- subset(LateData,Group == "YNG")
YNGWOData <- subset(WOData,Group == "YNG")

LateDataAll <- LateData


cat("********************SCA vs CTL vs YNG Late********************\n\n")


adplate.model <- lmer(ReachDir ~ Bout + Bout:Group + Group + (1|Subj),data = LateData,REML=FALSE)

adplate.nInt <- update(adplate.model,formula = . ~ . - Bout:Group)
adplate.nGroup <- update(adplate.model,formula = . ~ . - Group - Bout:Group)
adplate.nIntnGroup <- update(adplate.nInt,formula = . ~ . - Group)
adplate.nBout <- update(adplate.model,formula = . ~ . - Bout - Bout:Group)
adplate.nIntnBout <- update(adplate.nInt,formula = . ~ . - Bout)


print(summary(adplate.model))

cat("\n\n>>>Interaction effect ANOVA: \n\n")
intanova = anova(adplate.model,adplate.nInt)
print(intanova)

if (intanova$Pr[2] < 0.05) {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova = anova(adplate.nGroup,adplate.model)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA: \n\n")
  nBoutanova = anova(adplate.nBout,adplate.model)
  print(nBoutanova)
  
  
  
} else {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA (no interaction): \n\n")
  nGroupanova = anova(adplate.nIntnGroup,adplate.nInt)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA (no interaction): \n\n")
  nBoutanova = anova(adplate.nIntnBout,adplate.nInt)
  print(nBoutanova)
  

}

cat("\n\n\n")

ld <- LateData[complete.cases(LateData),]
ld <- aggregate(ld,by=list(bout = ld$Bout,grp = ld$Group, subj = ld$Subj),FUN = mean)
ld.b1 <- subset(ld,bout == 'A1')
ld.b1.sca <- subset(ld.b1,grp == 'SCA')
ld.b1.ctl <- subset(ld.b1,grp == 'CTL')
ld.b1.yng <- subset(ld.b1,grp == 'YNG')
 
ld.b2 <- subset(ld,bout == 'A2')
ld.b2.sca <- subset(ld.b2,grp == 'SCA')
ld.b2.ctl <- subset(ld.b2,grp == 'CTL')
ld.b2.yng <- subset(ld.b2,grp == 'YNG')

ld.b3 <- subset(ld,bout == 'A3')
ld.b3.sca <- subset(ld.b3,grp == 'SCA')
ld.b3.ctl <- subset(ld.b3,grp == 'CTL')
ld.b3.yng <- subset(ld.b3,grp == 'YNG')

cat("Bout\tGrp\tMean\tSEM")
cat("\tGrp\tMean\tSEM")
cat("\tGrp\tMean\tSEM\n")

cat(sprintf("1\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b1.sca$ReachDir), sd(ld.b1.sca$ReachDir)/sqrt(length(ld.b1.sca$ReachDir)),
            "CTL",mean(ld.b1.ctl$ReachDir), sd(ld.b1.ctl$ReachDir)/sqrt(length(ld.b1.ctl$ReachDir)),
            "YNG",mean(ld.b1.yng$ReachDir), sd(ld.b1.yng$ReachDir)/sqrt(length(ld.b1.yng$ReachDir))))

cat(sprintf("2\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b2.sca$ReachDir), sd(ld.b2.sca$ReachDir)/sqrt(length(ld.b2.sca$ReachDir)),
            "CTL",mean(ld.b2.ctl$ReachDir), sd(ld.b2.ctl$ReachDir)/sqrt(length(ld.b2.ctl$ReachDir)),
            "YNG",mean(ld.b2.yng$ReachDir), sd(ld.b2.yng$ReachDir)/sqrt(length(ld.b2.yng$ReachDir))))

cat(sprintf("3\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b3.sca$ReachDir), sd(ld.b3.sca$ReachDir)/sqrt(length(ld.b3.sca$ReachDir)),
            "CTL",mean(ld.b3.ctl$ReachDir), sd(ld.b3.ctl$ReachDir)/sqrt(length(ld.b3.ctl$ReachDir)),
            "YNG",mean(ld.b3.yng$ReachDir), sd(ld.b3.yng$ReachDir)/sqrt(length(ld.b3.yng$ReachDir))))
cat("\n\n")

cat("Grp\tN\n")
cat(sprintf("%s\t%d\n","SCA",length(ld.b1.sca$subj)))
cat(sprintf("%s\t%d\n","CTL",length(ld.b1.ctl$subj)))
cat(sprintf("%s\t%d\n","YNG",length(ld.b1.yng$subj)))



cat("\n\n\n")

LateData$Group_Bout <- interaction(LateData$Group, LateData$Bout)

adplate.model <- lmer(ReachDir ~ 0 + Group_Bout + (1|Subj),data = LateData,REML=FALSE)

K <- t(matrix(c(
                1,-1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1,-1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1,-1, 0,
                1, 0,-1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0,-1,
                0, 1,-1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1,-1,
                1/3,  1/3,  1/3, -1/3, -1/3, -1/3,    0,    0,    0,
                  0,    0,    0,  1/3,  1/3,  1/3, -1/3, -1/3, -1/3,
                1/3,  1/3,  1/3,    0,    0,    0, -1/3, -1/3, -1/3,
                1, 0, 0,-1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0,-1, 0, 0,
                1, 0, 0, 0, 0, 0,-1, 0, 0,
                0, 1, 0, 0,-1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0,-1, 0,
                0, 1, 0, 0, 0, 0, 0,-1, 0,
                0, 0, 1, 0, 0,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0,-1,
                0, 0, 1, 0, 0, 0, 0, 0,-1),9))



rownames(K) <- c("ctlVSsca Bout1","ctlVSsca Bout2","ctlVSsca Bout3","ctlVSyng Bout1","ctlVSyng Bout2","ctlVSyng Bout3","yngVSsca Bout1","yngVSsca Bout2","yngVSsca Bout3", "Bout1 vs Bout2","Bout2 vs Bout3","Bout1 vs Bout3", "Ctl Bout1v2", "Ctl Bout2v3", "Ctl Bout1v3","SCA Bout1v2", "SCA Bout2v3", "SCA Bout1v3", "YNG Bout1v2", "YNG Bout2v3", "YNG Bout1v3")

posthocs <- glht(adplate.model,linfct = K)

print(summary(posthocs, test = adjusted('holm')))



cat("\n\n\n")



cat("********************SCA vs CTL vs YNG Early********************\n\n")


EarlyData$Group_Bout <- interaction(EarlyData$Group, EarlyData$Bout)

adpearly.model <- lmer(ReachDir ~ Bout + Bout:Group + Group + (1|Subj),data = EarlyData,REML=FALSE)
adpearly.nInt <- lmer(ReachDir ~ Bout + Group + (1|Subj),data = EarlyData,REML=FALSE)
adpearly.nGroup <- update(adpearly.model,formula = . ~ . - Group - Bout:Group)
adpearly.nIntnGroup <- update(adpearly.nInt,formula = . ~ . - Group - Bout:Group)
adpearly.nBout <- update(adpearly.model,formula = . ~ . - Bout - Bout:Group)
adpearly.nIntnBout <- update(adpearly.nInt,formula = . ~ . - Bout)


print(summary(adpearly.model))

cat("\n\n>>>Interaction effect ANOVA: \n\n")
intanova = anova(adpearly.model,adpearly.nInt)
print(intanova)

if (intanova$Pr[2] < 0.05) {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova = anova(adpearly.nGroup,adpearly.model)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA: \n\n")
  nBoutanova = anova(adpearly.nBout,adpearly.model)
  print(nBoutanova)
  
  
  
} else {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA (no interaction): \n\n")
  nGroupanova = anova(adpearly.nIntnGroup,adpearly.nInt)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA (no interaction): \n\n")
  nBoutanova = anova(adpearly.nIntnBout,adpearly.nInt)
  print(nBoutanova)
  
  
}

cat("\n\n\n")


ld <- EarlyData[complete.cases(EarlyData),]
ld <- aggregate(ld,by=list(bout = ld$Bout,grp = ld$Group, subj = ld$Subj),FUN = mean)
ld.b1 <- subset(ld,bout == 'A1')
ld.b1.sca <- subset(ld.b1,grp == 'SCA')
ld.b1.ctl <- subset(ld.b1,grp == 'CTL')
ld.b1.yng <- subset(ld.b1,grp == 'YNG')

ld.b2 <- subset(ld,bout == 'A2')
ld.b2.sca <- subset(ld.b2,grp == 'SCA')
ld.b2.ctl <- subset(ld.b2,grp == 'CTL')
ld.b2.yng <- subset(ld.b2,grp == 'YNG')

ld.b3 <- subset(ld,bout == 'A3')
ld.b3.sca <- subset(ld.b3,grp == 'SCA')
ld.b3.ctl <- subset(ld.b3,grp == 'CTL')
ld.b3.yng <- subset(ld.b3,grp == 'YNG')

cat("Bout\tGrp\tMean\tSEM")
cat("\tGrp\tMean\tSEM")
cat("\tGrp\tMean\tSEM\n")

cat(sprintf("1\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b1.sca$ReachDir), sd(ld.b1.sca$ReachDir)/sqrt(length(ld.b1.sca$ReachDir)),
            "CTL",mean(ld.b1.ctl$ReachDir), sd(ld.b1.ctl$ReachDir)/sqrt(length(ld.b1.ctl$ReachDir)),
            "YNG",mean(ld.b1.yng$ReachDir), sd(ld.b1.yng$ReachDir)/sqrt(length(ld.b1.yng$ReachDir))))

cat(sprintf("2\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b2.sca$ReachDir), sd(ld.b2.sca$ReachDir)/sqrt(length(ld.b2.sca$ReachDir)),
            "CTL",mean(ld.b2.ctl$ReachDir), sd(ld.b2.ctl$ReachDir)/sqrt(length(ld.b2.ctl$ReachDir)),
            "YNG",mean(ld.b2.yng$ReachDir), sd(ld.b2.yng$ReachDir)/sqrt(length(ld.b2.yng$ReachDir))))

cat(sprintf("3\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\t%s\t%.3f\t%.3f\n",
            "SCA",mean(ld.b3.sca$ReachDir), sd(ld.b3.sca$ReachDir)/sqrt(length(ld.b3.sca$ReachDir)),
            "CTL",mean(ld.b3.ctl$ReachDir), sd(ld.b3.ctl$ReachDir)/sqrt(length(ld.b3.ctl$ReachDir)),
            "YNG",mean(ld.b3.yng$ReachDir), sd(ld.b3.yng$ReachDir)/sqrt(length(ld.b3.yng$ReachDir))))
cat("\n\n")

cat("Grp\tN\n")
cat(sprintf("%s\t%d\n","SCA",length(ld.b1.sca$subj)))
cat(sprintf("%s\t%d\n","CTL",length(ld.b1.ctl$subj)))
cat(sprintf("%s\t%d\n","YNG",length(ld.b1.yng$subj)))

cat("\n\n\n")




adpearly.model <- lmer(ReachDir ~ 0 + Group_Bout + (1|Subj),data = EarlyData,REML=FALSE)


K <- t(matrix(c(1,-1, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1,-1, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1,-1, 0,
                1, 0,-1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 0,-1,
                0, 1,-1, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1,-1,
                1/3,  1/3,  1/3, -1/3, -1/3, -1/3,    0,    0,    0,
                0,    0,    0,  1/3,  1/3,  1/3, -1/3, -1/3, -1/3,
                1/3,  1/3,  1/3,    0,    0,    0, -1/3, -1/3, -1/3,
                1, 0, 0,-1, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0,-1, 0, 0,
                1, 0, 0, 0, 0, 0,-1, 0, 0,
                0, 1, 0, 0,-1, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0,-1, 0,
                0, 1, 0, 0, 0, 0, 0,-1, 0,
                0, 0, 1, 0, 0,-1, 0, 0, 0,
                0, 0, 0, 0, 0, 1, 0, 0,-1,
                0, 0, 1, 0, 0, 0, 0, 0,-1 ),9))

rownames(K) <- c("ctlVSsca Bout1","ctlVSsca Bout2","ctlVSsca Bout3","ctlVSyng Bout1","ctlVSyng Bout2","ctlVSyng Bout3","yngVSsca Bout1","yngVSsca Bout2","yngVSsca Bout3", "Bout1 vs Bout2","Bout2 vs Bout3","Bout1 vs Bout3", "Ctl Bout1v2", "Ctl Bout2v3", "Ctl Bout1v3","SCA Bout1v2", "SCA Bout2v3", "SCA Bout1v3","YNG Bout1v2", "YNG Bout2v3", "YNG Bout1v3")


posthocs <- glht(adpearly.model,linfct = K)


print(summary(posthocs, test = adjusted('holm')))





cat("\n\n\n")


cat("********************SCA Early********************\n\n")

adpearly.model <- lmer(ReachDir ~ Bout + (1|Subj),data = SCAEarlyData,REML=FALSE)
adpearly.base <- lmer(ReachDir ~ 1 + (1|Subj),data = SCAEarlyData,REML=FALSE)

print(summary(adpearly.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adpearly.model,adpearly.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {

cat("\n\n>>>Pairwise Testing: \n\n")

  #To do pairwise testing with glht, we have three choices. We can do a Tukey test (which controls the familywise error rate),
  # but we can only test main effects (i.e., no interaction terms). Second, we can manually input a contrast matrix to
  # test any possible terms, but it assumes a z test instead of a t test because glht doesn't know how to handle the dof.
  # in either case, it is possible to  specify the method for adjusting the p values after that. Third, we could use the lsmeans
  # version of glht and create the contrast matrix using lsm, but again it will result in a z score.
  #Note, I suspect that doing a Tukey test is a workaround and not correct; in particular, because summary() still defaults 
  # to performing a multiple comparisons correction despite having supposedly applied a Tukey test. Indeed, this seems to
  # be the case from reading the glht documentation, that it *should* be computing a z score.
  #It seems the best solution may be to use lsmeans directly.
  
  #paircomp <- glht(adpearly.model,linfct=mcp(Bout="Tukey"))
  #model.pairs <- contrast(lsmeans(adpearly.model,~Bout),"pairwise",adjust = "holm")
  #print(model.pairs)
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adpearly.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
  
}




cat("\n\n\n")

cat("********************SCA Late********************\n\n")

adplate.model <- lmer(ReachDir ~ Bout + (1|Subj),data = SCALateData,REML=FALSE)
adplate.base <- lmer(ReachDir ~ 1 + (1|Subj),data = SCALateData,REML=FALSE)

print(summary(adplate.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adplate.model,adplate.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adplate.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))

}


cat("\n\n\n")

cat("********************CTL Early********************\n\n")

adpearly.model <- lmer(ReachDir ~ Bout + (1|Subj),data = CTLEarlyData,REML=FALSE)
adpearly.base <- lmer(ReachDir ~ 1 + (1|Subj),data = CTLEarlyData,REML=FALSE)

print(summary(adpearly.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adpearly.model,adpearly.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adpearly.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
}


cat("\n\n\n")

cat("********************CTL Late********************\n\n")

adplate.model <- lmer(ReachDir ~ Bout + (1|Subj),data = CTLLateData,REML=FALSE)
adplate.base <- lmer(ReachDir ~ 1 + (1|Subj),data = CTLLateData,REML=FALSE)

print(summary(adplate.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adplate.model,adplate.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adplate.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
}


cat("\n\n\n")

cat("********************YNG Early********************\n\n")

adpearly.model <- lmer(ReachDir ~ Bout + (1|Subj),data = YNGEarlyData,REML=FALSE)
adpearly.base <- lmer(ReachDir ~ 1 + (1|Subj),data = YNGEarlyData,REML=FALSE)

print(summary(adpearly.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adpearly.model,adpearly.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adpearly.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
}

cat("\n\n\n")

cat("********************YNG Late********************\n\n")

adplate.model <- lmer(ReachDir ~ Bout + (1|Subj),data = YNGLateData,REML=FALSE)
adplate.base <- lmer(ReachDir ~ 1 + (1|Subj),data = YNGLateData,REML=FALSE)

print(summary(adplate.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adplate.model,adplate.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adplate.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
}

cat("\n\n\n")

cat("********************YNG vs CTL Late********************\n\n")

datc <- subset(dat,Group != "SCA")
datc <- subset(datc, EarlyLate == 'Late')

datc$Group_Bout <- interaction(datc$Group, datc$Bout)

adplate.model <- lmer(ReachDir ~ Bout + Bout:Group + Group + (1|Subj),data = datc,REML=FALSE)
adplate.nInt <- update(adplate.model, formula = . ~ . - Bout:Group)
adplate.nGroup <- update(adplate.model, formula = . ~ . - Group)
adplate.nIntnGroup <- update(adplate.nInt, formula = . ~ . - Group)
adplate.nBout <- update(adplate.model, formula = . ~ . - Bout)
adplate.nIntnBout <- update(adplate.nInt, formula = . ~ . - Bout)

print(summary(adplate.model))

cat("\n\n>>>Interaction effect ANOVA: \n\n")
intanova = anova(adplate.model,adplate.nInt)
print(intanova)

if (intanova$Pr[2] < 0.05) {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova = anova(adplate.nGroup,adplate.model)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA: \n\n")
  nBoutanova = anova(adplate.nBout,adplate.model)
  print(nBoutanova)

} else {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA (no interaction): \n\n")
  nGroupanova = anova(adplate.nIntnGroup,adplate.nInt)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA (no interaction): \n\n")
  nBoutanova = anova(adplate.nIntnBout,adplate.nInt)
  print(nBoutanova)
  
}

cat("\n\n\n")


adplate.model <- lmer(ReachDir ~ 0 + Group_Bout + (1|Subj),data = datc,REML=FALSE)

K <- t(matrix(c(1,-1, 0, 0, 0, 0,
                0, 0, 1,-1, 0, 0,
                0, 0, 0, 0, 1,-1,
                0.5, 0.5,-0.5,-0.5,   0,   0,
                0,   0, 0.5, 0.5,-0.5,-0.5,
                0.5, 0.5,   0,   0,-0.5,-0.5),6))

rownames(K) <- c("Grp Bout1","Grp Bout2","Grp Bout3","Bout1 vs Bout2","Bout2 vs Bout3","Bout1 vs Bout3") #  "Ctl Bout1v2", "Ctl Bout2v3", "Ctl Bout1v3","SCA Bout1v2", "SCA Bout2v3", "SCA Bout1v3"

posthocs <- glht(adplate.model,linfct = K)

print(summary(posthocs, test = adjusted('holm')))



cat("\n\n\n")

cat("********************CTL vs SCA Late********************\n\n")

datc <- subset(dat,Group != "YNG")
datc <- subset(datc, EarlyLate == 'Late')

datc$Group_Bout <- interaction(datc$Group, datc$Bout)

adplate.model <- lmer(ReachDir ~ Bout + Bout:Group + Group + (1|Subj),data = datc,REML=FALSE)
adplate.nInt <- update(adplate.model, formula = . ~ . - Bout:Group)
adplate.nGroup <- update(adplate.model, formula = . ~ . - Group)
adplate.nIntnGroup <- update(adplate.nInt, formula = . ~ . - Group)
adplate.nBout <- update(adplate.model, formula = . ~ . - Bout)
adplate.nIntnBout <- update(adplate.nInt, formula = . ~ . - Bout)

print(summary(adplate.model))

cat("\n\n>>>Interaction effect ANOVA: \n\n")
intanova = anova(adplate.model,adplate.nInt)
print(intanova)

if (intanova$Pr[2] < 0.05) {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova = anova(adplate.nGroup,adplate.model)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA: \n\n")
  nBoutanova = anova(adplate.nBout,adplate.model)
  print(nBoutanova)
  
} else {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA (no interaction): \n\n")
  nGroupanova = anova(adplate.nIntnGroup,adplate.nInt)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA (no interaction): \n\n")
  nBoutanova = anova(adplate.nIntnBout,adplate.nInt)
  print(nBoutanova)
  
}

cat("\n\n\n")


adplate.model <- lmer(ReachDir ~ 0 + Group_Bout + (1|Subj),data = datc,REML=FALSE)

K <- t(matrix(c(1,-1, 0, 0, 0, 0,
                0, 0, 1,-1, 0, 0,
                0, 0, 0, 0, 1,-1,
                0.5, 0.5,-0.5,-0.5,   0,   0,
                0,   0, 0.5, 0.5,-0.5,-0.5,
                0.5, 0.5,   0,   0,-0.5,-0.5),6))

rownames(K) <- c("Grp Bout1","Grp Bout2","Grp Bout3","Bout1 vs Bout2","Bout2 vs Bout3","Bout1 vs Bout3") #  "Ctl Bout1v2", "Ctl Bout2v3", "Ctl Bout1v3","SCA Bout1v2", "SCA Bout2v3", "SCA Bout1v3"

posthocs <- glht(adplate.model,linfct = K)

print(summary(posthocs, test = adjusted('holm')))


dat2 <- rbind(SCAData, GENData)
dat2$Bout <- factor(dat2$Bout)
dat2$Group <- factor(dat2$Group)
dat2$EarlyLate <- factor(dat2$EarlyLate)
dat2$Subj <- factor(dat2$Subj)
dat2$Bout_Group <- interaction(dat2$Bout,dat2$Group)

LateData2 <- subset(dat2, EarlyLate == 'Late')
WOData2 <- subset(dat2, EarlyLate == 'Washout')
EarlyData2 <- subset(dat2, EarlyLate == 'Early')

GENLateData <- subset(LateData2,Group == 'GEN')


cat("\n\n\n")

cat("********************SCA vs GEN Late********************\n\n")

adplate.model <- lmer(ReachDir ~ Bout + Bout:Group + Group + (1|Subj),data = LateData2,REML=FALSE)
adplate.nInt <- lmer(ReachDir ~ Bout + Group + (1|Subj),data = LateData2,REML=FALSE)
adplate.nGroup <- lmer(ReachDir ~ Bout + Bout:Group + (1|Subj),data = LateData2,REML=FALSE)
adplate.nIntnGroup <- lmer(ReachDir ~ Bout + (1|Subj),data = LateData2,REML=FALSE)
adplate.nBout <- lmer(ReachDir ~ Bout:Group + Group + (1|Subj),data = LateData2,REML=FALSE)
adplate.nIntnBout <- lmer(ReachDir ~ Group + (1|Subj),data = LateData2,REML=FALSE)


print(summary(adplate.model))

cat("\n\n>>>Interaction effect ANOVA: \n\n")
intanova = anova(adplate.model,adplate.nInt)
print(intanova)

if (intanova$Pr[2] < 0.05) {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova = anova(adplate.nGroup,adplate.model)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA: \n\n")
  nBoutanova = anova(adplate.nBout,adplate.model)
  print(nBoutanova)
  
  
  
} else {
  
  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA (no interaction): \n\n")
  nGroupanova = anova(adplate.nIntnGroup,adplate.nInt)
  print(nGroupanova)
  
  cat("\n>>>Bout effect ANOVA (no interaction): \n\n")
  nBoutanova = anova(adplate.nIntnBout,adplate.nInt)
  print(nBoutanova)
  
  
}

cat("\n\n\n")



adplate.model <- lmer(ReachDir ~ 0 + Bout_Group + (1|Subj),data = LateData2,REML=FALSE)


cat("\n\n\n>>>Blockwise posthoc paired testing: \n\n")
K <- t(matrix(c(1, 0, 0, -1, 0, 0,
                0, 1, 0, 0,-1, 0,
                0, 0, 1, 0, 0,-1),6))
rownames(K) <- c("A1 vs A2", "A1 vs A3","A2 vs A3")

posthocs <- glht(adplate.model,linfct = K)


print(summary(posthocs, test = adjusted('holm')))


cat("\n\n-----Retention tests-----\n\n")

B1dat <- subset(LateData2,Bout == "A1")
B1dat <- subset(LateData2,Bout == "A1")
B1dat <- subset(B1dat,Subj != "S5" & Subj != "S10" & Subj != "S11" & Subj != "S14" & Subj != "S15")

B1B3dat <- subset(LateData2,(Bout_Group == "A3.SCA") | (Bout_Group == "A1.GEN"))
B1B3dat <- subset(B1B3dat,Subj != "S5" & Subj != "S10" & Subj != "S11" & Subj != "S14" & Subj != "S15")



cat("~~~Bout 1 comparison~~~\n\n")

B1mod <- lmer(ReachDir ~ Group + (1|Subj), data = B1dat)
B1mod.nGrp <- lmer(ReachDir ~ 1 + (1|Subj), data = B1dat)

print(summary(B1mod))

cat("\n\n---\n\n")

print(anova(B1mod,B1mod.nGrp))

cat("\n~~~Bout 1-Exp3 vs Bout3-Exp2 comparison~~~\n\n")

B31mod <- lmer(ReachDir ~ Group + (1|Subj), data = B1B3dat)
B31mod.nGrp <- lmer(ReachDir ~ 1 + (1|Subj), data = B1B3dat)

print(summary(B31mod))

cat("\n\n---\n\n")

print(anova(B31mod,B31mod.nGrp))



cat("\n\n\n")

cat("********************GEN Late********************\n\n")

adplate.model <- lmer(ReachDir ~ Bout + (1|Subj),data = GENLateData,REML=FALSE)
adplate.base <- lmer(ReachDir ~ 1 + (1|Subj),data = GENLateData,REML=FALSE)

print(summary(adplate.model))

cat("\n\n>>>Bout effect ANOVA: \n\n")
boutanova = anova(adplate.model,adplate.base)
print(boutanova)

if (boutanova$Pr[2] < 0.05) {
  
  cat("\n\n>>>Pairwise Testing: \n\n")
  contrast.matrix <- rbind("A1 vs A2" = c(0,1,0),
                           "A1 vs A3" = c(0,0,1),
                           "A2 vs A3" = c(0,1,-1))
  
  paircomp <- glht(adplate.model,linfct=contrast.matrix)
  print(summary(paircomp,test = adjusted("holm")))
  
}



sink()