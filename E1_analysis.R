library("lme4")
library("multcomp")

sink("R_Exp1.out")


DataAE <- read.table("E1-AEblocks_alternatevis.txt",na.strings = "NaN")
DataADP <- read.table("E1-adaptblocks_alternatevis.txt",na.strings = "NaN")
DataAEC <- read.table("E1-AEblocks_allvis.txt",na.strings = "NaN") #control group
DataADPC <- read.table("E1-adaptblocks_allvis.txt",na.strings = "NaN") #control group


cat("**************************** AE ReachDir ****************************\n\n")

DataAE$Group <- (rep(0,nrow(DataAE)))
DataAEC$Group <- (rep(1,nrow(DataAEC)))
DataAEAll <- rbind(DataAE,DataAEC)
DataAEAll$Block <- (DataAEAll$Drift + 2*DataAEAll$PostDrift)
DataAEAll$BlockF <- factor(DataAEAll$Drift + 2*DataAEAll$PostDrift)
DataAEAll$Adapt <- 1-DataAEAll$Drift
DataAEAll$Block_Group <- interaction(DataAEAll$BlockF, DataAEAll$Group)


cat("Full Model: ReachDir ~ Group*Block, (random effect of subject)\n\n")

AE.model <- lmer(ReachDir ~ Group*BlockF + (1|Subj), data = DataAEAll,REML = FALSE)
AE.nInt <- lmer(ReachDir ~ Group + BlockF + (1|Subj), data = DataAEAll,REML = FALSE)
AE.nGrp <- lmer(ReachDir ~ Group:Block + BlockF + (1|Subj), data = DataAEAll,REML = FALSE)
AE.nBlk <- lmer(ReachDir ~ Group + Group:BlockF + (1|Subj), data = DataAEAll,REML = FALSE)


print(summary(AE.model))


cat("\n\n>>>Interaction ANOVA: \n\n")
interactanova <- anova(AE.nInt,AE.model)
print(interactanova)

  cat("\n----------\n")
  
  cat("\n>>>Group effect ANOVA: \n\n")
  nGroupanova <- anova(AE.nGrp,AE.model)
  print(nGroupanova)
  
  cat("\n\n>>>Block effect ANOVA: \n\n")
  nBlockanova <- anova(AE.nBlk,AE.model)
  print(nBlockanova)
  

AE.model <- lmer(ReachDir ~ 0 + Block_Group + (1|Subj), data = DataAEAll,REML = FALSE)
  
  

cat("\n\n\n>>>Blockwise posthoc paired testing: \n\n")
K <- t(matrix(c(1, 0, 0,-1, 0, 0,
                0, 1, 0, 0,-1, 0,
                0, 0, 1, 0, 0,-1,
                0.5,-0.5,   0, 0.5,-0.5,   0,
                  0, 0.5,-0.5,   0, 0.5,-0.5,
                0.5,   0,-0.5, 0.5,   0,-0.5,
                  0, 1,-1, 0, 0, 0,
                  0, 0, 0, 0, 1,-1),6))
rownames(K) <- c("Grp AE1","Grp AE2","Grp AE3","AE1 vs AE2","AE2 vs AE3","AE1 vs AE3","AltVis AE2v3","AllVis AE2v3")

posthocs <- glht(AE.model,linfct = K)


print(summary(posthocs, test = adjusted('holm')))

  
cat("\n\n\n")



cat("**************************** ADAPT ReachDir ****************************\n\n")


DataADP$Group <- rep(0,nrow(DataADP))
DataADPC$Group <- rep(1,nrow(DataADPC))
DataADPAll <- rbind(DataADP,DataADPC)
#DataADPAll$Block <- DataADPAll$Drift + 2*DataADPAll$PostDrift
DataADPAll$Block <- DataADPAll$Adapt + 2*DataADPAll$Drift + 3*DataADPAll$PostDrift + 4*DataADPAll$Washout
DataADPAll$BlockF <- factor(DataADPAll$Block)
DataADPAll$Block_Group <- interaction(DataADPAll$BlockF, DataADPAll$Group)
DataADPAllnBase <- subset(DataADPAll,DataADPAll$Block > 0)
DataADPAllnBase$Block <- DataADPAllnBase$Block-1
DataADPAllnBase$BlockF <- factor(DataADPAllnBase$Block)
DataADPAllnBase$Block_Group <- interaction(DataADPAllnBase$BlockF, DataADPAllnBase$Group)

cat("Full Model: ReachDir ~ Group*Block, (random effect of subject)\n\n")


adp.model <- lmer(ReachDir ~ Group*BlockF + (1|Subj), data = DataADPAll,REML = FALSE)
adp.nInt <- lmer(ReachDir ~ Group + BlockF + (1|Subj), data = DataADPAll,REML = FALSE)
adp.nGrp <- lmer(ReachDir ~ Group:Block + BlockF + (1|Subj), data = DataADPAll,REML = FALSE)
adp.nBlk <- lmer(ReachDir ~ Group + Group:BlockF + (1|Subj), data = DataADPAll,REML = FALSE)


print(summary(adp.model))


cat("\n\n>>>Interaction ANOVA: \n\n")
interactanova <- anova(adp.nInt,adp.model)
print(interactanova)

cat("\n----------\n")

cat("\n>>>Group effect ANOVA: \n\n")
nGroupanova <- anova(adp.nGrp,adp.model)
print(nGroupanova)

cat("\n\n>>>Block effect ANOVA: \n\n")
nBlockanova <- anova(adp.nBlk,adp.model)
print(nBlockanova)


adp.model <- lmer(ReachDir ~ 0 + Block_Group + (1|Subj), data = DataADPAll,REML = FALSE)


cat("\n\n\n>>>Blockwise posthoc paired testing: \n\n")
K <- t(matrix(c(1, 0, 0, 0, 0,-1, 0, 0, 0, 0,
                0, 1, 0, 0, 0, 0,-1, 0, 0, 0,
                0, 0, 1, 0, 0, 0, 0,-1, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0,-1, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0,-1,
                0.5,-0.5, 0, 0, 0, 0.5,-0.5, 0, 0, 0),10))
rownames(K) <- c("Grp Baseline","Grp A1","Grp A2","Grp A3","Grp Washout","Base vs A1")

posthocs <- glht(adp.model,linfct = K)


print(summary(posthocs, test = adjusted('holm')))



cat("\n\n\n")





sink()
