library("lme4")
library("multcomp")
library("contrast")
library("lsmeans")
library("MASS")
library("car")

sink("R_adppsychR1.out")

psychdata <- read.table("adpscavis_psychR1.txt",na.strings = "NaN")
psychdatactl <- read.table("adpctlvis_psych2.txt",na.strings = "NaN")

psychdata$Subj <- factor(psychdata$Subj)
psychdatactl$Subj <- factor(psychdatactl$Subj)


psychdata$DigitSpan <- (psychdata$DigitBack + psychdata$DigitFwd)/2
psychdata$dA1A2 <- psychdata$A2-psychdata$A1
psychdata$dA1A3 <- psychdata$A3-psychdata$A1
psychdata$dA2A3 <- psychdata$A3-psychdata$A2
psychdata$dR1R2 <- psychdata$R2-psychdata$R1
psychdata$dR1R3 <- psychdata$R3-psychdata$R1
psychdata$dR2R3 <- psychdata$R3-psychdata$R2


psychdatactl$DigitSpan <- (psychdatactl$DigitBack + psychdatactl$DigitFwd)/2
psychdatactl$dA1A2 <- psychdatactl$A2-psychdatactl$A1
psychdatactl$dA1A3 <- psychdatactl$A3-psychdatactl$A1
psychdatactl$dA2A3 <- psychdatactl$A3-psychdatactl$A2
psychdatactl$dR1R2 <- psychdatactl$R2-psychdatactl$R1
psychdatactl$dR1R3 <- psychdatactl$R3-psychdatactl$R1
psychdatactl$dR2R3 <- psychdatactl$R3-psychdatactl$R2


cat("********************SCA Psych Data********************\n\n")

cat(">>>Adapt 1: Stepwise Regression\n\n")

psychdata <- subset(psychdata,!is.na(DigitFwd))

adp1.model <- lm(A1 ~ DigitSpan + ReyO + ICARS,data = psychdata)  #data = subset(psychdata,!is.na(DigitFwd))  #ReyO_Strategy +
adp1.step <- stepAIC(adp1.model,trace=0, direction = "both")

print(adp1.step$anova)
cat("\n\n")
print(summary(adp1.step))


cat("\n\n~~~~~~~~~~\n\n")
cat(">>>Adapt 2: Stepwise Regression\n\n")

adp2.model <- lm(A2 ~ A1 + DigitSpan + ReyO + ICARS,data = psychdata)
adp2.step <- stepAIC(adp2.model,trace=0, direction = "both")
print(adp2.step$anova)
cat("\n\n")
print(summary(adp2.step))

adp2.finalmodel <- lm(A2 ~ DigitSpan + ReyO, data = psychdata)
avPlots(adp2.finalmodel)


cat("\n\n~~~~~~~~~~\n\n")
cat(">>>Adapt 3: Stepwise Regression\n\n")

adp3.model <- lm(A3 ~ A1 + dA1A2 + DigitSpan + ReyO + ICARS,data = psychdata)

adp3.step <- stepAIC(adp3.model,trace=0, direction = "both")
print(adp3.step$anova)
cat("\n\n")
print(summary(adp3.step))

adp3.finalmodel <- lm(A3 ~ A1 + dA1A2 + DigitSpan, data = psychdata)
avPlots(adp3.finalmodel)




cat("\n\n\n\n\n")



cat("********************CTL Psych Data********************\n\n")

psychdatactl <- subset(psychdatactl,!is.na(DigitFwd))

cat(">>>Adapt 1: Stepwise Regression\n\n")

adp1.model <- lm(A1 ~ DigitSpan + ReyO,data = psychdatactl)
adp1.step <- stepAIC(adp1.model,trace=0, direction = "both")

print(adp1.step$anova)
cat("\n\n")
print(summary(adp1.step))


cat("\n\n~~~~~~~~~~\n\n")
cat(">>>Adapt 2: Stepwise Regression\n\n")

adp2.model <- lm(A2 ~ A1 + DigitSpan + ReyO,data = psychdatactl)
adp2.step <- stepAIC(adp2.model,trace=0, direction = "both")
print(adp2.step$anova)
cat("\n\n")
print(summary(adp2.step))

adp2ctl.finalmodel <- lm(A2 ~ DigitSpan, data = psychdatactl)
avp2c <- avPlots(adp2ctl.finalmodel)


cat("\n\n~~~~~~~~~~\n\n")
cat(">>>Adapt 3: Stepwise Regression\n\n")

adp3.model <- lm(A3 ~ A1 + dA1A2 + DigitSpan + ReyO,data = psychdatactl)
adp3.step <- stepAIC(adp3.model,trace=0, direction = "both")
print(adp3.step$anova)
cat("\n\n")
print(summary(adp3.step))

adp3ctl.finalmodel <- lm(A3 ~ A1  + DigitSpan, data = psychdatactl)
avp3c <- avPlots(adp3ctl.finalmodel)


cat("\n\n\n\n\n")


sink()