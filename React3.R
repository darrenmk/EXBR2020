setwd('D:\\Data\\Documents\\Research\\Tutorial\\')

subjects <- read.csv('React3Subjects.csv')
codes <- read.csv('React3Matrix.csv')

# Read data into a csv file
data_append <- {}

for (s in 1:dim(subjects)[1]){
  
  #print(s)
  file_append <- {}
  
  for (b in 1:4){
    
    filename <- paste('RawData\\react', subjects$sID[s], '_0', b, '_0.MLOG', sep="")
    #print(filename)
    
    file0 = read.table(filename, sep="\t", header=FALSE, fill=TRUE)
    
    # this is just some processing
    goodTrial1 <- file0[,3] != -2 # value of third column is NOT -2
    goodTrial2 <- is.na(file0[,6]) # is NA in column 6
    goodTrial3 <- !is.na(file0[,5]) # is NOT NA in column 5
    goodTrials <- goodTrial1 & goodTrial2 & goodTrial3
    file1 <- file0[goodTrials,]
    
    # make our data file
    file <- matrix(0,dim(file1)[1],5)
    file[,1] <- s
    file[,2] <- b
    file[,3] <- file1[,4] # this is our motion condition
    
    soa <- (file1[,5]+4)*(1000/60) # fame correction for motion simulator, this is our SOA
    
    file[,4] <- soa
    file[soa == -100,5] <- file1[soa < 0,1] + 100 # adjust reaction time
    file[soa >= 0,5] <- file1[soa >= 0,1] #keep reaction time
    
    file_append <- rbind(file_append,file)
    
  }
  
  data_append <- rbind(data_append,file_append)
}

colnames(data_append) <- c("Subject","Block","Condition","SOA","RT")

newCode <- rep(0,dim(data_append)[1])
for (i in 1:dim(data_append)[1]){
  newCode[i] <- which(data_append[i,3]==codes$code)
}

data_append[,3] <- newCode


## simple outlier removal

allOutliers <- as.logical(rep(FALSE,dim(data_append)[1])) # first make this all FALSE
for (o in 1:27){
  
  idx <- data_append[,3]==o
  
  logRT <- log(data_append[idx,5]) # index by condition number, get RT then calculate log
  Q1 <- quantile(logRT)[2]
  Q3 <- quantile(logRT)[4]
  IQR <- Q3-Q1
  
  lowOutlier <- logRT < Q1 - 3*IQR
  highOutlier <- logRT > Q3 + 3*IQR
  condOutliers <- lowOutlier | highOutlier
  
  idx[idx] <- condOutliers # idx now includes outliers for this condition only
  
  allOutliers <- allOutliers | idx # append after each loop
}

data <- data.frame(data_append[!allOutliers,])

grandMean <- tapply(data$RT,data$Condition,mean)
grandSd <- with(data,tapply(RT,Condition,sd))
grandSe <- grandSd / sqrt(35)
subjectMean <- tapply(data$RT,list(data$Subject,data$Condition),mean)
subjectBlockMean <- with(data,tapply(RT,list(Subject,Block,Condition),mean))

##### Now we do the statistics

options(contrasts=c("contr.sum","contr.poly") )

library(car) # for the repeated measures ANOVA
library(effects) # for the ANCOVA
library(tidyverse) # for the ggplot figures
library(ggpubr) # for arranging the panels

## Stats for Figure 2: SOA x Accel for each Vi (repeated measures ANOVA)
grandSoa <- c(rep("vestibular",9),rep("sync",9),rep("visual",9))
grandSoa <- factor(grandSoa,levels=c("vestibular","sync","visual"))
grandVi <- as.factor(rep(c(rep("3",3),rep("9",3),rep("15",3)),3))
grandVi <- factor(grandVi,levels=c(3,9,15)) # need to manually order
grandVf <- as.factor(rep(c(rep(c("21","27","33"),3)),3))
grandAccel <- as.factor(rep(c("54","72","90","36","54","72","18","36","54"),3))
grandAccel <- factor(grandAccel,levels=c("18","36","54","72","90")) # need to manually order
grandAccelNum <- as.numeric(as.character(grandAccel)) # convert to numeric values to treat as a continuous covariate
grandViSoa <- interaction(grandVi,grandSoa) # combine to make 9 levels / conditions

grandReact3 <- data.frame(grandMean,grandSe,grandSoa,grandVi,grandVf,grandAccel,grandAccelNum,grandViSoa)

# Vi=3, Accel=54,72,90
grandVi3Soa <- as.factor(c(rep("vestibular",3),rep("sync",3),rep("visual",3)))
grandVi3Accel <- as.factor(rep(c("54","72","90"),3))
grandVi3Factors <- interaction(grandVi3Soa,grandVi3Accel)
grandVi3iData <- data.frame(grandVi3Soa,grandVi3Accel)
grandVi3iDesign <- ~grandVi3Soa*grandVi3Accel
grandVi3Columns <- c(1L,2L,3L,1+9L,2+9L,3+9L,1+18L,2+18L,3+18L)
grandVi3SubjectMean <- subjectMean[,grandVi3Columns]
colnames(grandVi3SubjectMean) <- grandVi3Factors
grandVi3Mlm <- lm(grandVi3SubjectMean~1)
grandVi3Aov <- Anova(grandVi3Mlm,idata=grandVi3iData,idesign=grandVi3iDesign,type="III")
summary(grandVi3Aov,multivariate=FALSE)

# Vi=9, Accel=36,54,72
grandVi9Soa <- as.factor(c(rep("vestibular",3),rep("sync",3),rep("visual",3)))
grandVi9Accel <- as.factor(rep(c("36","54","72"),3))
grandVi9Factors <- interaction(grandVi9Soa,grandVi9Accel)
grandVi9iData <- data.frame(grandVi9Soa,grandVi9Accel)
grandVi9iDesign <- ~grandVi9Soa*grandVi9Accel
grandVi9Columns <- c(4L,5L,6L,4+9L,5+9L,6+9L,4+18L,5+18L,6+18L)
grandVi9SubjectMean <- subjectMean[,grandVi9Columns]
colnames(grandVi9SubjectMean) <- grandVi9Factors
grandVi9Mlm <- lm(grandVi9SubjectMean~1)
grandVi9Aov <- Anova(grandVi9Mlm,idata=grandVi9iData,idesign=grandVi9iDesign,type="III")
summary(grandVi9Aov,multivariate=FALSE)

# Vi=15, Accel=18,36,54
grandVi15Soa <- as.factor(c(rep("vestibular",3),rep("sync",3),rep("visual",3)))
grandVi15Accel <- as.factor(rep(c("18","36","54"),3))
grandVi15Factors <- interaction(grandVi15Soa,grandVi15Accel)
grandVi15iData <- data.frame(grandVi15Soa,grandVi15Accel)
grandVi15iDesign <- ~grandVi15Soa*grandVi15Accel
grandVi15Columns <- c(7L,8L,9L,7+9L,8+9L,9+9L,7+18L,8+18L,9+18L)
grandVi15SubjectMean <- subjectMean[,grandVi15Columns]
colnames(grandVi15SubjectMean) <- grandVi15Factors
grandVi15Mlm <- lm(grandVi15SubjectMean~1)
grandVi15Aov <- Anova(grandVi15Mlm,idata=grandVi15iData,idesign=grandVi15iDesign,type="III")
summary(grandVi15Aov,multivariate=FALSE)

## Figure 2: SOA x Accel for each Vi
grandVi3 <- ggplot(grandReact3 %>% filter(grandVi==3),aes(x=grandAccel,y=grandMean,group=grandSoa,colour=grandSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=grandSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=grandMean-grandSe, ymax=grandMean+grandSe), width=0.1, size=1) +
  labs(x = expression(paste("Acceleration (m/s" ^2, ")")), y = "Response Time (ms)") + 
  scale_y_continuous(breaks = seq(560, 760, by = 20), limits=c(560,860)) +
  guides(linetype=FALSE,shape=FALSE,colour=FALSE) +
  annotate("text",x=2,y=820,label="Initial Velocity: 3m/s", size=5) +
  theme_classic()

grandVi9 <- ggplot(grandReact3 %>% filter(grandVi==9),aes(x=grandAccel,y=grandMean,group=grandSoa,colour=grandSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=grandSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=grandMean-grandSe, ymax=grandMean+grandSe), width=0.1, size=1) +
  labs(x = expression(paste("Acceleration (m/s" ^2, ")")), y = "Response Time (ms)") + 
  scale_y_continuous(breaks = seq(560, 760, by = 20), limits=c(560,860)) +
  guides(linetype=FALSE,shape=FALSE,colour=FALSE) +
  annotate("text",x=2,y=820,label="Initial Velocity: 9m/s", size=5) +
  theme_classic()

grandVi15 <- ggplot(grandReact3 %>% filter(grandVi==15),aes(x=grandAccel,y=grandMean,group=grandSoa,colour=grandSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=grandSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=grandMean-grandSe, ymax=grandMean+grandSe), width=0.1, size=1) +
  labs(x = expression(paste("Acceleration (m/s" ^2, ")")), y = "Response Time (ms)") + 
  scale_colour_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  scale_shape_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  scale_y_continuous(breaks = seq(560, 760, by = 20), limits=c(560,860)) +
  annotate("text",x=2,y=820, label="Initial Velocity: 15m/s", size=5) +
  theme_classic()

## This is Figure 2 (save manually as .eps to avoid pixellation)
ggarrange(grandVi3, grandVi9, grandVi15, labels = c("A", "B", "C"),
          common.legend = TRUE, legend="top", nrow = 1, ncol = 3
)

## ---------------------------------------------------------- ##

## Stats for Figure 3: SOA x Vi with Acceleration as the covariate (ANCOVA)
ancovaCondition <- interaction(grandVi,grandSoa) # combine to make 9 levels / conditions

ancovaLM <- lm(grandMean~grandAccelNum+ancovaCondition) # use ancovaConditions to easily get alpha coefficients
coef(ancovaLM)
ancovaIntercept <- summary(ancovaLM)$coefficients[1]
ancovaBeta <- summary(ancovaLM)$coefficients[2]
ancovaAlphas <- c(summary(ancovaLM)$coefficients[3:10],0-sum(summary(ancovaLM)$coefficients[3:10])) 
ancovaTukey <- aov(grandMean~grandAccelNum+ancovaCondition)
TukeyHSD(ancovaTukey,which="ancovaCondition") # tukeyHSD to evaluate pairwise differences

ancovaResults <- lm(grandMean~grandAccelNum+grandSoa*grandVi) # ancova model
anova(ancovaResults) # anova table, report ancova and covariate (condition x rt): F, p

ancovaHomogeneity <- lm(grandMean~grandAccelNum+ancovaCondition+grandAccelNum:grandVi:grandSoa)
anova(ancovaHomogeneity) # report interaction: F, p ... if not significant = assumption is met

## Figure 3: SOA x Vi
ancovaPanelA <- ggplot(grandReact3, aes(x=grandAccelNum, y=grandMean, group=grandSoa, colour=grandSoa, linetype=grandVi)) +
  geom_point(position=position_jitter(0.5), aes(shape=grandSoa), size=4) +
  geom_vline(xintercept=54, colour="gray40", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[1], slope=ancovaBeta, colour="#F8766D", linetype="dashed", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[2], slope=ancovaBeta, colour="#F8766D", linetype="dotdash", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[3], slope=ancovaBeta, colour="#F8766D", linetype="dotted", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[4], slope=ancovaBeta, colour="#00BA38", linetype="dashed", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[5], slope=ancovaBeta, colour="#00BA38", linetype="dotdash", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[6], slope=ancovaBeta, colour="#00BA38", linetype="dotted", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[7], slope=ancovaBeta, colour="#619CFF", linetype="dashed", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[8], slope=ancovaBeta, colour="#619CFF", linetype="dotdash", size=1) +
  geom_abline(intercept=ancovaIntercept+ancovaAlphas[9], slope=ancovaBeta, colour="#619CFF", linetype="dotted", size=1) +
  labs(x = expression(paste("Acceleration (m/s" ^2, ")")), y = "Response Time (ms)") +
  scale_x_continuous(breaks=c(18,36,54,72,90), label=c("18","36","54","72","90"), limits=c(14,94)) +
  scale_y_continuous(breaks=seq(560, 740, by = 20), limits=c(560,740)) +
  scale_linetype_discrete(name="Vi",labels=c("3","9","15")) +
  guides(shape=FALSE,colour=FALSE) +
  theme_classic()

ancovaAdjMean <- as.numeric(effect(term="ancovaCondition",ancovaLM)$fit) # keep estimate of SEM?
ancovaAdjSd <- as.numeric(sqrt(tapply(grandSd^2,gl(9,3),mean))) # to get average of SDs
ancovaAdjSe <- ancovaAdjSd/sqrt(35*3) # 35 subjects * 3 conditions in each point
ancovaSoa <- as.factor(c(rep("vestibular",3),rep("sync",3),rep("visual",3)))
ancovaSoa <- factor(ancovaSoa,levels=c("vestibular","sync","visual"))
ancovaVi <- as.factor(rep(c("3","9","15"),3))
ancovaVi <- factor(ancovaVi,levels=c(3,9,15)) # need to manually order because for some reason doesn't do this on its own
ancovaReact3 <- data.frame(ancovaAdjMean,ancovaAdjSe,ancovaSoa,ancovaVi)

ancovaPanelB <- ggplot(ancovaReact3, aes(x=ancovaVi,y=ancovaAdjMean,group=ancovaSoa,colour=ancovaSoa)) +
  geom_line(position=position_dodge(width=0.1), size=1) + 
  geom_point(position=position_dodge(width=0.1), aes(shape=ancovaSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.1), aes(ymin=ancovaAdjMean-ancovaAdjSe, ymax=ancovaAdjMean+ancovaAdjSe), width=0.1, size=1) +
  labs(x = "Initial Velocity (m/s)", y = "Response Time (ms)") + 
  scale_y_continuous(breaks=seq(560, 740, by = 20), limits=c(560,740)) +
  scale_colour_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  scale_shape_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  theme_classic()

## This is Figure 3 (save manually as .eps to avoid pixellation)
ggarrange(ancovaPanelA, ancovaPanelB, labels = c("A", "B"),
          common.legend = TRUE, legend="right", nrow = 1, ncol = 2
)

## ---------------------------------------------------------- ##

## Stats for Figure 4: SOA x Vi x Block (repeated measures ANOVA)
blockSoa <- as.factor(c(rep("vestibular",12),rep("sync",12),rep("visual",12)))
blockVi <- as.factor(rep(c(rep("3",4),rep("9",4),rep("15",4)),3))
blockVi <- factor(blockVi,levels=c(3,9,15)) # need to manually order because for some reason doesn't do this on its own
blockBlock <- as.factor(rep(c("B1","B2","B3","B4"),3*3))
blockBlock <- factor(blockBlock,levels=c("B1","B2","B3","B4")) # need to manually order because for some reason doesn't do this on its own
blockFactors <- interaction(blockSoa,blockVi,blockBlock)
blockMatrixDim <- c(35,prod(dim(subjectBlockMean)[2:3]))
blockSubjectBlockMean <- matrix(subjectBlockMean,blockMatrixDim) # convert RTs from 3D to 2D format
blockSubjectInitBlockMean <- t(apply(blockSubjectBlockMean,1,tapply,gl(36,3),mean)) # average over Vf for each Vi
blockColumns <- c(1L,1+9L,1+18L,1+27L,2L,2+9L,2+18L,2+27L,3L,3+9L,3+18L,3+27L,4L,4+9L,4+18L,4+27L,5L,5+9L,5+18L,5+27L,6L,6+9L,6+18L,6+27L,7L,7+9L,7+18L,7+27L,8L,8+9L,8+18L,8+27L,9L,9+9L,9+18L,9+27L)
blockSubjectMean <- blockSubjectInitBlockMean[,blockColumns] # re-order
colnames(blockSubjectMean) <- blockFactors
blockDataframe <- data.frame(blockSoa,blockVi,blockBlock)
blockDesign <- ~blockSoa*blockVi*blockBlock
blockMlm <- lm(blockSubjectMean~1)
blockAov <- Anova(blockMlm,idata=blockDataframe,idesign=blockDesign,type="III")
summary(blockAov,multivariate=FALSE)

blockMean <- apply(blockSubjectMean,2,mean,na.rm=1)
blockSe <- apply(blockSubjectMean,2,sd,na.rm=1)/sqrt(35)
blockBlock<- as.factor(c(rep(c("B1","B2","B3","B4"),9)))
blockSoa <- as.factor(c(rep("vestibular",12),rep("sync",12),rep("visual",12)))
blockSoa <- factor(blockSoa,levels=c("vestibular","sync","visual"))
blockVi <- as.factor(rep(c(rep("3",4),rep("9",4),rep("15",4)),3))
blockVi <- factor(blockVi,levels=c(3,9,15)) # need to manually order because for some reason doesn't do this on its own

blockReact3 <- data.frame(blockMean,blockSe,blockBlock,blockSoa,blockVi)

## Figure 4: SOA x Vi x Block
blockVi3 <- ggplot(blockReact3 %>% filter(blockVi==3),aes(x=blockBlock,y=blockMean,group=blockSoa,colour=blockSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=blockSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=blockMean-blockSe, ymax=blockMean+blockSe), width=0.1, size=1) +
  labs(x = "Block", y = "Response Time (ms)") + 
  scale_y_continuous(breaks = seq(540, 760, by = 20), limits=c(540,760)) +
  guides(shape=FALSE,colour=FALSE) +
  annotate("text",x=2.5,y=760,label="Initial Velocity: 3m/s", size=5) +
  theme_classic()

blockVi9 <- ggplot(blockReact3 %>% filter(blockVi==9),aes(x=blockBlock,y=blockMean,group=blockSoa,colour=blockSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=blockSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=blockMean-blockSe, ymax=blockMean+blockSe), width=0.1, size=1) +
  labs(x = "Block", y = "Response Time (ms)") + 
  scale_y_continuous(breaks = seq(540, 760, by = 20), limits=c(540,760)) +
  guides(shape=FALSE,colour=FALSE) +
  annotate("text",x=2.5,y=760,label="Initial Velocity: 9m/s", size=5) +
  theme_classic()

blockVi15 <- ggplot(blockReact3 %>% filter(blockVi==15),aes(x=blockBlock,y=blockMean,group=blockSoa,colour=blockSoa)) +
  geom_line(position=position_dodge(width=0.15), size=1) +
  geom_point(position=position_dodge(width=0.15), aes(shape=blockSoa), size=4) +
  geom_errorbar(position=position_dodge(width=0.15), aes(ymin=blockMean-blockSe, ymax=blockMean+blockSe), width=0.1, size=1) +
  labs(x = "Block", y = "Response Time (ms)") + 
  scale_y_continuous(breaks = seq(540, 760, by = 20), limits=c(540,760)) +
  scale_colour_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  scale_shape_discrete(name="SOA",labels=c("Vestibular First (by 100ms)","In Sync","Visual First (by 100ms)")) +
  annotate("text",x=2.5,y=760, label="Initial Velocity: 15m/s", size=5) +
  theme_classic()

## This is Figure 4 (save manually as .eps to avoid pixellation)
ggarrange(blockVi3, blockVi9, blockVi15, labels = c("A", "B", "C"),
          common.legend = TRUE, legend="top", nrow = 1, ncol = 3
)

