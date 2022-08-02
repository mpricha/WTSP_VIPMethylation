
####Methylation data analysis
####Prichard, M.R.
####Last updated: 20220802

### General goals: 
### Differences in methylation
### Relationships between methylation and expression
### Relationships between methylation and probability of a difference

### Libraries used in this script ###
#stats and data management
library(dplyr)
library(plyr)
library(tidyr) 
library(e1071) #skewness
library(nlme) #lme
library(effectsize) #eta_squared
library(lsr) #cohensD
library(boot) #bootstrapping

#figures and data viz
library(ggplot2)
library(ggpubr)
library(ggtext) #wrap labels
library(gridExtra) #plot together
library(cowplot) #plot together
library(ggbreak) #axis breaks
library(ggnewscale) #heatmap with different color scales

#Colors key: tan blue red
#"#CDA434","#2142D7" ,"#CA3500"
here::i_am("r_script_2022_Github.R")
root <- here::here()
results <- paste0(root,"/results")

####  Methylation Data Frames and Data Cleaning         ###########

##  Here you would load in the bismark output files
##  The .csv read in just skips the step of reading each of those and stitching together the count data with 
##individual information for each bird.

###Load in data
setwd(root)
#use the .cov.gz files that Bismark outputs from the "methylation extraction" function
groups.bs <- read.table("files_manifest.txt", header=TRUE) #read in your sample manifest
#head(groups) #sanity check
#ddata <- c()
#for (x in list.files(pattern = "cov.gz$")) {
#  u<-read.table(x)
#  u$Label = factor(x)
#  ddata <- rbind(ddata, u)
#  cat(x, "\n ")
#}
#colnames(ddata)[c(1,2,3,4,5,6,7)]<-c("chromosome","position.1","position.2", "per.meth", 
#                                     "count.c", "count.t", "file")
#ddata2 <- merge(ddata, groups.bs, by="file")
#rm(u)
#rm(ddata)
#optional place to save all of your progress
#write.csv(ddata2, file = "methylation-count_data.csv") 

ddata2 <- read.csv("methylation-count_data.csv")

### Percent methylation and unshared CpGs
#The next steps make sure that CA's (that were artificially called CpGs previously) 
#are made to be 0% methylation.
#add total C's and T's as a measurement of coverage per CpG
ddata2$coverage <- ddata2$count.c + ddata2$count.t
ddata2$per.meth <- (ddata2$count.c/ddata2$coverage) #this is % methylation
#consider only data in the relevant sequence, remove CpGs not in your promoter 
#(would also work for amplicons). 
#Technically redundant with the "#select only CPGs in your CpG file" below but 
#I like doing both processing steps.
#Also filter out those with coverage less than desired amount (e.g. 10 for resolution of ~10%)
#For VIP the promoter range is 2257484 - 2259484
ddata2a <- subset(ddata2, position.1 < 2260000 & position.1 > 2258000) #all amplicons
ddata2b <- subset(ddata2a, coverage > 10)
unique(ddata2b$position.1)
###select only CpGs in your CpG file###
#cpg manifest based on the list used during the library preparation steps
cpgs <- read.table("cpg_manifest.txt", header=TRUE)
unique(cpgs$nw.pos1)
#head(cpgs)
dataCG1 <- ddata2b[ddata2b$position.1 %in% cpgs$NW.pos,] #match positions between the manifest and the library
cpgs$position.1 <- as.numeric(cpgs$NW.pos) #make the position values numeric
dataCG1 <-merge(ddata2b, cpgs, by="position.1") #bind with the rest of the information about the samples
dataCG1$position.1 <- as.factor(dataCG1$position.1) #change format of the position value
rm(ddata2, ddata2a, ddata2b) #remove unneeded data frames

###Set unshared CpG sites to = 0% methylation###
#Get rid of reads that are mis-assigned
#For me, these were called "unk" if they were assigned to "genome two" which was 
#supposed to be ZAL2m. There weren't many of these.
#I'm also removing the reads from the step before they were split to the two alleles. This data was probably
#slightly more robust (because the reads didn't have to have a SNP in them) but it will be simpler to take the averages.
dataCG1 <- subset(dataCG1, allele != "unk") 
dataCG1 <- subset(dataCG1, allele != "het" & allele != "homo")
#Check if your CpG sites are fitting to your data correctly, i.e. are all of these populated?
#DF with only shared CGs
shared1 <- subset(dataCG1, note == "shared")
#DFs for unshared CpGs only
unshared2 <- subset(dataCG1, note == "ZAL2_unshared" & allele == "ZAL2")
unshared2m <- subset(dataCG1, note == "ZAL2m_unshared" & allele == "ZAL2m")

###Set the polymorphic CpG sites on the wrong allele to zero 
unshared22m <- subset(dataCG1, note == "ZAL2_unshared" & allele == "ZAL2m")
unshared22m$per.meth <- 0
unshared2m2 <- subset(dataCG1, note == "ZAL2m_unshared" & allele == "ZAL2")
unshared2m2$per.meth <- 0 
unshared2m222m <- bind_rows(unshared22m,unshared2m2)
unshared <- bind_rows(unshared2m222m, unshared2)
unshared <- bind_rows(unshared, unshared2m)
dataCG <- bind_rows(unshared, shared1)
rm(unshared2, unshared2m, unshared22m, unshared2m2, unshared2m222m, unshared, shared1)
cpg.t<-unique(dataCG$position.1) #make a list of the cpgs in this data set to check that it's as expected
dataCG$position.a <- as.numeric(as.character(dataCG$position.1))
dataCG$coverage.a <- as.numeric(dataCG$coverage)
rm(dataCG1)
#stopping point. This is your final dataframe with all of the samples and % methylation appropriately calculated.



###Setting up a few data frames for statistical analysis

setwd(results)

#write.csv(dataCG,"completedata_bycpg.csv")

#Want averages by each allele within each individual and each individual
#first, nestlings were sequenced from two samples
#Combining brain regions will give a data frame of averages by each allele
#also want averages with and without polymorphic CpG sites
chicks <- subset(dataCG, age == "Chick") 
adults <- subset(dataCG, age == "Adult") #saving this to use to make a new whole df
chicksmean <- aggregate(per.meth ~ position.1+bird+allele, chicks, mean) #Calculate the new value
chicksa <- merge(chicks, chicksmean, by = c("position.1", "bird", "allele"), all = TRUE) #stitch it together
chicksa$per.meth <- chicksa$per.meth.y #rename variable
chicksa$per.meth.x <- NULL #remove variables
chicksa$per.meth.y <- NULL
chicksb <- ddply(chicksa, c("position.1","bird","allele"), head, 1)
dfcg <- rbind(chicksb, adults) #stitch the new data frame
rm(chicksb, chicksa, chicksmean, chicks, adults) #remove unneeded data frames
#now to average by individual
dmeansa <- aggregate(per.meth ~ bird+allele, dfcg, mean) #calculate the new value
dmeans <- merge(groups.bs, dmeansa, by=c("bird","allele")) #connect with the rest of the data set
dmeans <- subset(dmeans, region != "HYPv") #get rid of duplicate data from averaging the brain regions in the chicks
#I'm doing this a second time so I can have another dataframe with the unshared sites excluded.
#This is the shared only data frame.
shared <- subset(dataCG, note == "shared")
chicks_shared <- subset(shared, age == "Chick")
adults_shared <- subset(shared, age == "Adult")
chicksmean <- aggregate(per.meth ~ position.1+bird+allele, chicks_shared, mean)
chicksa <- merge(chicks_shared, chicksmean, by = c("position.1", "bird", "allele"), all = TRUE)
chicksa$per.meth <- chicksa$per.meth.y
chicksa$per.meth.x <- NULL
chicksa$per.meth.y <- NULL
chicksb <- ddply(chicksa, c("position.1","bird","allele"), head, 1)
dfcg_shared <- rbind(chicksb, adults_shared)
rm(chicksb, chicksa, chicksmean, chicks_shared, adults_shared)#cleanup
dmeans_shareda <- aggregate(per.meth ~ bird+allele, dfcg_shared, mean)
dmeans_shared <- merge(dmeans_shareda, groups.bs, by=c("bird","allele"))
dmeans_shared <- subset(dmeans_shared, region != "HYPv")        
#write a data frame that has both shared and all average methylation
dmeans_shared$per.meth.s <- dmeans_shared$per.meth
dmeans_shared$per.meth <- NULL
means <- merge(dmeans, dmeans_shared)

rm(dmeans_shared,dmeans_shareda,dmeansa,dmeans,dfcg_shared,shared)#cleanup

#Taking that data frame to be long form over shared and all means
ts <- subset(means,morph == "Tan")
ts2 <- pivot_longer(data = ts,
                    cols = c("per.meth","per.meth.s"),
                    names_to = "CpGs",
                    values_to = "per.meth")
ws <- subset(means, morph == "White")
ws2 <- pivot_longer(data = ws,
                    cols = c("per.meth","per.meth.s"),
                    names_to = "CpGs",
                    values_to = "per.meth")
means_allele <- rbind(ws2, ts2)
means_wide <- means

#write.csv(means_wide, "MeanMethyaltion_byallele_wide.csv")
#write.csv(means_allele, "MeanMethyaltion_byallele.csv")

rm(means)

#This section is for the "whole" data frame which includes the averages of the two alleles in WS birds
ws3 <- subset(means_wide, morph == "White")
wsmean <- aggregate(per.meth ~ bird, ws3, mean)
wsmean$CpGs <- "all"
wsmeans <- aggregate(per.meth.s ~ bird, ws3, mean)
wsmeans$CpGs <- "shared"
wsmeans$per.meth <- wsmeans$per.meth.s
wsmeans$per.meth.s <- NULL
ws4 <- rbind(wsmean, wsmeans)
ws5 <- subset(ws3, select = -per.meth)
ws6 <- merge(ws4, ws5, by = "bird")
ws6$per.meth.s <- NULL
ws6b <- ddply(ws6, c("bird","CpGs"), head, 1) #combine duplicates
means_whole <- rbind(ws6b, ts2) #add to ts df which only had one allele

#write.csv(means_whole,"MeanMethylation_bybird.csv")

rm(ts2, ts, ws, ws2, ws3, ws4,ws5,ws6,ws6b,wsmean,wsmeans)#cleanup

###Working Data frames
dfcg$distfromTSSn <- as.numeric(dfcg$NW.rel.pos)
dfcg <- subset(dfcg, bird != "WSF13-04") #this bird seems to have been mis-labeled as a WS, omitting
#log transformation
#adding a small constant to all of the values so there are no zeros 
#for the log transformation, constant is half of smallest value in df
dfcg$per.meth1 <- log(dfcg$per.meth+min(dfcg$per.meth[dfcg$per.meth > 0])/2) 
hist(dfcg$per.meth1,breaks=100)
skewness(dfcg$per.meth1)

#relevant and useful subset datasets
adult.dfcg <- subset(dfcg, age == "Adult")

chick.dfcg <- subset(dfcg, age == "Chick")

ws.adult <- subset(adult.dfcg, morph == "White")
ws.chick <- subset(chick.dfcg, morph == "White")
zal2.adult <- subset(adult.dfcg, allele == "ZAL2")
zal2.chick <- subset(chick.dfcg, allele == "ZAL2")

allele.adult <- subset(means_allele, age == "Adult")
allele.chick <- subset(means_allele, age == "Chick")

whole.adult <- subset(means_whole, age == "Adult")
whole.chick <- subset(means_whole, age == "Chick")

dists <- subset(dfcg, select = c('position.1','NW.rel.pos','nest'))
dists <- ddply(dists, "position.1",head,1)

#removing a bird with ambiguous notes on morph
allele.adult <- subset(allele.adult, bird != "WSF13-04")
whole.adult <- subset(whole.adult, bird != "WSF13-04")

#useful lists
#as a reminder, cpg.t is the whole list of CpGs
cpg.2<-c("2258417","2258431","2258616","2258795","2259319") #these are all of the ZAL2m CpGs, or the CpGs that aren't in the ZAL2 allele
cpg.2m<-c("2259224","2258582","2259007") #these are all of the ZAL2 CpGs, or the CpGs that aren't in the ZAL2m allele
cpg.poly<-c("2258417","2258431","2258616","2258795","2259319","2259224","2258582","2259007") #all polymorphic
cpg.tt<-setdiff(cpg.t, cpg.2) #CpG list without the ZAL2m specific CpGs
cpg.ttm <- setdiff(cpg.t, cpg.2m) #CpG list without the ZAL2 specific CpGs
cpg.shared <- setdiff(cpg.t, cpg.poly) #CpG list without the polymorphic CpGs

adult.dfcg_bymorph <- split(adult.dfcg, adult.dfcg$morph)
chick.dfcg_bymorph <- split(chick.dfcg, chick.dfcg$morph)
morph_list <- c("White","Tan")
adult.dfcg_by2.2m <- split(adult.dfcg, adult.dfcg$allele)
chick.dfcg_by2.2m <- split(chick.dfcg, chick.dfcg$allele)
by2.2m_list <- c("ZAL2","ZAL2m")
adult.dfcg_byallele <- split(adult.dfcg, adult.dfcg$allele.morph)
chick.dfcg_byallele <- split(chick.dfcg, chick.dfcg$allele.morph)
alleles_list <- c("ZAL2-TS","ZAL2-WS","ZAL2m-WS")


####  RNA-seq Data Frames and Cleaning                  ################

setwd(root)
#normsallele <- read.table("HypNormalizationallele.txt", header = TRUE)
#groups <- read.table("files_rna.txt", header=TRUE) #read in your sample manifest

### Data cleaning ##
#ddata <- c()
#for (x in list.files(pattern = ".count$")) {
#  u<-read.table(x, nrows = 16155)
#  u[u==""] <- NA
#  u <- na.omit(u)
#  u$Label = factor(x)
#  ddata <- rbind(ddata, u)
#  cat(x, "\n ")
#}
#colnames(ddata)[c(1,2,3,4)]<-c("gene_id","gene","count", "file")
#ddata2 <- merge(ddata, groups, by="file")
#Just VIP
#VIP <- subset(ddata2, ddata2$gene == "VIP")
#Just HYP
#VIP.hyp1 <- subset(VIP, region !="MeA")
#Get rid of ZAL2m from TS
#VIP.hypTSZAL2 <- subset(VIP.hyp1, VIP.hyp1$Morph == "TS") 
#VIP.hypTSZAL2 <- subset(VIP.hypTSZAL2, VIP.hypTSZAL2$allele == "ZAL2") 
#VIP.hypWSZAL2m <- subset(VIP.hyp1, Morph == "WS")
#VIP.hypWSZAL2m <- subset(VIP.hypWSZAL2m, allele == "ZAL2m")
#VIP.hypWSZAL2 <- subset(VIP.hyp1, Morph == "WS")
#VIP.hypWSZAL2 <- subset(VIP.hypWSZAL2, allele == "ZAL2")
#VIP.hyp <- bind_rows(VIP.hypTSZAL2, VIP.hypWSZAL2)
#VIP.hyp <- bind_rows(VIP.hyp, VIP.hypWSZAL2m)
#normalize
#VIP.hyp <- merge(VIP.hyp, normsallele, by="file")
#VIP.hyp$norm.count <- VIP.hyp$count/VIP.hyp$Hypnormalization
#write.csv(VIP.hyp, "RNAseq_VIPHypAllele.csv")

#rm(VIP, VIP.hyp1, VIP.hypTSZAL2, VIP.hypWSZAL2, VIP.hypWSZAL2m, ddata2,groups, u, ddata, normsallele)

VIP.hyp <- read.csv("RNAseq_VIPHypAllele.csv")

###   working data frames   ###
brcorr.means <- merge(VIP.hyp, means_allele, by = c("bird","allele"))
brcorr.means.a <- subset(brcorr.means, age == "Adult")
brcorr.means.c <- subset(brcorr.means, age == "Chick")
dfcg$position.1a <- as.numeric(as.character(dfcg$position.1))
brcorr.dfcg <- merge(dfcg, VIP.hyp, by = c("bird","allele"), all.x = TRUE)
brcorr.dfcg <- subset(brcorr.dfcg, !is.na(brcorr.dfcg$norm.count)) #there are 5 birds in the bs data that aren't in the rna ddata
brcorr.dfcg.a <- subset(brcorr.dfcg, age == "Adult")
brcorr.dfcg.c <- subset(brcorr.dfcg, age == "Chick")

brcorr.dfcg.a_byallele <- split(brcorr.dfcg.a, brcorr.dfcg.a$allele.morph, drop = TRUE)
brcorr.dfcg.c_byallele <- split(brcorr.dfcg.c, brcorr.dfcg.c$allele.morph, drop = TRUE)

setwd(results)

#write.csv(brcorr.dfcg.a, file ="adult-dfcg.csv")
#write.csv(brcorr.dfcg.c, file = "chick-dfcg.csv")



########    some summary statistics                       ######
#print means
mean(chick.dfcg$per.meth)
mean(whole.chick$per.meth)
mean(allele.chick$per.meth)
mean(adult.dfcg$per.meth)
mean(whole.adult$per.meth)
mean(allele.adult$per.meth)

#coverage after data cleaning
coverage <- function (x) {
  c(
    mean(x$coverage),
    sd(x$coverage),
    IQR(x$coverage)
  )
}

c<-do.call("rbind",lapply(split(dfcg,list(dfcg$morph,dfcg$age)),coverage))
colnames(c) <- c("mean","sd","IQR")
c <- as.data.frame(c)
stats_coverage.sum = tibble::rownames_to_column(c, "morph.age")

#adults
#setting up the data frame
adult.dfcg_byallele <- split(adult.dfcg, adult.dfcg$allele.morph)
chick.dfcg_byallele <- split(chick.dfcg, chick.dfcg$allele.morph)
alleles_list <- c("ZAL2-TS","ZAL2-WS","ZAL2m-WS")

summarize_dfcg <- function (x) {
  c(
  nrow(x),
  mean(x$per.meth),
  median(x$per.meth),
  sd(x$per.meth),
  mean(x$per.meth1),
  median(x$per.meth1),
  sd(x$per.meth1)
  )
}

stats_adult.summaries <- vector(mode = "list")
stats_chick.summaries <- vector(mode = "list")

for (i in alleles_list) {
  b <- adult.dfcg_byallele[[i]]
  c <- do.call("rbind",lapply(split(b,b$position.1),summarize_dfcg))
  colnames(c) <- c("n","mean_per.meth","median_per.meth","sd_per.meth","mean_per.meth1","median_per.meth1","sd_per.meth1")
  c <- as.data.frame(c)
  c = tibble::rownames_to_column(c, "position.1")
  stats_adult.summaries[[i]] <- as.data.frame(c)
  b <- chick.dfcg_byallele[[i]]
  c <- do.call("rbind",lapply(split(b,b$position.1),summarize_dfcg))
  colnames(c) <- c("n","mean_per.meth","median_per.meth","sd_per.meth","mean_per.meth1","median_per.meth1","sd_per.meth1")
  c <- as.data.frame(c)
  c = tibble::rownames_to_column(c, "position.1")
  stats_chick.summaries[[i]] <- as.data.frame(c)
  }

write.csv(stats_adult.summaries, "adult_summarystats.csv")
write.csv(stats_chick.summaries, "nestlings_summarystats.csv")

mean(adult.dfcg_byallele[["ZAL2-TS"]]$per.meth)
sd(adult.dfcg_byallele[["ZAL2-TS"]]$per.meth)
mean(adult.dfcg_byallele[["ZAL2-WS"]]$per.meth)
sd(adult.dfcg_byallele[["ZAL2-WS"]]$per.meth)
mean(adult.dfcg_byallele[["ZAL2m-WS"]]$per.meth)
sd(adult.dfcg_byallele[["ZAL2m-WS"]]$per.meth)

mean(chick.dfcg_byallele[["ZAL2-TS"]]$per.meth)
sd(chick.dfcg_byallele[["ZAL2-TS"]]$per.meth)
mean(chick.dfcg_byallele[["ZAL2-WS"]]$per.meth)
sd(chick.dfcg_byallele[["ZAL2-WS"]]$per.meth)
mean(chick.dfcg_byallele[["ZAL2m-WS"]]$per.meth)
sd(chick.dfcg_byallele[["ZAL2m-WS"]]$per.meth)

rm(b,c)


######    Adult: Level of allele                        ########################

#ANOVAs
stats_allelea.aov.sex <- lme(log(per.meth) ~ allele.morph + sex,
                             random = ~1|bird,
                             data = subset(allele.adult, CpGs == "per.meth.s"))
anova(stats_allelea.aov.sex)
eta_squared(stats_allelea.aov.sex)[2,2]

stats_allelea.aov.all <- lme(log(per.meth) ~ allele.morph,
                   random = ~1|bird,
                   data = subset(allele.adult, CpGs == "per.meth"))
anova(stats_allelea.aov.all)
eta_squared(stats_allelea.aov.all)[1,2]

stats_allelea.aov.shared <- lme(log(per.meth) ~ allele.morph,
                   random = ~1|bird,
                   data = subset(allele.adult, CpGs == "per.meth.s"))
anova(stats_allelea.aov.shared)
eta_squared(stats_allelea.aov.shared)[1,2]

#POSTHOCS
g <- subset(allele.adult, CpGs == "per.meth") %>%
  filter(allele.morph == "ZAL2-WS" | allele.morph == "ZAL2-TS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2-WS", "ZAL2-TS")))

hist(log(g$per.meth))
skewness(log(g$per.meth))
b <- t.test(log(per.meth) ~ allele.morph, data = g)
cohensD(log(per.meth) ~ allele.morph,
        data = g)

g <- subset(allele.adult, CpGs == "per.meth") %>%
  filter(allele.morph == "ZAL2m-WS" | allele.morph == "ZAL2-WS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2m-WS", "ZAL2-WS")))

hist(log(g$per.meth))
skewness(log(g$per.meth))
t.test(log(per.meth) ~ allele.morph,
       paired = TRUE,
       data = g)
cohensD(log(per.meth) ~ allele.morph,
        data = g)


g <- subset(allele.adult, CpGs == "per.meth.s") %>%
  filter(allele.morph == "ZAL2-TS" | allele.morph == "ZAL2-WS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2-TS", "ZAL2-WS")))

t.test(log(per.meth) ~ allele.morph,
       data = g)
cohensD(log(per.meth)~allele.morph,
        data = g)

g <- subset(allele.adult, CpGs == "per.meth.s") %>%
  filter(allele.morph == "ZAL2m-WS" | allele.morph == "ZAL2-WS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2m-WS", "ZAL2-WS")))

t.test(log(per.meth) ~ allele.morph,
       data = g,
       paired = TRUE)
cohensD(log(per.meth)~allele.morph,
        data = g)

rm(g)


######    chick: Level of allele                        ###########

#ANOVAs
stats_allelec.aov.sex <- lme(log(per.meth) ~ allele.morph + sex,
                             random = ~1|nest/bird,
                             data = subset(allele.chick, CpGs == "per.meth.s"))
anova(stats_allelec.aov.sex)
eta_squared(stats_allelec.aov.sex)[2,2]

stats_allelec.aov.all <- lme(log(per.meth) ~ allele.morph,
                             random = ~1|nest/bird,
                             data = subset(allele.chick, CpGs == "per.meth"))
anova(stats_allelec.aov.all)
eta_squared(stats_allelec.aov.all)[1,2]

stats_allelec.aov.shared <- lme(log(per.meth) ~ allele.morph,
                   random = ~1|nest/bird,
                   data = subset(allele.chick, CpGs == "per.meth.s"))
anova(stats_allelec.aov.shared)
eta_squared(stats_allelec.aov.shared)[1,2]


#POST HOC
g <- subset(allele.chick, CpGs == "per.meth") %>%
  filter(allele.morph == "ZAL2-WS" | allele.morph == "ZAL2-TS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2-WS", "ZAL2-TS")))

t.test(log(per.meth) ~ allele.morph,
       data = g)
cohensD(log(per.meth) ~ allele.morph,
        data = g)

g <- subset(allele.chick, CpGs == "per.meth") %>%
  filter(allele.morph == "ZAL2m-WS" | allele.morph == "ZAL2-WS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2m-WS", "ZAL2-WS")))

t.test(log(per.meth) ~ allele.morph,
       paired = TRUE,
       data = g)
cohensD(log(per.meth) ~ allele.morph,
        data = g)

g <- subset(allele.chick, CpGs == "per.meth.s") %>%
  filter(allele.morph == "ZAL2-WS" | allele.morph == "ZAL2-TS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2-WS", "ZAL2-TS")))


t.test(log(per.meth) ~ allele.morph,
       data = g)
cohensD(log(per.meth) ~ allele.morph,
        data = g)

g <- subset(allele.chick, CpGs == "per.meth.s") %>%
  filter(allele.morph == "ZAL2m-WS" | allele.morph == "ZAL2-WS") %>%
  mutate(allele.morph = factor(allele.morph, levels = c("ZAL2m-WS", "ZAL2-WS")))

t.test(log(per.meth) ~ allele.morph,
       data = g,
       paired = TRUE)
cohensD(log(per.meth) ~ allele.morph,
        data = g)

rm(g)

######    Individual CpGs                               ################

#Methylation across the promoter in order from Distance from TSS 
#Iterative ANOVAs to identify differentially methyalted CpGs
#adults

CpG_ANOVAs <- function (x) {
  a <- lme(per.meth1 ~ allele.morph, random = ~1|bird, data = x)
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[2,1],
    anova(a)[2,2]
  )
}

b <- do.call("rbind",lapply(split(adult.dfcg,adult.dfcg$position.1),CpG_ANOVAs))
colnames(b) = c("F-statistic","p-value","eta","n","df num","df den")
b <- as.data.frame(b)
stats_CpGs.ANOVAs.a = tibble::rownames_to_column(b, "Position")
stats_CpGs.ANOVAs.a.2 <- stats_CpGs.ANOVAs.a[stats_CpGs.ANOVAs.a$Position %in% cpg.shared,]
write.csv(stats_CpGs.ANOVAs.a.2, "stats_cpg-anovas-a.csv")

#cpgs that were p < 0.05
cpgs_posthoc.a <- stats_CpGs.ANOVAs.a.2[stats_CpGs.ANOVAs.a.2$`p-value` < 0.05, ]$Position

#nestlings
CpG_ANOVAs2 <- function (x) {
  a <- lme(per.meth1 ~ allele.morph, random = ~1|bird/nest, data = x)
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[2,1],
    anova(a)[2,2]
  )
}

b <- do.call("rbind",lapply(split(chick.dfcg,chick.dfcg$position.1),CpG_ANOVAs2))
colnames(b) = c("F-statistic","p-value","eta","n","df num","df den")
b <- as.data.frame(b)
stats_CpGs.ANOVAs.c = tibble::rownames_to_column(b, "Position")
stats_CpGs.ANOVAs.c.2 <- stats_CpGs.ANOVAs.c[stats_CpGs.ANOVAs.c$Position %in% cpg.shared,]
write.csv(stats_CpGs.ANOVAs.c.2, "stats_cpg-anovas-c.csv")

#cpgs that were p < 0.05
cpgs_posthoc.c <- stats_CpGs.ANOVAs.c.2[stats_CpGs.ANOVAs.c.2$`p-value` < 0.05, ]$Position

rm(b)


###   POSTHOC   ###

#adults

#post hoc test of allele within ws birds
CpG_ANOVAs <- function (x) {
  a <- lme(per.meth1 ~ allele, random = ~1|bird, data = x)
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[2,1],
    anova(a)[2,2]
  )
} 
b <- do.call("rbind",lapply(split(adult.dfcg_bymorph[["White"]],adult.dfcg_bymorph[["White"]]$position.1),CpG_ANOVAs))
colnames(b) = c("F-statistic","p-value","eta","n","df num","df den")
b <- as.data.frame(b)
stats_CpGs.posthoc.wsa = tibble::rownames_to_column(b, "Position")
stats_CpGs.posthoc.wsa2 <- stats_CpGs.posthoc.wsa[stats_CpGs.posthoc.wsa$Position %in% cpgs_posthoc.a, ]

#post hoc test of morph within zal2
CpG_ANOVAs3 <- function (x) {
  a <- lm(per.meth1 ~ morph, data = x)
  c(
    anova(a)[1,4],
    anova(a)[1,5],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[1,1],
    anova(a)[2,1]
  )
}

c <- do.call("rbind",lapply(split(adult.dfcg_by2.2m[["ZAL2"]],adult.dfcg_by2.2m[["ZAL2"]]$position.1),CpG_ANOVAs3))
colnames(c) = c("F.statistic","p.value","eta","n","df num","df den")
c <- as.data.frame(c)
stats_CpGs.posthoc.zal2a = tibble::rownames_to_column(c, "Position")
stats_CpGs.posthoc.zal2a2 <- stats_CpGs.posthoc.zal2a[stats_CpGs.posthoc.zal2a$Position %in% cpgs_posthoc.a, ]

write.csv(stats_CpGs.posthoc.wsa2, "stats_CpGs-posthoc-ws-a.csv")
write.csv(stats_CpGs.posthoc.zal2a2, "stats_CpGs-posthoc-zal2-a.csv")

#nestlings
#post hoc test of allele within ws birds
CpG_ANOVAs2 <- function (x) {
  a <- lme(per.meth1 ~ allele, random = ~1|bird, data = x)
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[2,1],
    anova(a)[2,2]
  )
} #this is the same as above

d <- do.call("rbind",lapply(split(chick.dfcg_bymorph[["White"]],chick.dfcg_bymorph[["White"]]$position.1),CpG_ANOVAs2))
colnames(d) = c("F-statistic","p-value","eta","n","df num","df den")
d <- as.data.frame(d)
stats_CpGs.posthoc.wsc = tibble::rownames_to_column(d, "Position")
stats_CpGs.posthoc.wsc2 <- stats_CpGs.posthoc.wsc[stats_CpGs.posthoc.wsc$Position %in% cpgs_posthoc.c, ]
#post hoc test of morph within zal2
CpG_ANOVAs4 <- function (x) {
  a <- lme(per.meth1 ~ morph, random = ~1|nest, data = x)
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    eta_squared(a)[1,2],
    nrow(x),
    anova(a)[2,1],
    anova(a)[2,2]
  )
}

l <- split(chick.dfcg_by2.2m[["ZAL2"]],chick.dfcg_by2.2m[["ZAL2"]]$position.1)
m <- l[names(l) %in% cpg.tt]

e <- do.call("rbind",lapply(m,CpG_ANOVAs4))
colnames(e) = c("F-statistic","p-value","eta","n","df num","df den")
e <- as.data.frame(e)
stats_CpGs.posthoc.zal2c = tibble::rownames_to_column(e, "Position")
stats_CpGs.posthoc.zal2c2 <- stats_CpGs.posthoc.zal2c[stats_CpGs.posthoc.zal2c$Position %in% cpgs_posthoc.c, ]

write.csv(stats_CpGs.posthoc.wsc2, "stats_CpGs-posthoc-ws-c.csv")
write.csv(stats_CpGs.posthoc.zal2c2, "stats_CpGs-posthoc-zal2-c.csv")

###   Comparing just the poly sites to zero  ###

WS.polyCpGs_t.tests <- function (x) {
  c(
    t.test(x$per.meth1, alternative = "less", )$statistic,
    t.test(x$per.meth1, alternative = "less")$p.value,
    cohensD(x$per.meth1),
    nrow(x)
  )
}

zed <- c("t.value","p.value","Cohens.d", "n")

a <- split(adult.dfcg_byallele[["ZAL2m-WS"]],adult.dfcg_byallele[["ZAL2m-WS"]]$position.1)
b <- a[names(a) %in% cpg.2]
c <- do.call("rbind",lapply(b,WS.polyCpGs_t.tests))
colnames(c) = zed 
stats_polyt.tests <- 
  c %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")

d <- split(adult.dfcg_byallele[["ZAL2-WS"]],adult.dfcg_byallele[["ZAL2-WS"]]$position.1)
e <- d[names(d) %in% cpg.2m]
f <- do.call("rbind",lapply(e,WS.polyCpGs_t.tests))
colnames(f) = zed 
g <- 
  f %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
stats_polyt.tests <- rbind(stats_polyt.tests, g)

### comparing ZAL2 specific CpGs between WS vs TS ###

ZAL2.polyCpGs_t.tests <- function (x) {
  c(
    t.test(per.meth1 ~ morph, data = x)$statistic,
    t.test(per.meth1 ~ morph, data = x)$p.value,
    cohensD(per.meth1 ~ morph, data = x),
    nrow(x)
  )
}

h <- split(adult.dfcg_by2.2m[["ZAL2"]],adult.dfcg_by2.2m[["ZAL2"]]$position.1)
i <- h[names(h) %in% cpg.2m]
j <- do.call("rbind",lapply(i,ZAL2.polyCpGs_t.tests))
colnames(j) = zed 
k <- 
  j %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
stats_polyt.tests$test <- "WS"
k$test <- "ZAL2"
stats_polyt.tests.a <- rbind(stats_polyt.tests, k)

###   again in nestlings   ###

a <- split(chick.dfcg_byallele[["ZAL2m-WS"]],chick.dfcg_byallele[["ZAL2m-WS"]]$position.1)
b <- a[names(a) %in% cpg.2]
c <- do.call("rbind",lapply(b,WS.polyCpGs_t.tests))
colnames(c) = zed 
stats_polyt.tests <- 
  c %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")

d <- split(chick.dfcg_byallele[["ZAL2-WS"]],chick.dfcg_byallele[["ZAL2-WS"]]$position.1)
e <- d[names(d) %in% cpg.2m]
f <- do.call("rbind",lapply(e,WS.polyCpGs_t.tests))
colnames(f) = zed 
g <- 
  f %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
stats_polyt.tests <- rbind(stats_polyt.tests, g)

h <- split(chick.dfcg_by2.2m[["ZAL2"]],chick.dfcg_by2.2m[["ZAL2"]]$position.1)
i <- h[names(h) %in% cpg.2m]
j <- do.call("rbind",lapply(i,ZAL2.polyCpGs_t.tests))
colnames(j) = zed 
k <- 
  j %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
stats_polyt.tests$test <- "WS"
k$test <- "ZAL2"
stats_polyt.tests.c <- rbind(stats_polyt.tests, k)

stats_polyt.tests.a$age = "adult"
stats_polyt.tests.c$age = "chick"
stats_polyt.tests <- rbind(stats_polyt.tests.a, stats_polyt.tests.c)
write.csv(stats_polyt.tests, file = "stats_poly-t-tests.csv")

rm(a,b,c,d,e,f,g,h,i,j,k,l,m, zed)

##########   Bootstrapping                      ###########

###   Is there more likely to be a pairwise difference when methyaltion is high? ###

#functions to pull relevant stats
diff.2 <- function(d1, i) {
  d = d1[i, ]
  Mean = tapply(X=d$per.meth1, INDEX = d$allele.morph, mean)
  Diff = Mean["ZAL2-TS"] - Mean["ZAL2-WS"]
  Diff
}

diff.ws <- function(d1, i) {
  d = d1[i, ]
  Mean = tapply(X=d$per.meth1, INDEX = d$allele.morph, mean)
  Diff = Mean["ZAL2-WS"] - Mean["ZAL2m-WS"]
  Diff
}


diff2.2 <- function(d1, i) {
  d = d1
  d$group <- d$allele.morph[i] #permutate "morph" as a group to randomly assign as the "null" condition
  Mean = tapply(X=d$per.meth1, INDEX = d$group, mean)
  Diff = Mean["ZAL2-TS"] - Mean["ZAL2-WS"]
  Diff
}

diff2.ws <- function(d1, i) {
  d = d1
  d$group <- d$allele.morph[i] #permutate "allele" as a group to randomly assign as the "null" condition
  Mean = tapply(X=d$per.meth1, INDEX = d$group, mean)
  Diff = Mean["ZAL2m-WS"] - Mean["ZAL2-WS"]
  Diff
}

#ZAL2 diffs
bootdiffs.2 <- function(x) { 
  results <- boot(data=x, statistic=diff.2, R=1000) #bootstrap differences
  results_null <- boot(data=x, statistic=diff2.2, R=1000) #bootstrap differences under NULL (morph scrambled)
  e <- quantile(results$t,c(0.025,0.975), na.rm = T)
  k <- table(cut(results_null$t,breaks=c(-Inf,-abs(results$t0),abs(results$t0)-0.000001,Inf)))/length(results_null$t)
  k <- as.data.frame(k)
  c(
    t.test(per.meth1 ~ allele.morph, x)$p.value,
    t.test(per.meth1~ allele.morph, x)$estimate,
    cohensD(per.meth1 ~ allele.morph, data = x),
    results$t0,
    quantile(results$t,c(0.025,0.975), na.rm= T), #95% CI for bootstrapped differences in means, if spans zero then α ≤ p-value, then do not reject H0
    ifelse(0 > e[1] & 0 < e[2], yes = "fail", no = "pass"), #pass = reject null, fail = fail to reject null
    quantile(results_null$t,c(0.025,0.975), na.rm = T), #95% CI for bootstrapped differences under "null" condition
    1-k[2,2] #tails of the "null" distribution beyond the observed difference; AKA probability of more extreme difference
  )
}

bootdiffs.ws <- function(x) { 
  results <- boot(data=x, statistic=diff.ws, R=1000) #bootstrap differences
  results_null <- boot(data=x, statistic=diff2.ws, R=1000) #bootstrap differences under NULL (morph scrambled)
  e<-quantile(results$t,c(0.025,0.975),na.rm = T)
  k <- table(cut(results_null$t,breaks=c(-Inf,-abs(results$t0),abs(results$t0)-0.000001,Inf)))/length(results_null$t)
  k <- as.data.frame(k)
  c(
    t.test(per.meth1 ~ allele.morph, x)$p.value,
    t.test(per.meth1~ allele.morph, x)$estimate,
    cohensD(per.meth1 ~ allele.morph, data = x),
    results$t0,
    quantile(results$t,c(0.025,0.975), na.rm=T),
    ifelse(0 > e[1] & 0 < e[2], yes = "fail", no = "pass"),
    quantile(results_null$t,c(0.025,0.975),na.rm=T),
    1-k[2,2]
  )
}

bootdiffs_columns <- c("t_pvalue","Alt_mean","ZAL2.WS_mean","cohensD",
                       "boot_t0","boot_0.025","boot_0.975","IQR_span0",
                       "H0_0.025","H0_0.975","prob_xtrme")

b <- split(adult.dfcg_by2.2m[["ZAL2"]],adult.dfcg_by2.2m[["ZAL2"]]$note)

c <- split(adult.dfcg_bymorph[["White"]],adult.dfcg_bymorph[["White"]]$note)

q <- do.call("rbind",lapply(split(b[["shared"]],b[["shared"]]$position.1,drop = TRUE), bootdiffs.2))
colnames(q) <- bootdiffs_columns
q <- q %>% as.data.frame() %>% tibble::rownames_to_column("position.1")
q$boot_test <- "TS2 vs WS2"

r <- do.call("rbind",lapply(split(c[["shared"]], c[["shared"]]$position.1, drop = TRUE), bootdiffs.ws))
colnames(r) <- bootdiffs_columns
r <- r %>% as.data.frame %>% tibble::rownames_to_column("position.1")
r$boot_test <- "WS2m vs WS2"

#long
bootdiffs.a2 <- rbind(r, q)
bootdiffs.a2$per.meth_mean.2 <- (as.numeric(q$Alt_mean) + as.numeric(q$ZAL2.WS_mean))/2
bootdiffs.a2$per.meth_mean.ws <- (as.numeric(r$Alt_mean) + as.numeric(r$ZAL2.WS_mean))/2
bootdiffs.a2$Alt_mean <- as.numeric(bootdiffs.a2$Alt_mean)
bootdiffs.a2$ZAL2.WS_mean <- as.numeric(bootdiffs.a2$ZAL2.WS_mean)
bootdiffs.a2$prob_xtrme <- as.numeric(bootdiffs.a2$prob_xtrme)
bootdiffs.a2$boot_t0 <- as.numeric(bootdiffs.a2$boot_t0)



#set up a wide format of the same data
bootdiffs.ws_columns <- c("position.1","t_pvalue.ws","Alt_mean.ws","ZAL2.WS_mean.ws","cohensD.ws",
                          "boot_t0.ws","boot_0.025.ws","boot_0.975.ws","IQR_span0.ws",
                          "H0_0.025.ws","H0_0.975.ws","prob_xtrme.ws","boot_test.ws")

colnames(r) <- bootdiffs.ws_columns
bootdiffs.a <- merge(r, q, by = "position.1")
write.csv(bootdiffs.a, file = "bootdiffs.a_shared.csv")

bootdiffs.a$per.meth_mean.2 <- (as.numeric(q$Alt_mean) + as.numeric(q$ZAL2.WS_mean))/2
bootdiffs.a$per.meth_mean.ws <- (as.numeric(r$Alt_mean) + as.numeric(r$ZAL2.WS_mean))/2
bootdiffs.a$prob_xtrme <- as.numeric(bootdiffs.a$prob_xtrme)
bootdiffs.a$prob_xtrme.ws <- as.numeric(bootdiffs.a$prob_xtrme.ws)


ggplot(bootdiffs.a, aes(x=prob_xtrme, y=per.meth_mean.2)) + 
  geom_point(show.legend = T, shape=16, size = 1 ) +
  geom_smooth(method = "lm", se = F, size = 1.5, show.legend = T) +
  theme_classic(base_size = 10)

bootdiffs.2.lm <- lm(per.meth_mean.2 ~ prob_xtrme, data = subset(bootdiffs.a2,boot_test == "TS2 vs WS2"))
summary(bootdiffs.2.lm)
anova(bootdiffs.2.lm)

ggplot(bootdiffs.a, aes(x=prob_xtrme.ws, y=per.meth_mean.ws)) + 
  geom_point(show.legend = T, shape=16, size = 1 ) +
  geom_smooth(method = "lm", se = F, size = 1.5, show.legend = T) +
  theme_classic(base_size = 10)

bootdiffs.ws.lm <- lm(per.meth_mean.ws ~ prob_xtrme, data = subset(bootdiffs.a2,boot_test == "WS2m vs WS2"))
summary(bootdiffs.ws.lm)
anova(bootdiffs.ws.lm)



### NESTLINGS ###
###   Is there more likely to be a pairwise difference when methyaltion is high?

#functions to pull relevant stats
diff.2 <- function(d1, i) {
  d = d1[i, ]
  Mean = tapply(X=d$per.meth1, INDEX = d$allele.morph, mean)
  Diff = Mean["ZAL2-TS"] - Mean["ZAL2-WS"]
  Diff
}

diff.ws <- function(d1, i) {
  d = d1[i, ]
  Mean = tapply(X=d$per.meth1, INDEX = d$allele.morph, mean)
  Diff = Mean["ZAL2-WS"] - Mean["ZAL2m-WS"]
  Diff
}


diff2.2 <- function(d1, i) {
  d = d1
  d$group <- d$allele.morph[i] #permutate "morph" as a group to randomly assign as the "null" condition
  Mean = tapply(X=d$per.meth1, INDEX = d$group, mean)
  Diff = Mean["ZAL2-TS"] - Mean["ZAL2-WS"]
  Diff
}

diff2.ws <- function(d1, i) {
  d = d1
  d$group <- d$allele.morph[i] #permutate "allele" as a group to randomly assign as the "null" condition
  Mean = tapply(X=d$per.meth1, INDEX = d$group, mean)
  Diff = Mean["ZAL2m-WS"] - Mean["ZAL2-WS"]
  Diff
}

#ZAL2 diffs
bootdiffs.2 <- function(x) { 
  results <- boot(data=x, statistic=diff.2, R=1000) #bootstrap differences
  results_null <- boot(data=x, statistic=diff2.2, R=1000) #bootstrap differences under NULL (morph scrambled)
  e <- quantile(results$t,c(0.025,0.975), na.rm = T)
  k <- table(cut(results_null$t,breaks=c(-Inf,-abs(results$t0),abs(results$t0)-0.00001,Inf)))/length(results_null$t)
  k <- as.data.frame(k)
  c(
    t.test(per.meth1 ~ allele.morph, x)$p.value,
    t.test(per.meth1~ allele.morph, x)$estimate,
    cohensD(per.meth1 ~ allele.morph, data = x),
    results$t0,
    quantile(results$t,c(0.025,0.975), na.rm= T), #95% CI for bootstrapped differences in means, if spans zero then α ≤ p-value, then do not reject H0
    ifelse(0 > e[1] & 0 < e[2], yes = "fail", no = "pass"), #pass = reject null, fail = fail to reject null
    quantile(results_null$t,c(0.025,0.975), na.rm = T), #95% CI for bootstrapped differences under "null" condition
    1-k[2,2] #tails of the "null" distribution beyond the observed difference; AKA probability of more extreme difference
  )
}

bootdiffs.ws <- function(x) { 
  results <- boot(data=x, statistic=diff.ws, R=1000) #bootstrap differences
  results_null <- boot(data=x, statistic=diff2.ws, R=1000) #bootstrap differences under NULL (morph scrambled)
  e<-quantile(results$t,c(0.025,0.975),na.rm = T)
  k <- table(cut(results_null$t,breaks=c(-Inf,-abs(results$t0),abs(results$t0)-0.00001,Inf)))/length(results_null$t)
  k <- as.data.frame(k)
  c(
    t.test(per.meth1 ~ allele.morph, x)$p.value,
    t.test(per.meth1~ allele.morph, x)$estimate,
    cohensD(per.meth1 ~ allele.morph, data = x),
    results$t0,
    quantile(results$t,c(0.025,0.975), na.rm=T),
    ifelse(0 > e[1] & 0 < e[2], yes = "fail", no = "pass"),
    quantile(results_null$t,c(0.025,0.975),na.rm=T),
    1-k[2,2]
  )
}

bootdiffs_columns <- c("t_pvalue","Alt_mean","ZAL2.WS_mean","cohensD",
                       "boot_t0","boot_0.025","boot_0.975","IQR_span0",
                       "H0_0.025","H0_0.975","prob_xtrme")

b <- split(chick.dfcg_by2.2m[["ZAL2"]],chick.dfcg_by2.2m[["ZAL2"]]$note)

c <- split(chick.dfcg_bymorph[["White"]],chick.dfcg_bymorph[["White"]]$note)


q <- do.call("rbind",lapply(split(b[["shared"]],b[["shared"]]$position.1,drop = TRUE), bootdiffs.2))
colnames(q) <- bootdiffs_columns
q <- q %>% as.data.frame() %>% tibble::rownames_to_column("position.1")
q$boot_test <- "TS2 vs WS2"

r <- do.call("rbind",lapply(split(c[["shared"]], c[["shared"]]$position.1, drop = TRUE), bootdiffs.ws))
colnames(r) <- bootdiffs_columns
r <- r %>% as.data.frame %>% tibble::rownames_to_column("position.1")
r$boot_test <- "WS2m vs WS2"

#long
bootdiffs.c2 <- rbind(r, q)
bootdiffs.c2$Alt_mean <- as.numeric(bootdiffs.c2$Alt_mean)
bootdiffs.c2$ZAL2.WS_mean <- as.numeric(bootdiffs.c2$ZAL2.WS_mean)
bootdiffs.c2$prob_xtrme <- as.numeric(bootdiffs.c2$prob_xtrme)
bootdiffs.c2$boot_t0 <- as.numeric(bootdiffs.c2$boot_t0)

#set up a wide format of the same data
bootdiffs.ws_columns <- c("position.1","t_pvalue.ws","Alt_mean.ws","ZAL2.WS_mean.ws","cohensD.ws",
                          "boot_t0.ws","boot_0.025.ws","boot_0.975.ws","IQR_span0.ws",
                          "H0_0.025.ws","H0_0.975.ws","prob_xtrme.ws","boot_test.ws")

colnames(r) <- bootdiffs.ws_columns
bootdiffs.c <- merge(r, q, by = "position.1")
write.csv(bootdiffs.c, file = "bootdiffs.c_shared.csv")

bootdiffs.c$per.meth_mean.2 <- (as.numeric(q$Alt_mean) + as.numeric(q$ZAL2.WS_mean))/2
bootdiffs.c$per.meth_mean.ws <- (as.numeric(r$Alt_mean) + as.numeric(r$ZAL2.WS_mean))/2
bootdiffs.c$prob_xtrme <- as.numeric(bootdiffs.c$prob_xtrme)
bootdiffs.c$prob_xtrme.ws <- as.numeric(bootdiffs.c$prob_xtrme.ws)

ggplot(bootdiffs.c, aes(x=prob_xtrme, y=per.meth_mean.2)) + 
  geom_point(show.legend = T, shape=16, size = 1 ) +
  geom_smooth(method = "lm", se = F, size = 1.5, show.legend = T) +
  theme_classic(base_size = 10)

stats_bootdiffs.zal2.lm.c <- lm(per.meth_mean.2 ~ prob_xtrme, data = subset(bootdiffs.c2,boot_test == "TS2 vs WS2"))
summary(stats_bootdiffs.zal2.lm.c)
anova(stats_bootdiffs.zal2.lm.c)

ggplot(bootdiffs.c, aes(x=prob_xtrme.ws, y=per.meth_mean.ws)) + 
  geom_point(show.legend = T, shape=16, size = 1 ) +
  geom_smooth(method = "lm", se = F, size = 1.5, show.legend = T) +
  theme_classic(base_size = 10)

stats_bootdiffs.ws.lm.c <- lm(per.meth_mean.ws ~ prob_xtrme, data = subset(bootdiffs.c2,boot_test == "WS2m vs WS2"))
summary(stats_bootdiffs.ws.lm.c)
anova(stats_bootdiffs.ws.lm.c)






################   METHYLATION FIGURES: Box plots       ####

### BOXPLOTS ###

setwd(results)

## Adults ##
# overall methylation: adults
FIG_allele.a <-  ggplot(allele.adult, aes(x=CpGs, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = T) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = T, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) + 
  scale_color_manual(values=c("black","black" ,"black"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,17.6))+
  geom_signif(y_position=c(15.3, 16.8), xmin=c(1, 2), xmax=c(1.25, 2.25), annotation=c("*","*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  scale_x_discrete(labels = c("per.meth" = "All CpGs","per.meth.s" = "Shared<br>CpGs Only")) +
  labs(x = " ", fill = "Allele", color = "Allele", title = "*VIP* CRE (Overall)") +
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(1.2)), axis.text.y = element_markdown(),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.9)), legend.title = element_markdown(size = rel(1)),
        text = element_text(size = 13, color = "black")
  )


# shared CpGs: adults

dfcg.zeda <- adult.dfcg[adult.dfcg$position.1 %in% cpgs_posthoc.a, ]

FIG_allele.sharedcpga <-  ggplot(dfcg.zeda, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = F, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500")) + 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,85))+
  geom_signif(y_position=c(31,21,41,80.5), xmin=c(3,4,6,7), xmax=c(3.25,4.25,6.25,7.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "Shared CpGs with a Significant Main Effect of Allele")+
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(size = rel(1.2)),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )


dfcg.polya2 <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2m, ]

FIG_allele.zal2.cpga <-  ggplot(dfcg.polya2, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = F, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500")) + 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,25))+
  geom_signif(y_position=c(13.2,14.4,10.2), xmin=c(1,2,3), xmax=c(1.25,2.25,3.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "VIP<sup>2</sup>-only CpGs") +
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(size = rel(1.2)),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )

dfcg.polya2m <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2, ]

FIG_allele.zal2m.cpga <-  ggplot(dfcg.polya2m, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = F, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500")) + 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,45))+
  geom_signif(y_position=c(11.9,25,40,11.3,5.7), xmin=c(1,2,3,4,5), xmax=c(1.25,2.25,3.25,4.25,5.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "VIP<sup>2m</sup>-only CpGs")+
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(size = rel(1.2)),
        axis.title.x = element_markdown(margin = margin(t=7)), axis.title.y = element_blank(),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )

#alternative method of printing multiple plots together
FIG_adult.boxes <- 
  arrangeGrob(FIG_allele.a, FIG_allele.sharedcpga, FIG_allele.zal2.cpga, FIG_allele.zal2m.cpga,
            ncol = 10, nrow = 3,
            layout_matrix = rbind(c(1,1,1,1,1,1,1,NA,NA,NA), c(2,2,2,2,2,2,2,2,2,2), c(3,3,3,3,4,4,4,4,4,4))) %>% 
  as_ggplot() + draw_plot_label(label = c("A", "B", "C","d"), size = 15, x = c(0, 0, 0, .4), y = c(1, 0.6667, 0.328, 0.328)) # Add labels 
  
ggsave("FIG_adultbox.tiff", plot = FIG_adult.boxes, units="in", width=6.5, height=8.25, dpi=300, compression = 'lzw')

## Nestlings ##
#Overall methylation: nestlings
FIG_allele.c <-  ggplot(allele.chick, aes(x=CpGs, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = T) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = T, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_color_manual(values=c("black","black" ,"black"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,11.5))+
  geom_signif(y_position=c(11.0), xmin=c(2), xmax=c(2.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  scale_x_discrete(labels = c("per.meth" = "All CpGs","per.meth.s" = "Shared<br>CpGs Only")) +
  labs(x = " ", fill = "Allele", color = "Allele", size = 14,  title = "VIP CRE (Overall)") +
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(1.2)), axis.text.y = element_markdown(),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.9)), legend.title = element_markdown(size = rel(1)),
        text = element_text(size = 13, color = "black")
  )

# shared CpGs: nestlings
dfcg.zedc <- chick.dfcg[chick.dfcg$position.1 %in% cpgs_posthoc.c, ]

FIG_allele.sharedcpgc <- ggplot(dfcg.zedc, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = FALSE, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500")) + 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,80))+
  geom_signif(y_position=c(47,23,55,74,11,10), xmin=c(1,2,3,4,5,6), xmax=c(1.25,2.25,3.25,4.25,5.25,6.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "Shared CpGs with a Significant Main Effect of Allele")+
  theme(strip.background = element_rect(fill = "white", linetype = NULL), 
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )
  
dfcg.polyc2 <- chick.dfcg[chick.dfcg$position.1 %in% cpg.2m, ]

FIG_allele.zal2.cpgc <-  ggplot(dfcg.polyc2, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = FALSE, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500"))+ 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(name="%5mC Methylation", limits = c(0,15))+
  geom_signif(y_position=c(8.8,7.9,10.9), xmin=c(1,2,3), xmax=c(1.25,2.25,3.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "VIP<sup>2</sup>-only CpGs")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(),
        axis.title.x = element_markdown(margin = margin(t=7)),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )
  
dfcg.polyc2m <- chick.dfcg[chick.dfcg$position.1 %in% cpg.2, ]

FIG_allele.zal2m.cpgc <-  ggplot(dfcg.polyc2m, aes(x=position.1, y=per.meth*100, color = allele.morph)) + 
  geom_boxplot(aes(fill=allele.morph), outlier.shape = NA, position = position_dodge(preserve = "single"), show.legend = F) +
  geom_point(position = position_jitterdodge(dodge.width = 0.75), aes(colour=allele.morph, fill = factor(allele.morph)), show.legend = FALSE, shape=21, size = 2, alpha = 1) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500")) + 
  scale_color_manual(values=c("black","black" ,"black")) +
  scale_y_continuous(limits = c(0,50), breaks = c(5,10,15,26.25,49.5)) +
  scale_y_break(c(17,25.75,26.75,49), ticklabels = NULL) +
  geom_signif(y_position=c(11,13.6,9.3,15.7,5.8), xmin=c(1,2,3,4,5), xmax=c(1.25,2.25,3.25,4.25,5.25), annotation=c("*"), tip_length=0, vjust = 0.4, color = "black", textsize=8)+ 
  labs(x = "CpG Position", title = "VIP<sup>2m</sup>-only CpGs")+
  theme(panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"),
        plot.background = element_rect(fill = "white"),
        plot.title = element_markdown(hjust = 0.2, size = rel(1.2)),
        axis.text.x = element_markdown(size = rel(.95)), axis.text.y = element_markdown(),
        axis.title.x = element_markdown(margin = margin(t=0),hjust = .6), axis.title.y = element_blank(), 
        axis.text.y.right = element_blank(), axis.ticks.y.right = element_blank(), axis.line.y.right = element_blank(),
        legend.text = element_markdown(size = rel(.8)), legend.title = element_markdown(size = rel(.9)),
        text = element_text(size = 13, color = "black")
  )

#stitch the plots together
FIG_nestling.boxes <- 
  arrangeGrob(FIG_allele.c, FIG_allele.sharedcpgc, print(FIG_allele.zal2.cpgc), print(FIG_allele.zal2m.cpgc),
                  ncol = 10, nrow = 3, layout_matrix = rbind(c(1,1,1,1,1,1,1,NA,NA,NA), c(2,2,2,2,2,2,2,2,2,2), c(3,3,3,3,4,4,4,4,4,4))) %>%
  as_ggplot() + draw_plot_label(label = c("A", "B", "C","D"), size = 15,  x = c(0, 0, 0, .4), y = c(1, 0.66, 0.34, 0.34)) # Add labels

ggsave("FIG_nestbox.tiff", plot = FIG_nestling.boxes, units="in", width=6.5, height=8.25, dpi=300, compression = 'lzw')

#cleanup
rm(dfcg.zeda, dfcg.zedc, dfcg.polya2m, dfcg.polya2, dfcg.polyc2, dfcg.polyc2m,
   gt, p)


################   METHYLATION FIGURES: Landscape Maps  ####

setwd(results)

#lolipops
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.shared,]
q <- unique(r$distfromTSSn)
k <- data.frame(
  x = q,
  y = rep(-2.8, 66)
)
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2,]
q <- unique(r$distfromTSSn)
l <- data.frame(
  x = q,
  y = rep(-2.8, 5)
)
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2m,]
q <- unique(r$distfromTSSn)
m <- data.frame(
  x = q,
  y = rep(-2.8, 3)
)

#ALL three alleles on the same map
#removing the allele specific, zero values for this one

## adults ##
#zal2 specific cpgs
a <- subset(adult.dfcg, position.1 == "2258582" & allele != "ZAL2m")
b <- subset(adult.dfcg, position.1 == "2259007" & allele != "ZAL2m")
c <- subset(adult.dfcg, position.1 == "2259224" & allele != "ZAL2m")

#zal2m speific cpgs
d <- subset(adult.dfcg, position.1 == "2258417" & allele != "ZAL2")
e <- subset(adult.dfcg, position.1 == "2258431" & allele != "ZAL2")
f <- subset(adult.dfcg, position.1 == "2258616" & allele != "ZAL2")
g <- subset(adult.dfcg, position.1 == "2258795" & allele != "ZAL2")
h <- subset(adult.dfcg, position.1 == "2259319" & allele != "ZAL2")

j <- subset(adult.dfcg, note == "shared")

adult.dfcg.r <- rbind(a,b,c,d,e,f,g,h,j)
mean(adult.dfcg.r$per.meth)

FIG_adultmapthree <- ggplot(data = adult.dfcg.r, aes(x=distfromTSSn, y=per.meth*100))+
  stat_smooth(data = adult.dfcg.r,method="loess",span=0.3, se=F, aes(fill=allele.morph, color = allele.morph), alpha=0.4) +
  geom_point(data = adult.dfcg.r,shape=21, size = 1, mapping=aes(fill=allele.morph, color = allele.morph), show.legend = T) +  
  scale_color_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-7,100), expand = c(0, 0))+
  labs(x = "Relative position", y = "%5mC Methylation") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"), legend.position = c(0.9,0.7), legend.text = element_markdown(),
        axis.text.x = element_markdown(), axis.text.y = element_markdown(),
        text = element_text(size = 17, color = "black")) +
  geom_rug(sides="b", color= "black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 2, color = "black") +
  geom_point(data = l, aes(x=x, y=y), shape = 19, size = 2, color = "#CA3500") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 2, color = "#2142D7") +
  annotate("segment", x = 0, y = -7, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0) ,clip = 'off')  # This keeps the labels from disappearing

ggsave("FIG_adultmapthree.tiff", plot = FIG_adultmapthree, units="in", width=12, height=4.2, dpi=300, compression = 'lzw')


## Nestlings ##
#zal2 specific cpgs
a <- subset(chick.dfcg, position.1 == "2258582" & allele != "ZAL2m")
b <- subset(chick.dfcg, position.1 == "2259007" & allele != "ZAL2m")
c <- subset(chick.dfcg, position.1 == "2259224" & allele != "ZAL2m")

#zal2m specific cpgs
d <- subset(chick.dfcg, position.1 == "2258417" & allele != "ZAL2")
e <- subset(chick.dfcg, position.1 == "2258431" & allele != "ZAL2")
f <- subset(chick.dfcg, position.1 == "2258616" & allele != "ZAL2")
g <- subset(chick.dfcg, position.1 == "2258795" & allele != "ZAL2")
h <- subset(chick.dfcg, position.1 == "2259319" & allele != "ZAL2")

j <- subset(chick.dfcg, note == "shared")

chick.dfcg.r <- rbind(a,b,c,d,e,f,g,h,j)
mean(chick.dfcg.r$per.meth)

FIG_chickmapthree <- ggplot(data = chick.dfcg.r, aes(x=distfromTSSn, y=per.meth*100))+
  stat_smooth(data = chick.dfcg.r,method="loess",span=0.3, se=F, aes(fill=allele.morph, color = allele.morph), alpha=0.4) +
  geom_point(data = chick.dfcg.r,shape=21, size = 1, mapping=aes(fill=allele.morph, color = allele.morph), show.legend = T) +  
  scale_color_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_fill_manual(values=c("#CDA434","#2142D7" ,"#CA3500"), name = "Allele", labels = c(
    "TS-VIP<sup>2</sup>",
    "WS-VIP<sup>2</sup>",
    "WS-VIP<sup>2m</sup>"
  )) +
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-7,100), expand = c(0, 0))+
  labs(x = "Relative position", y = "%5mC Methylation",  fill = "Allele", color = "Allele") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        legend.background = element_rect(color = "black", fill = "white"), legend.position = c(0.9,0.7), legend.text = element_markdown(),
        axis.text.x = element_markdown(), axis.text.y = element_markdown(),
        text = element_text(size = 17, color = "black")) +
  geom_rug(sides="b", color= "black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 2, color = "black") +
  geom_point(data = l, aes(x=x, y=y), shape = 19, size = 2, color = "#CA3500") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 2, color = "#2142D7") +
  annotate("segment", x = 0, y = -7, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0) ,clip = 'off')  # This keeps the labels from disappearing

ggsave("FIG_nestmapthree.tiff", plot = FIG_chickmapthree, units="in", width=12, height=4.2, dpi=300, compression = 'lzw')

rm(a,b,c,d,e,f,g,h,i,j)


################   METHYLATION FIGURES: Heatmaps        ####

setwd(results)

#compililng cohen's d's
#need to note direction of effect
#doing two pairwise tests per age group

cds_WS <- function (x) {
  c(
    cohens_d(per.meth1 ~ allele, data = x)
  )
}

cds_ZAL2 <- function (x) {
  c(
    cohens_d(per.meth1 ~ morph, data = x)
  )
}

a <- do.call("rbind",lapply(split(adult.dfcg_bymorph[["White"]],adult.dfcg_bymorph[["White"]]$position.1), cds_WS))
b <- a %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
b$age = "Adult"
b$test = "ZAL2-WS vs ZAL2m-WS"

i <- split(adult.dfcg_by2.2m[["ZAL2"]],adult.dfcg_by2.2m[["ZAL2"]]$position.1)
i <- i[names(i) %in% cpg.tt]

c <- do.call("rbind",lapply(i, cds_ZAL2))
d <- c %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
d$age = "Adult"
d$test = "ZAL2-WS vs ZAL2-TS"

e <- do.call("rbind",lapply(split(chick.dfcg_bymorph[["White"]],chick.dfcg_bymorph[["White"]]$position.1), cds_WS))
f <- e %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
f$age = "Nestling"
f$test = "ZAL2-WS vs ZAL2m-WS"

i <- split(chick.dfcg_by2.2m[["ZAL2"]],chick.dfcg_by2.2m[["ZAL2"]]$position.1)
i <- i[names(i) %in% cpg.tt]

g <- do.call("rbind",lapply(i, cds_ZAL2))
h <- g %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Position")
h$age = "Nestling"
h$test = "ZAL2-WS vs ZAL2-TS"

cds <- rbind (b,d,f,h)

j <- as.data.frame(cpg.2)
colnames(j) = "Position"
j$Cohens_d <- NA
j$CI <- NA
j$CI_low <- NA
j$CI_high <- NA
j$test <- "ZAL2-WS vs ZAL2-TS"
k <- j
j$age = "Adult"
k$age = "Nestling"
cds <- rbind(cds,j,k)
cds$Cohens_d <- as.numeric(cds$Cohens_d)
cds$test <- as.factor(cds$test)


#set up facets
cds$facets = factor(cds$test, levels = c("ZAL2-WS vs ZAL2-TS","ZAL2-WS vs ZAL2m-WS"), labels = c(
  "WS-VIP<sup>2</sup><br>vs. TS-VIP<sup>2</sup>",
  "WS-VIP<sup>2</sup><br>vs. WS-VIP<sup>2m</sup>"
))

sigcds <- data.frame( 
  facets = c(rep("WS-VIP<sup>2</sup><br>vs. WS-VIP<sup>2m</sup>", times = 26)), 
  x = c("Adult","Adult","Adult","Adult",
        "Nestling","Nestling","Nestling","Nestling","Nestling","Nestling",
        "Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult",
        "Nestling","Nestling","Nestling","Nestling","Nestling","Nestling","Nestling","Nestling"),
  y = c("2258863","2258890","2259049","2259083",
        "2258341","2259083","2259341","2259049","2258440","2259236",
        "2258417","2258431","2258616","2258795","2259319","2258582","2259007","2259224",
        "2258417","2258431","2258616","2258795","2259319","2258582","2259007","2259224"),
  text = c("*","*","*","*",
           "*","*","*","*","*","*",
           "*","*","*","*","*","*","*","*",
           "*","*","*","*","*","*","*","*"))

sigcds2 <- data.frame(
  facets = "WS-VIP<sup>2</sup><br>vs. WS-VIP<sup>2m</sup>",
  x = "Nestling",
  y = "2259007",
  text = "*")

labscds.2m <- data.frame(
  facets = "WS-VIP<sup>2</sup><br>vs. TS-VIP<sup>2</sup>",
  x = c(-.55),
  y = cpg.2,
  xend = c(.4),
  yend = cpg.2
)

labscds.2 <- data.frame(
  facets = "WS-VIP<sup>2</sup><br>vs. TS-VIP<sup>2</sup>",
  x = c(-.55),
  y = cpg.2m,
  xend = c(.4),
  yend = cpg.2m
)


#Heatmap of Cohen's D
FIG_cpg.heatmap <- ggplot(data = cds, mapping = aes(x = age, y = Position)) +
  facet_wrap(~facets, labeller = label_value) +
  geom_tile(data = subset(cds, test == "ZAL2-WS vs ZAL2-TS"), aes(fill = Cohens_d), color = "white", show.legend = TRUE) +
  scale_fill_gradient2(low = "#CDA434", mid = "white", high = "#2142D7", limits = c(-1.5,1), na.value = NA) +
  labs(x = "Age", y = "Position", fill = "Cohen's d<br>(left)")+
  #to get different color scales for each allele
  new_scale("fill") +
  geom_tile(data = subset(cds, test == "ZAL2-WS vs ZAL2m-WS"), aes(fill = Cohens_d), color = "white", show.legend = TRUE) +
  scale_fill_gradient2(low = "#CA3500", mid = "white", high = "#2142D7", na.value = NA) +
  scale_y_discrete(limits=rev)+
  labs(x = "Age", y = "Position", fill = "Cohen's d<br>(right)") +
  #asterisks at significant post-hoc tests
  geom_text(data = sigcds, aes(x = x, y = y, label = text), nudge_y = -.35, size = 10, color = "black") +
  geom_text(data = sigcds2, aes(x = x, y = y, label = text), nudge_y = -.35, size = 10, color = "white") +
  #highlights for polymorphic sites
  coord_cartesian(xlim = c(1,2) ,clip = 'off') +  # This keeps the labels from disappearing
  geom_segment(data = labscds.2m, aes(x = x, y = y, xend = xend, yend = yend), color = "#CA3500", size = 5, alpha = 0.5) + 
  geom_segment(data = labscds.2, aes(x = x, y = y, xend = xend, yend = yend), color = "#2142D7", size = 5, alpha = 0.5) +
  theme(strip.background = element_rect(fill = "white", linetype = NULL), strip.text = element_markdown(size = rel(1.5)),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed"), 
        panel.border = element_rect(color = "black", fill = NA, linetype = "solid"),
        plot.background = element_rect(fill = "white"), plot.margin = unit(c(1,1,1,2), "lines"),
        axis.text.y = element_markdown(size = rel(1.3)), axis.title.x = element_blank(),
        axis.text.x = element_markdown(size = rel(1.75)), axis.title.y = element_markdown(size = rel(1.75)), 
        legend.title = element_markdown(size = rel(1.2)), legend.text = element_text(size = rel(1)),
        text = element_text(size = 12.5, color = "black")
        )
  
ggsave("FIG_cpg.heatmap.tiff", plot = FIG_cpg.heatmap, units="in", width=6.5, height=16.5, dpi=300, compression = 'lzw')

rm(a,b,c,d,e,f,g,h,i,j,k)
rm(labscds.2, labscds.2m, sigcds, sigcds2)

####    Methylation and Expression Correlations       #####################
setwd(results)

###   Correlation by allele in Adults   ###
lmemeans_byallele.a <- function (x) {
  x <- subset(x, CpGs == "per.meth")
  a <- lm(log2(norm.count)~log(per.meth) + sex,
           data = x)
  c(
    anova(a)[1,4],
    anova(a)[1,5],
    summary(a)$r.squared,
    anova(a)[1,1],
    nrow(x)
  )
}
zed <- c("F-stat","p-value","R2","df","n")

c <- do.call("rbind",lapply(split(brcorr.means.a, brcorr.means.a$allele.morph),lmemeans_byallele.a))
colnames(c) = zed
stats_regs.adults <- c %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Allele")

###   Nestlings   ###
lmemeans_byallele.c <- function (x) {
  x <- subset(x, CpGs == "per.meth")
  a <- lme(log2(norm.count)~log(per.meth) + sex,
           random = ~1|nest,
           data = x)
  b <- performance::r2(a)[["R2_conditional"]]
  attributes(b) <- NULL
  c(
    anova(a)[2,3],
    anova(a)[2,4],
    b,
    anova(a)[2,1],
    anova(a)[2,2],
    nrow(x)
  )
}
zed <- c("F-stat","p-value","R2","df_num","df_den","n")

c <- do.call("rbind",lapply(split(brcorr.means.c, brcorr.means.c$allele.morph),lmemeans_byallele.c))
colnames(c) = zed
stats_regs.chicks <- c %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Allele")

rm(zed,c)

######    Regressions by each Cpg                       ###########

cpgR_adults <- function (x) {
  corr <- lm(log2(norm.count)~per.meth1 + sex, data = x)
  w <- summary(corr)
  w2 <- pf(w$fstatistic[1],w$fstatistic[2],w$fstatistic[3],lower.tail = F)
  attributes(w2) <- NULL
  c (
    mean(x$per.meth), #add actual percent methylation
    mean(x$per.meth1), #add transformed percent methylation
    as.data.frame(w$coefficients)[2,1], #coefficient estimate for percent methylation
    anova(corr)[1,4], #F-value
    anova(corr)[1,5], #p-value
    w$r.squared, #R2 for model
    nrow(x)
  )
}
cpgR_columns <- c("per.meth","per.meth1","coef_est","permeth_F","permeth_p","R2","n")

d <- split(brcorr.dfcg.a_byallele[["ZAL2-TS"]], brcorr.dfcg.a_byallele[["ZAL2-TS"]]$position.1, drop = TRUE) 
e <- d[names(d) %in% cpg.tt]
cpgR_tsa <- do.call("rbind",lapply(e, cpgR_adults))
f <- split(brcorr.dfcg.a_byallele[["ZAL2-WS"]], brcorr.dfcg.a_byallele[["ZAL2-WS"]]$position.1, drop = TRUE)
g <- f[names(f) %in% cpg.tt]
cpgR_ws2a <- do.call("rbind",lapply(g, cpgR_adults))
h <- split(brcorr.dfcg.a_byallele[["ZAL2m-WS"]], brcorr.dfcg.a_byallele[["ZAL2m-WS"]]$position.1, drop = TRUE)
i <- h[names(h) %in% cpg.ttm]
cpgR_ws2ma <- do.call("rbind",lapply(i, cpgR_adults)) 

colnames(cpgR_tsa) <- cpgR_columns
cpgR_tsa <- cpgR_tsa %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("position.1")

colnames(cpgR_ws2a) <- cpgR_columns
cpgR_ws2a <- cpgR_ws2a %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("position.1")

colnames(cpgR_ws2ma) <- cpgR_columns
cpgR_ws2ma <- cpgR_ws2ma %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("position.1")

#And now in chicks

cpgR_chicks <- function (x) {
  corr <- lme(log2(norm.count)~per.meth1 + sex, random = ~1|nest , data = x)
  w <- summary(corr)
  r2 <- performance::r2(corr)[["R2_conditional"]]
  attributes(r2) = NULL
  c (
    mean(x$per.meth), #add actual percent methylation
    mean(x$per.meth1), #add transformed percent methylation
    as.data.frame(coef(w))[2,1], #coefficient estimate per.meth
    anova(corr)[2,3], #F-stat
    anova(corr)[2,4], #p-value
    r2, #conditional R2 for model
    nrow(x) #n
  )
}

j <- split(brcorr.dfcg.c_byallele[["ZAL2-TS"]], brcorr.dfcg.c_byallele[["ZAL2-TS"]]$position.1, drop = TRUE)
k <- j[names(j) %in% cpg.tt]
cpgR_tsc <- do.call("rbind",lapply(k, cpgR_chicks))
j <- split(brcorr.dfcg.c_byallele[["ZAL2-WS"]],brcorr.dfcg.c_byallele[["ZAL2-WS"]]$position.1, drop = TRUE)
k <- j[names(j) %in% cpg.tt]
cpgR_ws2c <- do.call("rbind",lapply(k, cpgR_chicks))
j <- split(brcorr.dfcg.c_byallele[["ZAL2m-WS"]],brcorr.dfcg.c_byallele[["ZAL2m-WS"]]$position.1, drop = TRUE)
k <- j[names(j) %in% cpg.ttm]
cpgR_ws2mc <- do.call("rbind",lapply(k, cpgR_chicks))

colnames(cpgR_tsc) <- cpgR_columns
cpgR_tsc <- cpgR_tsc %>% as.data.frame() %>% tibble::rownames_to_column("position.1")
colnames(cpgR_ws2c) <- cpgR_columns
cpgR_ws2c <- cpgR_ws2c %>% as.data.frame() %>% tibble::rownames_to_column("position.1")
colnames(cpgR_ws2mc) <- cpgR_columns
cpgR_ws2mc <- cpgR_ws2mc%>% as.data.frame() %>% tibble::rownames_to_column("position.1")

a <- lme(log2(norm.count)~per.meth1 + sex, random = ~1|nest , data = j[["2258839"]])
anova(a)

#combine dfs and tidy up

cpgR_tsa3 <- merge(cpgR_tsa, dists, by = "position.1", all.x = F)
cpgR_tsa3$distfromTSSn <- as.numeric(cpgR_tsa3$NW.rel.pos)

cpgR_tsc3 <- merge(cpgR_tsc, dists, by = "position.1", all.x = F)
cpgR_tsc3$distfromTSSn <- as.numeric(cpgR_tsc3$NW.rel.pos)

cpgR_ws2a3 <- merge(cpgR_ws2a, dists, by = "position.1", all.x = F)
cpgR_ws2a3$distfromTSSn <- as.numeric(cpgR_ws2a3$NW.rel.pos)

cpgR_ws2c3 <- merge(cpgR_ws2c, dists, by = "position.1", all.x = F)
cpgR_ws2c3$distfromTSSn <- as.numeric(cpgR_ws2c3$NW.rel.pos)

cpgR_ws2ma3 <- merge(cpgR_ws2ma, dists, by = "position.1", all.x = F)
cpgR_ws2ma3$distfromTSSn <- as.numeric(cpgR_ws2ma3$NW.rel.pos)

cpgR_ws2mc3 <- merge(cpgR_ws2mc, dists, by = "position.1", all.x = F)
cpgR_ws2mc3$distfromTSSn <- as.numeric(cpgR_ws2mc3$NW.rel.pos)

write.csv(cpgR_tsa3, "stats_cpgR_tsa.csv")
write.csv(cpgR_tsc3, "stats_cpgR_tsc.csv")
write.csv(cpgR_ws2a3, "stats_cpgR_ws2a.csv")
write.csv(cpgR_ws2c3, "stats_cpgR_ws2c.csv")
write.csv(cpgR_ws2ma3, "stats_cpgR_ws2ma.csv")
write.csv(cpgR_ws2mc3, "stats_cpgR_ws2mc.csv")

rm(d,e,f,g,h,i,j,k,cpgR_tsa,cpgR_tsc,cpgR_ws2a,cpgR_ws2c,cpgR_ws2ma,cpgR_ws2mc)

###############   REGRESSIONS FIGURES: Landscape Maps   ####
#lolipops
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.shared,]
q <- unique(r$distfromTSSn)
k <- data.frame(
  x = q,
  y = rep(-.045, 66)
)
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2,]
q <- unique(r$distfromTSSn)
l <- data.frame(
  x = q,
  y = rep(-.045, 5)
)
r <- adult.dfcg[adult.dfcg$position.1 %in% cpg.2m,]
q <- unique(r$distfromTSSn)
m <- data.frame(
  x = q,
  y = rep(-.045, 3)
)

R_sigs <- data.frame(
  x = c(-374,-228,-188,-184,-172,-146,-143,-62),
  y = c(0.87,0.9,0.91,1.01,1.0,1.03,1.0,0.95),
  text = "*")


## ADULTS ##
FIG_R2.tsmap.a <- ggplot(cpgR_tsa3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.4, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1.05), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "TS-VIP^2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
        )+
  #signifiance Rs
  geom_text(data = R_sigs, aes(x = x, y = y, label = text), size = 7, color = "black") +
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing


FIG_R2.wsmap.a <- ggplot(cpgR_ws2a3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.3, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "WS-VIP^2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
  )+
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing

FIG_R2.mwsmap.a <- ggplot(cpgR_ws2ma3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.3, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "WS-VIP^2m") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
  )+
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing

FIG_adult.R2map <- 
  arrangeGrob(FIG_R2.tsmap.a, FIG_R2.wsmap.a, FIG_R2.mwsmap.a, 
              ncol = 1, nrow = 3,
              layout_matrix = rbind(1,2,3)) %>% 
  as_ggplot() + draw_plot_label(label = c("A", "B", "C"), size = 20, x = c(0, 0, 0), y = c(1.005, 0.67, 0.338)) # Add labels 

ggsave("FIG_regsmap.adult.tiff", plot = FIG_adult.R2map, units="in", width=8.5, height=11, dpi=300, compression = 'lzw')


## NESTLINGS ##
FIG_R2.tsmap.c <- ggplot(cpgR_tsc3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.3, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1.05), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "TS-VIP^2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
  )+
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing


FIG_R2.wsmap.c <- ggplot(cpgR_ws2c3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.3, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "WS-VIP^2") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
  )+
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing


FIG_R2.mwsmap.c <- ggplot(cpgR_ws2mc3, aes(x=distfromTSSn, y=R2))+
  geom_bar(stat = "identity")+
  stat_smooth(method="loess",span=0.3, se=F, alpha=0.4) +
  scale_color_manual(values=c("#375590","#CF0A27","black")) + 
  scale_fill_manual(values=c("#375590","#CF0A27", "black"))+
  scale_x_continuous(breaks = c(-2000, -1500, -1000, -500, 0)) +
  scale_y_continuous(limits=c(-0.1,1), expand = c(0, 0), breaks = seq(-.5, 1, 0.5))+
  labs(x = NULL, y = "R^2", title = "WS-VIP^2m") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.line = element_line(color = "black"), axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,1,0,0), "lines"), # This widens the right margin
        plot.title = element_markdown(), axis.title.y = element_markdown(),
        text = element_text(colour = "black", size = 15)
  )+
  #lolipops
  geom_rug(sides="b",color="black") +
  geom_point(data = k, aes(x=x, y=y), shape = 21, size = 3, color = "black") +
  geom_point(data = m, aes(x=x, y=y), shape = 19, size = 3, color = "#2142D7") +
  #TSS note
  annotate("segment", x = 0, y = -0.1, xend = 0, yend = 0) +
  annotate("segment", x = 0, y = 0, xend = 30, yend = 0, arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
  annotate("text", x=25, y=0, label = "TSS", hjust = -.1, color = "black", size = 4, fontface = "bold") +
  coord_cartesian(xlim = c(-1250,0), ylim = c(-0.1,1.05), clip = 'off')  # This keeps the labels from disappearing


FIG_nestling.R2map <- 
  arrangeGrob(FIG_R2.tsmap.c, FIG_R2.wsmap.c, FIG_R2.mwsmap.c,  
              ncol = 1, nrow = 3,
              layout_matrix = rbind(1,2,3)) %>% 
  as_ggplot() + draw_plot_label(label = c("A", "B", "C"), size = 20, x = c(0, 0, 0), y = c(1.005, 0.67, 0.338)) # Add labels 

ggsave("FIG_regsmap.nest.tiff", plot = FIG_nestling.R2map, units="in", width=8.5, height=11, dpi=300, compression = 'lzw')




###############   REGRESSIONS FIGURES: Heatmaps of R2   #############

cpgR_tsa3$age <- "Adult"
cpgR_ws2a3$age <- "Adult"
cpgR_ws2ma3$age <- "Adult"

cpgR_tsc3$age <- "Nestling"
cpgR_ws2c3$age <- "Nestling"
cpgR_ws2mc3$age <- "Nestling"

cpgR_tsa3$allele <- "TS-ZAL2"
cpgR_ws2a3$allele <- "WS-ZAL2"
cpgR_ws2ma3$allele <- "WS-ZAL2m"

cpgR_tsc3$allele <- "TS-ZAL2"
cpgR_ws2c3$allele <- "WS-ZAL2"
cpgR_ws2mc3$allele <- "WS-ZAL2m"

stats_cpg.R2a <- rbind(cpgR_tsa3,
                 cpgR_ws2a3,
                 cpgR_ws2ma3)

stats_cpg.R2c <- rbind(cpgR_tsc3,
                 cpgR_ws2c3,
                 cpgR_ws2mc3)

stats_cpg.R2 <- rbind(stats_cpg.R2c,stats_cpg.R2a)

stats_cpg.R2$R2_dir <- (stats_cpg.R2$R2*stats_cpg.R2$coef_est)/abs(stats_cpg.R2$coef_est)

stats_cpg.R2$facets = factor(stats_cpg.R2$allele, labels = c(
  "TS-VIP^2",
  "WS-VIP^2",
  "WS-VIP^2m"
))

r2_sigs <- data.frame(
  facets = c("TS-VIP^2","TS-VIP^2","TS-VIP^2","TS-VIP^2","TS-VIP^2","TS-VIP^2","TS-VIP^2","TS-VIP^2"),
  x = c("Adult","Adult","Adult","Adult","Adult","Adult","Adult","Adult"),
  y = c("2259300","2259338","2259312","2259341","2259422","2259110","2259296","2259256" ),
  text = c("*","*","*","*","*","*","*","*")
)

labscds.2m <- data.frame(
  facets = "TS-VIP^2",
  x = c(-1.4),
  y = cpg.2,
  xend = c(.4),
  yend = cpg.2
)

labscds.2 <- data.frame(
  facets = "TS-VIP^2",
  x = c(-1.4),
  y = cpg.2m,
  xend = c(.4),
  yend = cpg.2m
)

FIG_R2.heatmap <- ggplot(data = stats_cpg.R2, mapping = aes(x = age, y = position.1)) +
  geom_tile(data = stats_cpg.R2, aes(fill = R2_dir), color = "white") +
  facet_grid(~facets, labeller = label_value)+
  labs(y = "Postion", fill = "R^2") +
  scale_fill_gradient2(low = "#025669", mid = "white", high = "#ED760E", limits = c(-1,1)) +
  scale_y_discrete(limits=rev)+
  geom_text(data = r2_sigs, aes(x = x, y = y, label = text), nudge_y = -.4, size = 10) +
  theme(strip.background = element_rect(fill = "white", linetype = NULL), strip.text = element_markdown(size = rel(1.1)),
        panel.background = element_rect(fill = "white"), plot.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "grey", linetype = "dashed"), 
        panel.border = element_rect(color = "black", fill = NA, linetype = "solid"), plot.margin = unit(c(0,0,0,1), "lines"),
        axis.text.y = element_markdown(size = rel(1.3)), axis.title.y = element_markdown(size = rel(1.5)),
        axis.text.x = element_markdown(size = rel(1.2), angle = 45, hjust = 1), axis.title.x = element_blank(),
        legend.title = element_markdown(size = rel(1.5)), legend.text = element_text(size = rel(.9)),
        text = element_text(size = 15, color = "black")
        ) +
  coord_cartesian(xlim = c(1,2) ,clip = 'off') +  # This keeps the labels from disappearing
  geom_segment(data = labscds.2m, aes(x = x, y = y, xend = xend, yend = yend), color = "#CA3500", size = 6.5, alpha = 0.5) + 
  geom_segment(data = labscds.2, aes(x = x, y = y, xend = xend, yend = yend), color = "#2142D7", size = 7, alpha = 0.5)
    
ggsave("FIG_Rheatmap.tiff", plot = FIG_R2.heatmap, units="in", width=6, height=15, dpi=300, compression = 'lzw')






#the end.