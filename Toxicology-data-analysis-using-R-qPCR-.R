library(tidyverse)
library(readxl)
library(magrittr)
library(tibble)
library(multcomp)


Chemical_A.geneA <- read_excel("./Raw data example (qPCR).xlsx")
Chemical_A.geneA

## For Shapiro & Bartlette and Statistics ##

dCt.Chemical_A.geneA <- Chemical_A.geneA$`delCt (geneA-housekeeping gene)`
sample.Chemical_A.geneA <- Chemical_A.geneA$"Sample Name"
df.Chemical_A.geneA <- data.frame(sample = sample.Chemical_A.geneA, dCt = dCt.Chemical_A.geneA)
df.Chemical_A.geneA$conc <- factor(rep(c(0, 1, 10),times = c(4, 4, 4)))

## For visualization ##

v.Chemical_A.geneA <- na.omit(Chemical_A.geneA)
names(v.Chemical_A.geneA) <- c("sample", "dCt")

#####################################################################################################

## Shapiro-Wilk analysis ##
# If p-value > 0.05, then the dataset has normality

# First, make groups: Control, Chemical_A_1uM, and Chemical_A_10uM

df.Chemical_A.geneA

for(i in 1:nrow(df.Chemical_A.geneA)){
  df.Chemical_A.geneA[i, 1] <- gsub(".{2}$", "", df.Chemical_A.geneA[i, 1])
}

# Then perform Shapiro-Wilk test by a group

distinct_names <- df.Chemical_A.geneA %>% dplyr::distinct(sample)

for(j in 1:nrow(distinct_names)){
  print(distinct_names[j, ])
  print(shapiro.test((df.Chemical_A.geneA %>% dplyr::filter(sample == distinct_names[j, ]))$dCt))
}
    

## Homogeneity of variance (Kruskal-Wallis) ##
# If p-value > 0.05, then the dataset has homogeneity of variance

print(bartlett.test(dCt ~ conc, data = df.Chemical_A.geneA))


## Statistical significance (Kruskal-Wallis) ##
# If p-value < 0.05, then the dataset has significance

print(kruskal.test(dCt ~ conc, data = df.Chemical_A.geneA))

#####################################################################################################

## ddCt calculation ##
# For visualization, we need to convert dCt to ddCt

v.Chemical_A.geneA

v.Chemical_A.geneA$conc <- factor(rep(c(0, 1, 10),times = c(4, 4, 4))) 
v.Chemical_A.geneA

dct.mean.DMSO <- mean(v.Chemical_A.geneA$dCt[v.Chemical_A.geneA$conc=="0"])
v.Chemical_A.geneA$ddCt <- v.Chemical_A.geneA$dCt -dct.mean.DMSO
v.Chemical_A.geneA$express <- 2^-(v.Chemical_A.geneA$ddCt)

v.Chemical_A.geneA <- v.Chemical_A.geneA %>% group_by(conc) %>% mutate(ddctmean=mean(ddCt),ddctse=plotrix::std.error(ddCt))
v.Chemical_A.geneA$ddctse_max <- v.Chemical_A.geneA$ddctmean + v.Chemical_A.geneA$ddctse
v.Chemical_A.geneA$ddctse_min <- v.Chemical_A.geneA$ddctmean - v.Chemical_A.geneA$ddctse
v.Chemical_A.geneA$expr.mean <- 2^-(v.Chemical_A.geneA$ddctmean) 
v.Chemical_A.geneA$expr.se_max <- 2^-(v.Chemical_A.geneA$ddctse_max) 
v.Chemical_A.geneA$expr.se_min <- 2^-(v.Chemical_A.geneA$ddctse_min)

v.Chemical_A.geneA.bar <- v.Chemical_A.geneA %>% dplyr::distinct(conc,expr.mean,expr.se_min,expr.se_max)
v.Chemical_A.geneA.se <- v.Chemical_A.geneA %>% dplyr::distinct(conc,expr.mean,ddctse)
v.Chemical_A.geneA.se
v.Chemical_A.geneA.se$expr.mean <- signif(v.Chemical_A.geneA.se$expr.mean, 2)
v.Chemical_A.geneA.se$ddctse <- signif(v.Chemical_A.geneA.se$ddctse, 2)
v.Chemical_A.geneA.se

v.Chemical_A.geneA.bar_2 <- v.Chemical_A.geneA.bar 
options(pillar.sigfig = 2)
as_tibble(v.Chemical_A.geneA.bar_2)


## Visualization ##

ggplot() +
  xlab("Chemical_A concentration")+
  ylab("relative expression of geneA")+
  layer(
    data=v.Chemical_A.geneA.bar, 
    mapping=aes(x=conc, y=expr.mean), 
    geom="bar",
    params = list(fill = "forestgreen", alpha = 0.5),
    stat="identity", 
    position="identity"
  )+
  layer(
    data=v.Chemical_A.geneA.bar, 
    mapping=aes(x=conc, ymin = expr.se_max, ymax = expr.se_min, width = 0.5), 
    geom="errorbar",
    stat = "identity",
    position="identity"
  )+
  layer(
    data=v.Chemical_A.geneA, 
    mapping=aes(x=conc, y=express), 
    geom="point", 
    stat="identity", 
    position="identity",
    params = list(shape = 21, fill = "white", color = "red" ,size = 2.5, stroke = 1.5)
  )+
  ylim(0,max(round(v.Chemical_A.geneA$express, digits = -1))+10)

#####################################################################################################

## Statistics to check asterisk ##
# If p-value < 0.05, then the dataset has significance

vx <- df.Chemical_A.geneA$dCt
fx <- df.Chemical_A.geneA$conc

summary(aov(vx ~ fx))

summary(multcomp::glht(aov(vx~fx),linfct=multcomp::mcp(fx="Dunnett")))