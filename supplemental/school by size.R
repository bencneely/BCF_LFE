#--------------------------------------------------------------
#Ben Neely
#07/24/2025
#Check out whether BCF school by size
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
## Checks if package is installed, installs if not, activates for current session
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)

if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)

if("car" %in% rownames(installed.packages()) == FALSE) {install.packages("car")}
library(car)

if("vegan" %in% rownames(installed.packages()) == FALSE) {install.packages("vegan")}
library(vegan)

if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)

if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)

## Set seed to make results reproducible
RNGkind(kind="Mersenne-Twister",normal.kind="Inversion",sample.kind="Rejection")
set.seed(724)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22,color="black",face="bold"),
        axis.text=element_text(size=18,color="black"),
        legend.position="inside",
        legend.position.inside=c(0.99,0.99),
        legend.justification=c("right","top"),
        legend.text=element_text(size=20),
        legend.title=element_blank())
options(scipen=999)

## Read in data with import
fish=import("bcf sampling data.xlsx",which="fish")%>%
  expandCounts(~count)

################################################################################
## Use permutation tests to determine if observed variation in fish length
## in grids differs from a random distribution of lengths
## using as a surrogate to test if Blue Catfish school by size.
## The thought is that if observed variability < random variability
## there is evidence that fish school by size
################################################################################

################################################################################
## Clinton
cltr=filter(fish,impd=="CLTR")

## Calculate CV of mean length for each grid
cltr_obs_cv=cltr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
cltr_obs_mean_cv=mean(cltr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
cltr_perm_mean_cvs=replicate(n_perm, {
  cltr$tl_perm=sample(cltr$tl)
  cltr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
cltr_perm_mean=mean(cltr_perm_mean_cvs)
cltr_perm_sd=sd(cltr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
cltr_obs_z=(cltr_obs_mean_cv-cltr_perm_mean)/cltr_perm_sd
cltr_perm_z=(cltr_perm_mean_cvs-cltr_perm_mean)/cltr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(cltr_p=mean(cltr_perm_z <= cltr_obs_z))

## Plot mean CVs and observed CV
cltr_plot=ggplot()+
  geom_histogram(aes(x=cltr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#4E79A7")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression(""))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression("Proportion of permutations"))+
  geom_vline(aes(xintercept=cltr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Clinton"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)<0.001),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## El Dorado
eldr=filter(fish,impd=="ELDR")

eldr_obs_cv=eldr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
eldr_obs_mean_cv=mean(eldr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
eldr_perm_mean_cvs=replicate(n_perm, {
  eldr$tl_perm=sample(eldr$tl)
  eldr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
eldr_perm_mean=mean(eldr_perm_mean_cvs)
eldr_perm_sd=sd(eldr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
eldr_obs_z=(eldr_obs_mean_cv-eldr_perm_mean)/eldr_perm_sd
eldr_perm_z=(eldr_perm_mean_cvs-eldr_perm_mean)/eldr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(eldr_p=mean(eldr_perm_z <= eldr_obs_z))

## Plot mean CVs and observed CV
eldr_plot=ggplot()+
  geom_histogram(aes(x=eldr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#F28E2B")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression(""))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression(""))+
  geom_vline(aes(xintercept=eldr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("El Dorado"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)==0.026),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Melvern
melr=filter(fish,impd=="MELR")

melr_obs_cv=melr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
melr_obs_mean_cv=mean(melr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
melr_perm_mean_cvs=replicate(n_perm, {
  melr$tl_perm=sample(melr$tl)
  melr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
melr_perm_mean=mean(melr_perm_mean_cvs)
melr_perm_sd=sd(melr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
melr_obs_z=(melr_obs_mean_cv-melr_perm_mean)/melr_perm_sd
melr_perm_z=(melr_perm_mean_cvs-melr_perm_mean)/melr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(melr_p=mean(melr_perm_z <= melr_obs_z))

## Plot mean CVs and observed CV
melr_plot=ggplot()+
  geom_histogram(aes(x=melr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#E15759")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression(""))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression(""))+
  geom_vline(aes(xintercept=melr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Melvern"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)==0.257),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Milford
milr=filter(fish,impd=="MILR")

milr_obs_cv=milr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
milr_obs_mean_cv=mean(milr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
milr_perm_mean_cvs=replicate(n_perm, {
  milr$tl_perm=sample(milr$tl)
  milr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
milr_perm_mean=mean(milr_perm_mean_cvs)
milr_perm_sd=sd(milr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
milr_obs_z=(milr_obs_mean_cv-milr_perm_mean)/milr_perm_sd
milr_perm_z=(milr_perm_mean_cvs-milr_perm_mean)/milr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(milr_p=mean(milr_perm_z <= milr_obs_z))

## Plot mean CVs and observed CV
milr_plot=ggplot()+
  geom_histogram(aes(x=milr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#76B7B2")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression("Z-score mean length CV"))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression(""))+
  geom_vline(aes(xintercept=milr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Milford"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)==0.042),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Perry
perr=filter(fish,impd=="PERR")

perr_obs_cv=perr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
perr_obs_mean_cv=mean(perr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
perr_perm_mean_cvs=replicate(n_perm, {
  perr$tl_perm=sample(perr$tl)
  perr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
perr_perm_mean=mean(perr_perm_mean_cvs)
perr_perm_sd=sd(perr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
perr_obs_z=(perr_obs_mean_cv-perr_perm_mean)/perr_perm_sd
perr_perm_z=(perr_perm_mean_cvs-perr_perm_mean)/perr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(perr_p=mean(perr_perm_z <= perr_obs_z))

## Plot mean CVs and observed CV
perr_plot=ggplot()+
  geom_histogram(aes(x=perr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#59A14F")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression("Z-score mean length CV"))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression("Proportion of permutations"))+
  geom_vline(aes(xintercept=perr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Perry"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)==0.733),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Tuttle Creek
tcrr=filter(fish,impd=="TCRR")

tcrr_obs_cv=tcrr%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
tcrr_obs_mean_cv=mean(tcrr_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
tcrr_perm_mean_cvs=replicate(n_perm, {
  tcrr$tl_perm=sample(tcrr$tl)
  tcrr%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
tcrr_perm_mean=mean(tcrr_perm_mean_cvs)
tcrr_perm_sd=sd(tcrr_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
tcrr_obs_z=(tcrr_obs_mean_cv-tcrr_perm_mean)/tcrr_perm_sd
tcrr_perm_z=(tcrr_perm_mean_cvs-tcrr_perm_mean)/tcrr_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(tcrr_p=mean(tcrr_perm_z <= tcrr_obs_z))

## Plot mean CVs and observed CV
tcrr_plot=ggplot()+
  geom_histogram(aes(x=tcrr_perm_z,y=after_stat(density)),binwidth=0.2,fill="#EDC948")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression("Z-score mean length CV"))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression(""))+
  geom_vline(aes(xintercept=tcrr_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Tuttle Creek"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)<0.001),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Wolf Creek
wlfc=filter(fish,impd=="WLFC")

wlfc_obs_cv=wlfc%>%
  group_by(grid)%>%
  summarise(cv=sd(tl)/mean(tl))%>%
  ungroup()

## Observed mean CV across all grids
wlfc_obs_mean_cv=mean(wlfc_obs_cv$cv,na.rm=TRUE)

## Run permutation test
n_perm=1000
wlfc_perm_mean_cvs=replicate(n_perm, {
  wlfc$tl_perm=sample(wlfc$tl)
  wlfc%>%
    group_by(grid)%>%
    summarise(cv=sd(tl_perm)/mean(tl_perm))%>%
    ungroup()%>%
    summarise(mean_cv=mean(cv,na.rm=TRUE))%>%
    pull(mean_cv)
})

## Mean and standard deviation of mean CVs
wlfc_perm_mean=mean(wlfc_perm_mean_cvs)
wlfc_perm_sd=sd(wlfc_perm_mean_cvs)

## Z-score standardized obs_mean_cv and perm_mean_cvs
wlfc_obs_z=(wlfc_obs_mean_cv-wlfc_perm_mean)/wlfc_perm_sd
wlfc_perm_z=(wlfc_perm_mean_cvs-wlfc_perm_mean)/wlfc_perm_sd

## Calculate p-value: proportion of permutations with mean CV ≤ observed
(wlfc_p=mean(wlfc_perm_z <= wlfc_obs_z))

## Plot mean CVs and observed CV
wlfc_plot=ggplot()+
  geom_histogram(aes(x=wlfc_perm_z,y=after_stat(density)),binwidth=0.2,fill="#B07AA1")+
  scale_x_continuous(breaks=seq(-4,4,1),
                     name=expression("Z-score mean length CV"))+
  scale_y_continuous(breaks=seq(0,0.5,0.05),
                     name=expression(""))+
  geom_vline(aes(xintercept=wlfc_obs_z),linetype="dashed",linewidth=1.5)+
  annotate("text",label=expression("Wolf Creek"),x=4,y=0.5,hjust=1,vjust=1,size=8)+
  annotate("text",label=expression(italic(P)<0.001),x=4,y=0.46,hjust=1,vjust=1,size=6)+
  coord_cartesian(xlim=c(-4.05,4.05),
                  ylim=c(0,0.51),
                  expand=F)+
  pubtheme

################################################################################
## Combine and export plot
cvplot=wrap_plots(cltr_plot,eldr_plot,melr_plot,milr_plot,
                  perr_plot,tcrr_plot,wlfc_plot,plot_spacer(),
                  ncol=4)

ggsave(plot=cvplot,"cv.png",width=18,height=10,units="in",bg="white")