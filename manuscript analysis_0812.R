#--------------------------------------------------------------
# Ben Neely
# 08/12/2025
# Assess number of sites needed for Blue Catfish CPE precision and size structure
#--------------------------------------------------------------

## Clear R
cat("\014")  
rm(list=ls())

## Install and load packages
if("FSA" %in% rownames(installed.packages()) == FALSE) {install.packages("FSA")}
library(FSA)
if("rio" %in% rownames(installed.packages()) == FALSE) {install.packages("rio")}
library(rio)
if("patchwork" %in% rownames(installed.packages()) == FALSE) {install.packages("patchwork")}
library(patchwork)
if("tidyverse" %in% rownames(installed.packages()) == FALSE) {install.packages("tidyverse")}
library(tidyverse)
if("future.apply" %in% rownames(installed.packages()) == FALSE) {install.packages("future.apply")}
library(future.apply)

## Set up parallel processing
plan(sequential)

## Set seed for reproducibility
RNGkind(kind="Mersenne-Twister", normal.kind="Inversion", sample.kind="Rejection")
set.seed(812)

## Set ggplot theme
pubtheme=theme_classic()+
  theme(panel.grid=element_blank(), 
        panel.background=element_blank(),
        plot.background=element_blank(),
        panel.border=element_rect(fill="transparent"),
        axis.title=element_text(size=22, color="black", face="bold"),
        axis.text=element_text(size=18, color="black"),
        legend.position="inside",
        legend.position.inside=c(0.99, 0.99),
        legend.justification=c("right", "top"),
        legend.text=element_text(size=20),
        legend.title=element_blank())
options(scipen=999)

## Read in data
samp=import("bcf sampling data.xlsx", which="samp")%>%
  mutate(grid=as.numeric(grid))
fish=import("bcf sampling data.xlsx", which="fish")%>%
  expandCounts(~count)%>%
  mutate(grid=as.numeric(grid))

################################################################################
## Objective 1: Summary statistics for table
################################################################################
fish1=fish%>%
  mutate(tot=case_when(tl>=1 ~ 1, TRUE ~ 0),
         stock=case_when(tl>=300 ~ 1, TRUE ~ 0),
         qual=case_when(tl>=510 ~ 1, TRUE ~ 0),
         over=case_when(tl>=760 ~ 1, TRUE ~ 0))%>%
  group_by(impd, grid)%>%
  summarize(tot_catch=sum(tot),
            st_catch=sum(stock),
            qual_catch=sum(qual),
            over_catch=sum(over))%>%
  ungroup()

## Merge fish and samples for CPE
out1=samp%>%
  left_join(fish1,by=c("impd","grid"))%>%
  replace_na(list(tot_catch=0,st_catch=0,qual_catch=0,over_catch=0))%>%
  mutate(tot_cpe=tot_catch/eff,
         st_cpe=st_catch/eff,
         over_cpe=over_catch/eff)%>%
  select(impd,grid,tot_catch,tot_cpe,st_catch,st_cpe,qual_catch,over_catch,over_cpe)

sstats=out1%>%
  group_by(impd)%>%
  summarize(n=n(),
            tot_n1=sum(tot_catch),
            tot_cpe1=mean(tot_cpe),
            tot_sd1=sd(tot_cpe),
            tot_se1=tot_sd1/sqrt(n),
            tot_rse1=(tot_se1/tot_cpe1)*100,
            st_n1=sum(st_catch),
            st_cpe1=mean(st_cpe),
            st_sd1=sd(st_cpe),
            st_se1=st_sd1/sqrt(n),
            st_rse1=(st_se1/st_cpe1)*100,
            qu_n1=sum(qual_catch),
            over_n1=sum(over_catch),
            over_cpe1=mean(over_cpe),
            over_sd1=sd(over_cpe),
            over_se1=over_sd1/sqrt(n),
            over_rse1=(over_se1/over_cpe1)*100)%>%
  ungroup()%>%
  mutate(psdq=(qu_n1/st_n1)*100)%>%
  select(impd,n,tot_n1,tot_cpe1,tot_rse1,over_n1,over_cpe1,over_rse1,psdq)

export(sstats, "summary stats.xlsx")

################################################################################
## View data for reporting in results
sstats

## Effort
sum(sstats$n)
sum(sstats$n)*0.08333
mean(sstats$n)
sd(sstats$n)

## Total fish
sum(sstats$tot_n1)
sum(sstats$over_n1)

## CPE
mean(sstats$tot_cpe1)
sd(sstats$tot_cpe1)
mean(sstats$over_cpe1)
sd(sstats$over_cpe1)

################################################################################
## Objective 2: Precision of CPEtot and CPE760
################################################################################
fish2=fish%>%
  mutate(tot=case_when(tl>=1 ~ 1, TRUE ~ 0),
         over=case_when(tl>=760 ~ 1, TRUE ~ 0))%>%
  group_by(impd,grid)%>%
  summarize(tot_catch=sum(tot),
            over_catch=sum(over))%>%
  ungroup()

out2=samp%>%
  left_join(fish2,by=c("impd","grid"))%>%
  replace_na(list(tot_catch=0,over_catch=0))%>%
  mutate(tot_cpe=tot_catch/eff,
         over_cpe=over_catch/eff)%>%
  select(impd,grid,tot_catch,tot_cpe,over_catch,over_cpe)

# Function to calculate RSE
calc_rse=function(x) {
  mu=mean(x)
  n=length(x)
  if(mu==0 || n<2) return(Inf)
  se=sd(x)/sqrt(n)
  return(se/mu)
}

# Function to calculate RSE percentiles
rsefunc=function(data,reps=1000,sample_sizes=seq(5,200,by=5)) {
  set.seed(812)
  data=na.omit(data)
  if(length(data)<2) return(NULL)
  
  res_list=lapply(sample_sizes, function(n) {
    rse_vals=replicate(reps, {
      sample_data=sample(data, size=min(n, length(data)), replace=TRUE)
      calc_rse(sample_data)
    })
    tibble(n=n,
           rse_70=quantile(rse_vals, 0.70, na.rm=TRUE),
           rse_80=quantile(rse_vals, 0.80, na.rm=TRUE),
           rse_90=quantile(rse_vals, 0.90, na.rm=TRUE))
  })
  
  bind_rows(res_list)
}

# Apply to CPEtot
rse_tot=out2%>%
  group_by(impd)%>%
  group_map(~{
    rsefunc(.x$tot_cpe)%>%
      select(n,rse_70,rse_80,rse_90)%>%
      mutate(impd=.y$impd)
  })%>%
  bind_rows()%>%
  mutate(across(starts_with("rse_"), ~ ifelse(is.infinite(.), NA, .)))

# Apply to CPE760
rse_760=out2%>%
  group_by(impd)%>%
  group_map(~{
    rsefunc(.x$over_cpe)%>%
      select(n,rse_70,rse_80,rse_90)%>%
      mutate(impd=.y$impd)
  })%>%
  bind_rows()%>%
  mutate(across(starts_with("rse_"), ~ ifelse(is.infinite(.), NA, .)))

# Plot for CPEtot
cpetot_plot=ggplot(rse_tot)+
  geom_line(aes(x=n,y=rse_80,color=impd),linewidth=1.5)+
  geom_ribbon(aes(x=n,ymin=rse_70,ymax=rse_90,fill=impd),alpha=0.3)+
  geom_hline(yintercept=0.25, linetype="dashed", color="black")+
  scale_color_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                     values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_fill_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                    values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_y_continuous(breaks=seq(0,1,0.1),
                     labels=scales::label_percent(accuracy=1),
                     name=expression("RSE of "*CPE[tot]))+
  scale_x_continuous(breaks=seq(0,200,20),
                     name="")+
  coord_cartesian(xlim=c(0,202),ylim=c(0,1.01),expand=F)+
  pubtheme

# Plot for CPE760
cpe760_plot=ggplot(rse_760)+
  geom_line(aes(x=n, y=rse_80, color=impd), linewidth=1.5)+
  geom_ribbon(aes(x=n, ymin=rse_70, ymax=rse_90, fill=impd), alpha=0.3)+
  geom_hline(yintercept=0.25, linetype="dashed", color="black")+
  scale_color_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                     values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_fill_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                    values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     labels=scales::label_percent(accuracy=1),
                     name=expression("RSE of "*CPE[760]))+
  scale_x_continuous(breaks=seq(0, 200, 20),
                     name=expression("Sample sites"))+
  coord_cartesian(xlim=c(0, 202), ylim=c(0, 1.01), expand=F)+
  pubtheme+
  theme(legend.position="none")

# Combine and export plots
rseplot=cpetot_plot/cpe760_plot
ggsave(plot=rseplot,"Figure 2.png",width=10,height=10,units="in",bg="white")

################################################################################
## Find where RSE <= 0.25 for each population
rse_tot%>%
  filter(rse_80<=0.25)%>%
  group_by(impd)%>%
  slice_min(n,n=1)%>%
  select(impd,n,rse_80)

rse_760%>%
  filter(rse_80<=0.25)%>%
  group_by(impd)%>%
  slice_min(n,n=1)%>%
  select(impd,n,rse_80)

## Find RSE760 at 200 sites
rse_760%>%
  filter(n==200)%>%
  select(impd,n,rse_80)
  
################################################################################
## Objective 3: Sample size needed for representative size structure
################################################################################
fish$site_id=paste(fish$impd,fish$grid,sep="_")
samp$site_id=paste(samp$impd,samp$grid,sep="_")

all_sites=unique(samp[, c("impd", "grid", "site_id")])
lengths_by_site=split(fish$tl, fish$site_id)
missing_sites=setdiff(all_sites$site_id, names(lengths_by_site))
for(ms in missing_sites) lengths_by_site[[ms]]=numeric(0)
lengths_by_site=lengths_by_site[all_sites$site_id]
all_sites$lengths=lengths_by_site
impd_lists=split(all_sites, all_sites$impd)
full_lengths_by_impd=lapply(impd_lists, function(df) unlist(df$lengths))

# Function to run resampling comparison
run_resample_sim=function(site_df, n_sites, n_reps=1000) {
  site_lengths=site_df$lengths
  full_dist=unlist(site_lengths)
  impd=unique(site_df$impd)
  
  p_vals=future_replicate(n_reps, {
    sample_sites=sample(site_lengths, min(n_sites, length(site_lengths)), replace=FALSE)
    sample_dist=unlist(sample_sites)
    if(length(sample_dist) < 5) return(0)  # p = 0 if too few fish
    suppressWarnings(ks.test(sample_dist, full_dist)$p.value)
  },future.seed=812)
  
  successes=sum(p_vals > 0.05, na.rm=TRUE)
  ci=binom.test(successes, n_reps, conf.level=0.95)$conf.int
  prop_p_gt_0.05=mean(p_vals > 0.05, na.rm=TRUE)
  
  data.frame(impd=impd, n_sites=n_sites, prop_p_gt_0.05=prop_p_gt_0.05, lci=ci[1], uci=ci[2])
}

# Apply to all impoundments
sample_sizes=seq(5, 200, 5)
out1=do.call(rbind, lapply(names(impd_lists), function(imp) {
  do.call(rbind, lapply(sample_sizes, function(n) {
    run_resample_sim(impd_lists[[imp]], n)
  }))
}))

# Plot with confidence intervals
sizestructure_plot=ggplot(out1)+
  geom_line(aes(x=n_sites, y=prop_p_gt_0.05, color=impd), linewidth=1.5)+
  geom_ribbon(aes(x=n_sites, ymin=lci, ymax=uci, fill=impd), alpha=0.3)+
  geom_hline(yintercept=0.80, linetype="dashed", color="black")+
  scale_color_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                     values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_fill_manual(labels=c("Clinton", "El Dorado", "Melvern", "Milford", "Perry", "Tuttle Creek", "Wolf Creek"),
                    values=c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948", "#B07AA1"))+
  scale_x_continuous(breaks=seq(0, 200, 20),
                     name=expression("Sample sites"))+
  scale_y_continuous(breaks=seq(0, 1, 0.1),
                     name = expression("Pr("*italic(P)*" ≥ 0.05 from KS test)"))+
  coord_cartesian(xlim=c(0, 205), 
                  ylim=c(0, 1.01), 
                  expand=F)+
  pubtheme+
  theme(legend.position=c(0.99,0.01),
        legend.justification=c("right","bottom"))

# Export plot
ggsave(plot=sizestructure_plot, "Figure 3.png", width=7, height=4.5, units="in", bg="white")

################################################################################
## Find where size structure is adequately estimated for each population
out1%>%
  filter(prop_p_gt_0.05>=0.80)%>%
  group_by(impd)%>%
  slice_min(n_sites,n=1)%>%
  select(impd,n_sites,prop_p_gt_0.05)

## Clean up
plan(sequential)