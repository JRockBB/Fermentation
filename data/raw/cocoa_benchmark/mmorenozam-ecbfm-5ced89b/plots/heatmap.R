#######################################################################
#                                                                     #
# R Script for plotting Figure S30 in supplementary material of:      #
#                                                                     #
# "Exploring cocoa bean fermentation processes by kinetic modelling"  #
#                                                                     #
# Function getSimsTable written by Charles Margossian in "Differen-   #                                                             #
# tial Equation Based Models in Stan", doi:10.5281/zenodo.1284264     #
#                                                                     #
# Code made by Mauricio Moreno-Zambrano                               #
# 2021                                                                #
#                                                                     #
#######################################################################

#loading necessary R packages

library(tidyverse)
library(latex2exp)
library(cowplot)

setwd("/home")      #set working directory accordingly
dpath<- '/home/'    #character string with the directory where Stan outputs are saved
outpth <- '/home/'  #character string with the directory where plots are going to be saved  

########################################################################
# FIGURE S30 IN THE SUPPLEMENTARY MATERIAL                             #
#                                                                      #
# explanation of functions:                                            #
#                                                                      #
# pmt_fun = function that returns a dataframe with re-scaled parameter #
#           estimates for kinetic parameters compared. Arguments:      #
#           fnam = name of Stan output                                 #
#           ivv = vector with maximum values for each time series      #
#           dsn = dataset label                                        #
#                                                                      #
# lt = function that returns a list of output names of a given model   #
#      Arguments:                                                      #
#      patt = pattern of model outputs to be used                      #
#                                                                      #
# Accesory functions that do not require input arguments directly by   #
# the user: getSimsTable.                                              #
#                                                                      #
########################################################################

getSimsTable <- function(x, ...){
  require(dplyr)
  nChains <- dim(x)[2]
  nPost <- dim(x)[1]
  x %>%
    as.data.frame(...) %>%
    mutate(chain = rep(1:nChains, ea = nPost),
           iteration = rep(1:nPost, nChains))
}

pmt_fun <- function(fnam,ivv,dsn){
  #parameter names in latex
  mu1 <- '$\\mu_{max}^{Y_{Glc}}$'
  mu2 <- '$\\mu_{max}^{Y_{Fru}}$'
  mu3 <- '$\\mu_{max}^{LAB_{Glc}}$'
  mu4 <- '$\\mu_{max}^{AAB_{EtOH}}$'
  mu5 <- '$\\mu_{max}^{AAB_{LA}}$'
  ks1 <- '$K_{Glc}^{Y}$'
  ks2 <- '$K_{Fru}^{Y}$'
  ks3 <- '$K_{Glc}^{LAB}$'
  ks4 <- '$K_{EtOH}^{AAB}$'
  ks5 <- '$K_{LA}^{AAB}$'
  k1  <- '$k_{Y}$'
  k2  <- '$k_{LAB}$'
  k3  <- '$k_{AAB}$'
  y1  <- '$Y_{Glc|Y}$'
  y2  <- '$Y_{Glc|LAB}$'
  y3  <- '$Y_{Fru|Y}$'
  y4  <- '$Y^{Glc}_{EtOH|Y}$'
  y5  <- '$Y^{Fru}_{EtOH|Y}$'
  y6  <- '$Y_{EtOH|AAB}$'
  y7  <- '$Y^{Glc}_{LA|LAB}$'
  y8  <- '$Y_{LA|AAB}$'
  y9  <- '$Y^{Glc}_{Ac|LAB}$'
  y10 <- '$Y^{EtOH}_{Ac|AAB}$'
  y11 <- '$Y^{LA}_{Ac|AAB}$'
  #parn <- c(mu1,mu2,mu3,mu4,mu5,ks1,ks2,ks3,ks4,ks5,k1,k2,k3,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)
  parn <- c(y6,y8,y11)
  rnams <- c('0','1','2','3','4','1,2','1,3','1,4','2,3','2,4','3,4','1,2,3','1,2,4',
             '1,3,4','2,3,4','1,2,3,4')
  dummy_df <- data.frame(cbind(rnams,gsub(',','',paste0('mm',rnams))))
  dummy_df$rnams <- factor(dummy_df$rnams,levels=dummy_df$rnams)
  mlo <- list()
  for (i in 1:length(fnam)){
    load(paste0(dpath,fnam[i]))
    mnam <- gsub('.Rsave','',fnam[i])
    mnam <- sub('.*_','',mnam)
    df <- getSimsTable(fit)
    df$yc6 = df$yc6*(ivv[3]/ivv[8]) #ethanol
    df$yc8 = df$yc8*(ivv[4]/ivv[8]) #lactic acid
    df$yc11 = df$yc11*(ivv[5]/ivv[8]) #acetic acid

    mdf = c(round(mean(df$yc6),6),
      round(mean(df$yc8),6),
      round(mean(df$yc11),6))
    
    dfout <- data.frame(cbind(mnam,dsn,parn,mdf))
    mlo[[i]] <- dfout
  }
  dfout <- tibble(do.call('rbind',mlo))
  dfout <- dfout %>% right_join(.,dummy_df,by=c('mnam'='V2'))
  dfout <- dfout[complete.cases(dfout),]
  return(dfout)
}

lt <- function(patt){
  lmm <- list.files(dpath,pattern=patt)
  return(lmm)
}


#maximum values of each dataset required to re-scaled posterior estimates to their real units
ghhp1 <- c(51.963000000,57.741000000,22.492000000,8.593000000,6.139000000,0.955193281,0.566121975,0.001914953)
ecwb1 <- c(42.083000000,53.899000000,17.847000000,5.062000000,7.842000000,0.841571964,0.647008540,0.002723692)
ecwb2 <- c(36.674000000,32.884000000,12.517000000,5.152000000,24.097000000,0.662355671,0.641076730,0.005459565)
brpb1 <- c(1.023420e+02,9.505800e+01,6.532600e+01,3.142500e+01,6.872600e+01,5.483922e-01,3.356681e-01,6.041680e-04)

# values for ticks on x axis for Figure S30:
rnams <- c('0','1','2','3','4','1,2','1,3','1,4','2,3','2,4','3,4','1,2,3','1,2,4',
           '1,3,4','2,3,4','1,2,3,4')

# extracting re-scaled parameters for datasets ghhp1, ecwb1, ecwb2 and brpb1
ghhp1_df <- pmt_fun(lt('ghhp1'),ghhp1,'ghhp1')
ecwb1_df <- pmt_fun(lt('ecwb1'),ecwb1,'ecwb1')
ecwb2_df <- pmt_fun(lt('ecwb2'),ecwb2,'ecwb2')
brpb1_df <- pmt_fun(lt('brpb1'),brpb1,'brpb1')

# plotting panels for FigureS30: p1<-(a), p2<-(b), p3<-(c), p4<-(d)
p1<-ggplot(data=ghhp1_df,aes(y=parn,x=rnams,fill=(as.numeric(mdf))))+
  geom_tile()+
  scale_fill_gradient(low='white',high='red',na.value=NA,name=expression(paste('mg g(pulp) '^{-1})))+
  scale_y_discrete(labels=c('$Y_{EtOH|AAB}$'=parse(text=TeX('$Y_{EtOH|AAB}$')),
                            '$Y_{LA|AAB}$'=parse(text=TeX('$Y_{LA|AAB}$')),
                            '$Y^{LA}_{Ac|AAB}$'=parse(text=TeX('$Y^{LA}_{Ac|AAB}$'))))+
  scale_x_discrete(limits=rnams)+
  theme_minimal()+ 
  theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank()
  )+
  ylab(NULL)+
  coord_fixed(clip='on')+
  xlab(NULL)+
  labs(title = "(a)")

p2<-ggplot(data=ecwb1_df ,aes(y=parn,x=rnams,fill=(as.numeric(mdf))))+
  geom_tile()+
  scale_fill_gradient(low='white',high='red',na.value=NA,name=expression(paste('mg g(pulp) '^{-1})))+
  scale_y_discrete(labels=c('$Y_{EtOH|AAB}$'=parse(text=TeX('$Y_{EtOH|AAB}$')),
                            '$Y_{LA|AAB}$'=parse(text=TeX('$Y_{LA|AAB}$')),
                            '$Y^{LA}_{Ac|AAB}$'=parse(text=TeX('$Y^{LA}_{Ac|AAB}$'))))+
  scale_x_discrete(limits=rnams)+
  theme_minimal()+ 
  theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank()
  )+
  ylab(NULL)+
  coord_fixed(clip='on')+
  xlab(NULL)+
  labs(title = "(b)")

p3<-ggplot(data=ecwb2_df ,aes(y=parn,x=rnams,fill=(as.numeric(mdf))))+
  geom_tile()+
  scale_fill_gradient(low='white',high='red',na.value=NA,name=expression(paste('mg g(pulp) '^{-1})))+
  scale_y_discrete(labels=c('$Y_{EtOH|AAB}$'=parse(text=TeX('$Y_{EtOH|AAB}$')),
                            '$Y_{LA|AAB}$'=parse(text=TeX('$Y_{LA|AAB}$')),
                            '$Y^{LA}_{Ac|AAB}$'=parse(text=TeX('$Y^{LA}_{Ac|AAB}$'))))+
  scale_x_discrete(limits=rnams)+
  theme_minimal()+ 
  theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_blank()
  )+
  ylab(NULL)+
  coord_fixed(clip='on')+
  xlab(NULL)+
  labs(title = "(c)")

p4<-ggplot(data=brpb1_df,aes(y=parn,x=rnams,fill=(as.numeric(mdf))))+
  geom_tile()+
  scale_fill_gradient(low='white',high='red',na.value=NA,name=expression(paste('mg g(pulp) '^{-1})))+
  scale_y_discrete(labels=c('$Y_{EtOH|AAB}$'=parse(text=TeX('$Y_{EtOH|AAB}$')),
                            '$Y_{LA|AAB}$'=parse(text=TeX('$Y_{LA|AAB}$')),
                            '$Y^{LA}_{Ac|AAB}$'=parse(text=TeX('$Y^{LA}_{Ac|AAB}$'))))+
  scale_x_discrete(limits=rnams)+
  theme_minimal()+ 
  theme(plot.title = element_text(size=18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=14),
        legend.title = element_text(size=18),
        axis.title=element_text(size=16),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90,vjust=0.2,hjust=0.95)
  )+
  ylab(NULL)+
  coord_fixed(clip='on')+
  xlab('MI ( )')+
  labs(title = "(d)")


plt <- plot_grid(p1,p2,p3,p4,ncol=1,align='hv',axis='bt')
save_plot(paste0(outpth,'FigureS30.pdf'),plt,base_height = 13,base_width =13)