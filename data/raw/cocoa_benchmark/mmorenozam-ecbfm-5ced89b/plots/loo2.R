#######################################################################
#                                                                     #
# R Script for plotting Figure S31 in supplementary material of:      #
#                                                                     #
# "Exploring cocoa bean fermentation processes by kinetic modelling"  #
#                                                                     #
# Code made by Mauricio Moreno-Zambrano                               #
# 2021                                                                #
#                                                                     #
#######################################################################

#loading necessary R packages

library(loo)
library(tidyverse)
library(cowplot)

setwd("/home")      #set working directory accordingly
dpath<- '/home/'    #character string with the directory where Stan outputs are saved
outpth <- '/home/'  #character string with the directory where plots are going to be saved  

########################################################################
# FIGURE S31 IN THE SUPPLEMENTARY MATERIAL                             #
#                                                                      #
# explanation of functions:                                            #
#                                                                      #
# loo_ext = Function to extract psis-loo estimates and standard errors #
#           It returns a dataframe ready for plotting. Arguments:      #
#           patt = pattern for creating a list with all outputs of a   #
#                  dataset in particular                               #
#           patt1 = labelling for extracted estimated parameters       #
#                                                                      #
########################################################################


loo_ext<-function(patt,patt1){
  fname <- list.files(dpath,pattern=patt)
  looic <- c()
  looic_se <- c()
  mi <- c()
  for (i in 1:length(fname)){
    load(paste0(dpath,fname[i]))
    loo_stm <- loo(fit)
    looic[i] <- loo_stm$estimates[3]
    looic_se[i] <- loo_stm$estimates[6]
    mi[i] <- gsub('.Rsave','',gsub(patt,'',fname[i]))
  }
  out_df <- data.frame(cbind(round(looic,2),round(looic_se,2),mi))
  out_df$V1 <- as.numeric(out_df$V1)
  out_df$V2 <- as.numeric(out_df$V2)
  out_df$data <- patt1
  
  dat_lev <- c('brwb1','brwb2','mywb3')
  
  rnams <- c('0','1','2','3')
  
  dummy_df <- data.frame(cbind(rnams,gsub(',','',paste0('mm',rnams))))
  dummy_df$rnams <- factor(dummy_df$rnams,levels=dummy_df$rnams)
  out_df<-  out_df %>% right_join(.,dummy_df,by=c('mi'='V2'))  
  out_df <- out_df[complete.cases(out_df),]
  out_df$data <- factor(out_df$data,levels=dat_lev)
  
  return(out_df)
}


#creating a dataframe with all psis-loo values from all fitted datasets

alles <- rbind(loo_ext("ghhp1_",'ghhp1'),                             
               loo_ext("brwb1_",'brwb1'),
               loo_ext("brwb2_",'brwb2'))

#plotting Figure S31

l_plt<-ggplot(alles,aes(y=V1,x=rnams))+
  geom_point(color="#D43F3AFF")+
  geom_errorbar(aes(ymin=V1-V2,ymax=V1+V2),color='#EEA2367F',size=0.4)+
  facet_wrap(~data,scales='free_y',ncol=3)+
  xlab("MI( )")+
  ylab('PSIS-LOO')+
  theme_bw()+
  theme(strip.text.x = element_text(hjust = -0.01),
        panel.grid.major = element_line(colour = "grey90",linetype='dashed',size=0.2),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=16),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12),
        axis.text.x = element_text(angle = 90,vjust=0.2,hjust=0.95)
  )

save_plot('psis_loo_temp.pdf',l_plt,base_height = 14.5,base_width = 12)