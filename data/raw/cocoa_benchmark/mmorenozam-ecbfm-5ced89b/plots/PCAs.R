#######################################################################
#                                                                     #
# R Script for plotting Figures S27, S28 and S29 in supplementary ma- #
# terial and Figures 3 and 4 in main manuscript of:                   #
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

library(dplyr)
library(HDMD)
library(tidyverse)
library(cowplot)

setwd("/home")      #set working directory accordingly
dpath<- '/home/'    #character string with the directory where Stan outputs are saved
outpth <- '/home/'  #character string with the directory where plots are going to be saved  

########################################################################
# FIGURE S27 IN THE SUPPLEMENTARY MATERIAL                             #
#                                                                      #
# explanation of functions:                                            #
#                                                                      #
# dg_list_all, dg_list_kss, dg_list_mus, dg_list_ycf, dg_list_yest,    #
# dg_list_lab, dg_list_aab = Series of functions to create lists of    #
#                            dataframes with output from function      #
#                            sfunc. Function dg_list_all needs to be   #
#                            called first in order for the rest to     #
#                            work. Argument in common to all:          #
#                            ncods = string character vector with mo-  #
#                                    del codes                         #
#                                                                      #
# pcag_fun = Function that returns a dataframe with heatmap results of #
#            Figure S27 ready to plot. Arguments:                      #
#            ncod = vector with numeric part of models' code           #
#            ls_rd = output from dg_list_ functions                    #
#            feat = feature code name to be object of PCA              #
#                   'cult' = Cultivar                                  #
#                   'meth' = Fermentation method                       #
#                   'country' = Country of origin                      #
#                   'temp' = Controlled temperature                    #
#                   'turning' = Turning of fermenting mass             #
#            vars_group = label for subgroup of parameters subjet to   #
#                         PCA                                          #
#                                                                      #
# plt_fun = Plotting function for dataframe returned by pcag_fun. Arg- #
#           uments:                                                    #
#           df = dataframe resulted from pcag_fun                      #
#                                                                      #
# Accesory functions that do not require input arguments directly by   #
# the user: getSimsTable, sfunc, gn_fun, pca.md.                       #
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

sfunc <- function(lnames,patk,vdc){
  pttn <- paste0('_',patk,'.Rsave')
  dli <- list()
  for (i in 1:length(lnames)){
    load(paste0(dpath,lnames[i]))
    prd <- getSimsTable(fit)
    ch <- prd$chain
    prd <- prd[,!grepl("^x",colnames(prd))]
    prd <- prd[,!grepl("^l",colnames(prd))]
    prd <- prd[,!grepl("^ks",colnames(prd))]
    prd <- prd[,c(1:(ncol(prd)-2))]
    lnames[i] <- str_remove(lnames[i],pttn)
    prd$chain <- ch
    prd$data <- lnames[i]
    prd <- prd %>%  mutate(data = recode(data, #aut,trialcode,country,cultivar,method,starter,turn,temp,simul
                                         `brwb1` = 'papa_pbx1_BR_cf_wb_st0_tr_nctrl_si0',
                                         `brwb2` = 'papa_pbx2_BR_cf_wb_st0_tr_nctrl_si0',
                                         `ghhp1` = 'camu_h5_GH_cf_hp_st0_ntr_nctrl_si0',
                                         `ghhp2` = 'camu_h10_GH_NA_hp_st0_tr_nctrl_si0',
                                         `ghhp3` = 'camu_h11_GH_NA_hp_st0_ntr_nctrl_si0',
                                         `ghhp4` = 'camu_h12_GH_NA_hp_st0_tr_nctrl_si0',
                                         `ghhp5` = 'camu_h13_GH_NA_hp_st0_ntr_nctrl_si0',
                                         `ecpb1` = 'lee_lpc_EC_c_pc_st0_NA_ctrl_si1',
                                         `brpb1` = 'melo_mpc_BR_NA_pc_st0_tr_ctrl_si0',
                                         `brst1` = 'melo_mst12_BR_NA_st_st0_tr_ctrl_si0',
                                         `brwb3` = 'melo_mbx1_BR_NA_wb_st0_tr_nctrl_si0',
                                         `brwb4` = 'melo_mbx2_BR_NA_wb_st0_tr_nctrl_si0',
                                         `brst2` = 'melo_mst13_BR_NA_st_st0_tr_nctrl_si0',
                                         `ecwb1` = 'papa_pbx1ec_EC_nt_wb_st0_tr_nctrl_si0',
                                         `ecwb2` = 'papa_pbx2ec_EC_nt_wb_st0_tr_nctrl_si0',
                                         `ecpt1` = 'papa_pt1_EC_nt_pt_st0_ntr_nctrl_si0',
                                         `ecpt2` = 'papa_pt2_EC_nt_pt_st0_ntr_nctrl_si0',
                                         `mywb3` = 'papa_pbx13_MY_NA_wb_st0_tr_nctrl_si0',
                                         `niwb1` = 'papa_ngu_NI_nugu_wb_st0_tr_nctrl_si0',
                                         `niwb2` = 'papa_opa_NI_opayo_wb_st0_tr_nctrl_si0',
                                         `brwb7` = 'veiga_ph16_BR_ph16_wb_st0_NA_nctrl_si0',
                                         `dowb1` = 'lagunes_lbx_DO_t_wb_st0_tr_nctrl_si0',
                                         `hnwb1` = 'romanens_rbx_HN_IUU_wb_st0_tr_nctrl_si0'))
    
    prd$sigma <- NULL
    
    if (!is.null(vdc)){
      prd <- dplyr::select(prd,c(all_of(vdc),'chain','data'))
    }
    
    df <- prd
    set.seed(100785)
    df <- df %>%
      gather(par,value,-data,-chain)%>%
      group_by(par,data,chain)%>%
      mutate(ql=quantile(value,0.025,na.rm=T),qu=quantile(value,0.975,na.rm=T),row=row_number()) %>%
      dplyr::filter(value <= qu & value >=ql )%>%
      dplyr::select(-qu,-ql)%>%
      sample_n(1900) %>%
      mutate(row=seq(1,1900,1)) %>%
      spread(par,value) %>%
      dplyr::select(-row)%>% 
      group_by(data,chain) %>% 
      separate(data,into=c('aut','ident','country','cult','meth','start','turn','temp','simul'),sep='_')
    
    dli[[i]] <- df
  }
  
  dli <- do.call('rbind',dli)
  
  return(dli)
}

gn_fun <- function(patt){#function to read bulks of files per model
  mm <- list.files(dpath,pattern=paste0(patt,'.Rsave'))
  return(mm)
}

dg_list_all <- function(ncods){#function to creat dataframes with all parameters from a given list of models
  ls_rd <- list()
  for (i in 1:length(ncods)){
    ls_rd[[i]] <- sfunc(gn_fun(ncods[i]),ncods[i],NULL)
  }
  return(ls_rd)
}

dg_list_kss <- function(ncods){#function to creat dataframes with mortality parameters from a given dataframe
  ls_rd <- list()
  for (i in 1:length(ncods)){
    ls_rd[[i]] <- tibble(cbind(ncods[[i]][,c(1:10)],ncods[[i]][,grepl("^k",colnames(ncods[[i]]))]))
  }
  return(ls_rd)
}

dg_list_mus <- function(ncods){#function to creat dataframes with growth parameters from a given dataframe
  ls_rd <- list()
  for (i in 1:length(ncods)){
    ls_rd[[i]] <- tibble(cbind(ncods[[i]][,c(1:10)],ncods[[i]][,grepl("^mu",colnames(ncods[[i]]))]))
  }
  return(ls_rd)
}

dg_list_ycf <- function(ncods){#function to creat dataframes with yield coefficients from a given dataframe
  ls_rd <- list()
  for (i in 1:length(ncods)){
    ls_rd[[i]] <- tibble(cbind(ncods[[i]][,c(1:10)],ncods[[i]][,grepl("^y",colnames(ncods[[i]]))]))
  }
  return(ls_rd)
}

dg_list_yst <- function(ncods){#function to creat dataframes with yeast related parameters from a given dataframe
  ls_rd <- list()
  prd <- list()
  
  for (i in 1:length(ncods)){
    vdc <- c('mu1','mu2','yc1','yc3','yc4','yc5','k1')
    if (i==1|i==2|i==3|i==6){#mm0,mm1,mm2,mm12
      vdc <- vdc
    }
    if (i==4|i==7|i==9|i==12){#mm3,mm13,mm23,mm123
      vdc <- c(vdc,'yc12','yc13') 
    }
    if (i==5|i==8|i==10|i==13){#mm4,mm14,mm24,mm124
      vdc <- c(vdc,'mu6','yc12','yc13')
    }
    if (i==11|i==14|i==15|i==16){#mm34,mm134.mm234,mm1234
      vdc <- c(vdc,'mu6','yc12','yc13','yc14','yc15')
    }
    prd[[i]] <- dplyr::select(ncods[[i]],'aut','ident','country','cult','meth','start',
                              'turn','temp','simul','chain',
                              c(all_of(vdc)))
  }
  return(prd)
}

dg_list_lab <- function(ncods){#function to creat dataframes with LAB related parameters from a given dataframe
  ls_rd <- list()
  prd <- list()
  
  for (i in 1:length(ncods)){
    vdc <- c('mu3','k2','yc2','yc7','yc9')
    if (i==1|i==2|i==4|i==5|i==7|i==8|i==11|i==14){#mm0,mm1,mm3,mm4,mm13,mm14.mm34,mm134
      vdc <- vdc
    }
    if (i==3|i==6|i==9|i==10|i==12|i==13|i==15|i==16){#mm2,mm12,mm23,mm24,mm123,mm124,mm234,mm1234
      vdc <- c(vdc,'mu6','yc12','yc13','yc14','yc15','yc16')
    }
    prd[[i]] <- dplyr::select(ncods[[i]],'aut','ident','country','cult','meth','start',
                              'turn','temp','simul','chain',
                              c(all_of(vdc)))
  }
  return(prd)
}

dg_list_aab <- function(ncods){#function to creat dataframes with AAB related parameters from a given dataframe
  ls_rd <- list()
  prd <- list()
  for (i in 1:length(ncods)){
    vdc <- c('mu4','mu5','k3','yc6','yc8','yc10','yc11')
    prd[[i]] <- dplyr::select(ncods[[i]],'aut','ident','country','cult','meth','start',
                              'turn','temp','simul','chain',
                              c(all_of(vdc)))
  }
  return(prd)
}

pca.md <- function(df,feat,mm){ #function to compute Malahanobis distance
  fl <- names(which(table(df[,feat])==4))
  if (length(fl)!=0){
    for (j in 1:length(fl)){
      df <- subset(df,get(feat)!=fl[j])
    }
  }
  
  if (feat=='cult'){
    df <- subset(df,cult!='NA')
  }
  
  if (feat=='turn'){
    df <- subset(df,turn!='NA')
  }
  
  if (feat=='meth'){
    df <- subset(df,meth!='NA')
  }
  
  if (feat=='country'){
    df <- subset(df,country!='NA')
  }
  
  if (feat=='temp'){
    df <- subset(df,temp!='NA')
  }
  
  my.pca <- prcomp(df[c(11:(ncol(df)))],retx=T,center=T,scale=F)
  pca.ind <- get_pca_ind(my.pca)
  df.sc <- cbind(pca.ind$coord[,c(1,2)],df[,feat])
  md <-  pairwise.mahalanobis(df.sc[,1:2],grouping=df.sc[,feat],cov=cov(df.sc[,1:2]))
  mdd <- md$distance
  fmd <- mdd[upper.tri(mdd,diag=F)]
  out <- cbind(round(median(fmd),2),feat,mm)
  return(out)
}

pcag_fun <- function(ncod,ls_rd,feat,vars_group){
  ls_md <- list()
  for (i in 1:length(ncod)){
    ls_md[[i]] <- pca.md(ls_rd[[i]],feat,ncod[i])
  }
  
  ls_md<- data.frame(do.call('rbind',ls_md))
  ls_md$group <- vars_group
  ls_md$V1 <- as.numeric(ls_md$V1)
  ls_md$mm <- factor(ls_md$mm,levels=ncod)
  ls_md$group <- factor(ls_md$group,levels=c('ALL','MSGR','MR','YC','Y-related','LAB-related','AAB-related'))
  return(ls_md)
}

ncod_pl <- c('0','1','2','3','4','1,2','1,3','1,4','2,3','2,4','3,4',
             '1,2,3','1,2,4','1,3,4','2,3,4','1,2,3,4')

ncodes <- c('mm0','mm1','mm2','mm3','mm4','mm12','mm13','mm14','mm23','mm24','mm34',
            'mm123','mm124','mm134','mm234','mm1234')

alles_raw <- dg_list_all(ncodes)
deded_raw <- dg_list_kss(alles_raw)
muses_raw <- dg_list_mus(alles_raw)
yield_raw <- dg_list_ycf(alles_raw)
yeast_raw <- dg_list_yst(alles_raw)
labbb_raw <- dg_list_lab(alles_raw)
aabbb_raw <- dg_list_aab(alles_raw)

alles_heat <- rbind(pcag_fun(ncod_pl,alles_raw,'cult','ALL'),
                    pcag_fun(ncod_pl,muses_raw,'cult','MSGR'),
                    pcag_fun(ncod_pl,yield_raw,'cult','YC'),
                    pcag_fun(ncod_pl,deded_raw,'cult','MR'),
                    pcag_fun(ncod_pl,yeast_raw,'cult','Y-related'),
                    pcag_fun(ncod_pl,labbb_raw,'cult','LAB-related'),
                    pcag_fun(ncod_pl,aabbb_raw,'cult','AAB-related'))

alles_heat1 <- rbind(pcag_fun(ncod_pl,alles_raw,'meth','ALL'),
                     pcag_fun(ncod_pl,muses_raw,'meth','MSGR'),
                     pcag_fun(ncod_pl,yield_raw,'meth','YC'),
                     pcag_fun(ncod_pl,deded_raw,'meth','MR'),
                     pcag_fun(ncod_pl,yeast_raw,'meth','Y-related'),
                     pcag_fun(ncod_pl,labbb_raw,'meth','LAB-related'),
                     pcag_fun(ncod_pl,aabbb_raw,'meth','AAB-related'))

alles_heat2 <- rbind(pcag_fun(ncod_pl,alles_raw,'country','ALL'),
                     pcag_fun(ncod_pl,muses_raw,'country','MSGR'),
                     pcag_fun(ncod_pl,yield_raw,'country','YC'),
                     pcag_fun(ncod_pl,deded_raw,'country','MR'),
                     pcag_fun(ncod_pl,yeast_raw,'country','Y-related'),
                     pcag_fun(ncod_pl,labbb_raw,'country','LAB-related'),
                     pcag_fun(ncod_pl,aabbb_raw,'country','AAB-related'))

alles_heat3 <- rbind(pcag_fun(ncod_pl,alles_raw,'temp','ALL'),
                     pcag_fun(ncod_pl,muses_raw,'temp','MSGR'),
                     pcag_fun(ncod_pl,yield_raw,'temp','YC'),
                     pcag_fun(ncod_pl,deded_raw,'temp','MR'),
                     pcag_fun(ncod_pl,yeast_raw,'temp','Y-related'),
                     pcag_fun(ncod_pl,labbb_raw,'temp','LAB-related'),
                     pcag_fun(ncod_pl,aabbb_raw,'temp','AAB-related'))

alles_heat4 <- rbind(pcag_fun(ncod_pl,alles_raw,'turn','ALL'),
                     pcag_fun(ncod_pl,muses_raw,'turn','MSGR'),
                     pcag_fun(ncod_pl,yield_raw,'turn','YC'),
                     pcag_fun(ncod_pl,deded_raw,'turn','MR'),
                     pcag_fun(ncod_pl,yeast_raw,'turn','Y-related'),
                     pcag_fun(ncod_pl,labbb_raw,'turn','LAB-related'),
                     pcag_fun(ncod_pl,aabbb_raw,'turn','AAB-related'))


fdf <- rbind(alles_heat,alles_heat1,alles_heat2,alles_heat3,alles_heat4)

fdf <- fdf %>%  mutate(feat = recode(feat, 
                                       `country` = '(a)',
                                       `cult` = '(b)',
                                       `turn` = '(d)',
                                       `meth` = '(c)',
                                       `temp` = '(e)'))

plt_fun <- function(df){
  p1<-ggplot(data=df,aes(y=group,x=mm,fill=V1))+
    geom_tile()+
    scale_fill_gradient(low='white',high='red',na.value=NA,name=expression(paste(tilde('D')['M'])))+
    facet_wrap(~feat,ncol=1)+
    theme_minimal()+ 
    theme(strip.text.x = element_text(hjust = -0.01),
          strip.text = element_text(size=18),
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
    scale_y_discrete(limits = rev(levels(df$group)))
  return(p1)
}


p1<-plt_fun(test)

save_plot(paste0(outpth,'figS27.pdf'),p1,base_height = 20.5,base_width =10.2)


########################################################################
# FIGURES 3 and 4 in Main Manuscript and Figures S28 and S29 in Supple-#
# mentary Material                                                     #
#                                                                      #
# explanation of functions:                                            #
#                                                                      #
# plt_func = Function that returns a list with panels (a) and (b) for  #
#            PCA figures. Arguments:                                   #
#            df = dataframe resulted from any dg_list_ function        #
#            feat = feature code name to be object of PCA              #
#                   'cult' = Cultivar                                  #
#                   'meth' = Fermentation method                       #
#                   'country' = Country of origin                      #
#                   'temp' = Controlled temperature                    #
#                   'turning' = Turning of fermenting mass             #
#            lab and lab2 = labels to left and right panel respective- #
#                           ly. (a) and (b)                            #
#            mdl = model numeric code                                  #
#                                                                      #
# splt = Function to plot output from plt_func. Arguments:             #
#        lpair = output from plt_func                                  #
#        ext_name = name of final figure                               #
#                                                                      #
########################################################################

plt_func <- function(df,feat,lab,lab2,mdl){
  df <- do.call('rbind',df)
  fl <- names(which(table(df[,feat])==7600))
  if (length(fl)!=0){
    for (j in 1:length(fl)){
      df <- subset(df,get(feat)!=fl[j])
    }
  }
  
  if (feat=='cult'){
    df <- subset(df,cult!='NA')
    cols <- c("cf"="#D43F3AFF",
              "nt"="#EEA236FF",
              "un1"="#5CB85CFF",
              "un2"="#357EBDFF")
  }
  
  if (feat=='turn'){
    df <- subset(df,turn!='NA')
    df$turn <- factor(df$turn,levels=c('tr','ntr'),labels=c('tr','ntr'))
    cols <- c("tr"="#D43F3AFF",
              "ntr"="#EEA236FF")
  }
  
  if (feat=='meth'){
    df <- subset(df,meth!='NA')
    cols <- c("hp"="#D43F3AFF",
              "pt"="#EEA236FF",
              "st"="#5CB85CFF",
              "wb"="#357EBDFF")
  }
  
  if (feat=='country'){
    df <- subset(df,country!='NA')
    cols <- c("BR"="#D43F3AFF",
              "EC"="#EEA236FF",
              "GH"="#5CB85CFF")
  }
  
  if (feat=='temp'){
    df <- subset(df,temp!='NA')
    cols <- c("ctrl"="#D43F3AFF",
              "nctrl"="#EEA236FF")
  }
  
  pca_df <- prcomp(df[c(11:(ncol(df)))],retx=T,center=T,scale=F)
  eig.val <- get_eigenvalue(pca_df)
  PC1.expl <- round(eig.val[1,2],2)
  PC2.expl <- round(eig.val[2,2],2)
  res.ind <- get_pca_ind(pca_df)
  PC1.ind <- res.ind$coord[,1]
  PC2.ind <- res.ind$coord[,2]
  res.var <- get_pca_var(pca_df)
  PC1.var <- res.var$coord[,1]
  PC2.var <- res.var$coord[,2]
  labs.var <- rownames(res.var$coord)
  
  dfp <- data.frame(PC1.ind,PC2.ind,df[,feat])
  colnames(dfp)[3] <- 'feat'
  
  set.seed(100785)
  dfp_sub <- dfp %>%
    sample_n(760)
  
  dfl <- data.frame(PC1.var,PC2.var,labs.var)
  
  if (mdl==2){
    lbs_df <- data.frame(cbind(c('mu1','mu2','mu3','mu4','mu5','mu6','k1','k2','k3','yc1','yc2','yc3','yc4','yc5',
                                 'yc6','yc7','yc8','yc9','yc10','yc11','yc12','yc13','yc14','yc15','yc16'),
                               c('$\\mu_{max}^{Y_{Glc}}$','$\\mu_{max}^{Y_{Fru}}$','$\\mu_{max}^{LAB_{Glc}}$',
                                 '$\\mu_{max}^{AAB_{EtOH}}$','$\\mu_{max}^{AAB_{LA}}$','$\\mu_{max}^{LAB_{Fru}}$',
                                 '$k_{Y}$','$k_{LAB}$','$k_{AAB}$','$Y_{Glc|Y}$','$Y_{Glc|LAB}$',
                                 '$Y_{Fru|Y}$','$Y^{Glc}_{EtOH|Y}$','$Y^{Fru}_{EtOH|Y}$','$Y_{EtOH|AAB}$',
                                 '$Y^{Glc}_{LA|LAB}$','$Y_{LA|AAB}$','$Y^{Glc}_{Ac|LAB}$','$Y^{EtOH}_{Ac|AAB}$',
                                 '$Y^{LA}_{Ac|AAB}$','$Y_{Fru|LAB}$','$Y^{Glc}_{EtOH|LAB}$','$Y^{Fru}_{EtOH|LAB}$',
                                 '$Y^{Fru}_{LA|LAB}$','$Y^{Fru}_{Ac|LAB}$')))  
  }
  
  if (mdl==23){
    lbs_df <- data.frame(cbind(c('mu1','mu2','mu3','mu4','mu5','mu6','k1','k2','k3','yc1','yc2','yc3','yc4','yc5',
                                 'yc6','yc7','yc8','yc9','yc10','yc11','yc12','yc13','yc14','yc15','yc16',
                                 'yc17','yc18'),
                               c('$\\mu_{max}^{Y_{Glc}}$','$\\mu_{max}^{Y_{Fru}}$','$\\mu_{max}^{LAB_{Glc}}$',
                                 '$\\mu_{max}^{AAB_{EtOH}}$','$\\mu_{max}^{AAB_{LA}}$','$\\mu_{max}^{LAB_{Fru}}$',
                                 '$k_{Y}$','$k_{LAB}$','$k_{AAB}$','$Y_{Glc|Y}$','$Y_{Glc|LAB}$',
                                 '$Y_{Fru|Y}$','$Y^{Glc}_{EtOH|Y}$','$Y^{Fru}_{EtOH|Y}$','$Y_{EtOH|AAB}$',
                                 '$Y^{Glc}_{LA|LAB}$','$Y_{LA|AAB}$','$Y^{Glc}_{Ac|LAB}$','$Y^{EtOH}_{Ac|AAB}$',
                                 '$Y^{LA}_{Ac|AAB}$','$Y_{Fru|LAB}$','$Y^{Glc}_{EtOH|LAB}$','$Y^{Fru}_{EtOH|LAB}$',
                                 '$Y^{Fru}_{LA|LAB}$','$Y^{Fru}_{Ac|LAB}$','$Y^{Glc}_{Ac|Y}$','$Y^{Fru}_{Ac|Y}$')))  
  }
  
  if (mdl==134){
    lbs_df <- data.frame(cbind(c('mu1','mu2','mu3','mu4','mu5','mu6','k1','k2','k3','yc1','yc2',
                                 'yc3','yc4','yc5','yc6','yc7','yc8','yc9','yc10','yc11',
                                 'yc12','yc13','ev1','ev2','ev3','yc14','yc15'),
                               c('$\\mu_{max}^{Y_{Glc}}$','$\\mu_{max}^{Y_{Fru}}$','$\\mu_{max}^{LAB_{Glc}}$',
                                 '$\\mu_{max}^{AAB_{EtOH}}$','$\\mu_{max}^{AAB_{LA}}$','$\\mu_{max}^{Y_{LA}}$',
                                 '$k_{Y}$','$k_{LAB}$','$k_{AAB}$','$Y_{Glc|Y}$','$Y_{Glc|LAB}$',
                                 '$Y_{Fru|Y}$','$Y^{Glc}_{EtOH|Y}$','$Y^{Fru}_{EtOH|Y}$','$Y_{EtOH|AAB}$',
                                 '$Y^{Glc}_{LA|LAB}$','$Y_{LA|AAB}$','$Y^{Glc}_{Ac|LAB}$','$Y^{EtOH}_{Ac|AAB}$',
                                 '$Y^{LA}_{Ac|AAB}$','$Y^{Glc}_{Ac|Y}$','$Y^{Fru}_{Ac|Y}$','$b_{EtOH}$',
                                 '$b_{LA}$','$b_{Ac}$','$Y_{LA|Y}$','$Y^{LA}_{EtOH|Y}$')))  
    
  }
  
  
  
  dfl <- left_join(dfl,lbs_df,by=c('labs.var'='X1'))
  
  plp1 <- ggplot(data=dfp,aes(x=PC1.ind,PC2.ind,colour=feat))+
    geom_point(data=dfp_sub,aes(shape=feat))+
    theme_bw() +
    theme(strip.text.x = element_text(hjust = -0.01),
          panel.grid.major = element_line(colour = "grey90",linetype='dashed'),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines'),
          strip.background = element_blank(),
          strip.text = element_text(size=16),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14),
          legend.position = "bottom"
    )+
    ggtitle(lab)+
    scale_colour_manual(values=cols)+
    scale_fill_manual(values=cols)+
    xlab(paste('Scores PC1 (',PC1.expl,'%',')', sep ='')) +
    ylab(paste('Scores PC2 (',PC2.expl,'%',')', sep ='')) +
    labs(fill='Groups',colour='Groups',shape='Groups')+
    stat_ellipse(aes(fill=feat),type='norm',level=0.90,geom='polygon',alpha=0.2,linetype='blank')
  
  plp2 <- ggplot(data=dfl,aes(x=PC1.var,y=PC2.var,label=TeX(X2,output='character')))+
    geom_point(size=2,colour="#D43F3AFF")+
    geom_text_repel(parse=T,size=4)+
    theme_bw() +
    theme(strip.text.x = element_text(hjust = -0.01),
          panel.grid.major = element_line(colour = "grey90",linetype='dashed'),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines'),
          strip.background = element_blank(),
          strip.text = element_text(size=16),
          axis.text=element_text(size=14),
          axis.title=element_text(size=14)
    )+
    ggtitle(lab2)+
    xlab(paste('Loadings PC1 (',PC1.expl,'%',')', sep ='')) +
    ylab(paste('Loadings PC2 (',PC2.expl,'%',')', sep ='')) 
  
  plt <- list(plp1,plp2)
  return(plt)
  
}

splt <- function(lpair,ext_name){
  plt <- plot_grid(lpair[[1]],lpair[[2]],ncol=2,align='hv',axis='bt')
  save_plot(paste0(outpth,ext_name),plt,base_height = 6.5,base_width = 11.7)
}

# creating dataframes to compute PCA
m23_alles <- dg_list_all('mm23') 
m2_alles <- dg_list_all('mm2')
m134_alles <- dg_list_all('mm134')
m134_lab <- dg_list_lab(m134_alles)

# creating panels for figures
pair1 <- plt_func(m23_alles,'cult','(a)','(b)',23) 
pair2 <- plt_func(m23_alles,'temp','(a)','(b)',23)
pair3 <- plt_func(m2_alles,'meth','(a)','(b)',2)
pair4 <- plt_func(m134_lab,'country','(a)','(b)',134)

# plotting figures
splt(pair1,'figure3.pdf') #Figure 3 main manuscript
splt(pair2,'figure4.pdf') #Figure 4 main manuscript
splt(pair3,'figureS28.pdf')   #Figure S28 supplementary material
splt(pair4,'figureS29.pdf') #Figure S29 supplementary material