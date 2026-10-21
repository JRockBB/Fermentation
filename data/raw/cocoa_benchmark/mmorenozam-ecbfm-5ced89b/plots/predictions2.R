#######################################################################
#                                                                     #
# R Script for plotting posterior predictions in supplementary mate-  #
# rial of:                                                            #
#                                                                     #
# "Exploring cocoa bean fermentation processes by kinetic modelling"  #
#                                                                     #
# This script contains the functions for plotting Figures S32 to S34  #
# of Supplementary Material, Section 4.                               #
#                                                                     #
# Code made by Mauricio Moreno-Zambrano                               #
# 2021                                                                #
#                                                                     #
#######################################################################

#loading necessary R packages

library(cowplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

setwd("/home")      #set working directory accordingly
dpath<- '/home/'    #character string with the directory where Stan outputs are saved
outpth <- '/home/'  #character string with the directory where plots are going to be saved  

########################################################################
# FIGURES S32 TO S34 OF SUPPLEMENTARY MATERIAL                         #
#                                                                      #
# explanation of functions:                                            #
#                                                                      #
# swf = plotting function for Figures S32 to S34. Arguments:           #
#       dsn = dataset to extract posterior predictions:                #
#              1: brwb1                                                #
#              2: brwb2                                                #
#              3: ghhp1                                                #
#        fitnam = name of Rsave file where Stan output is stored.      #
#        mm = model iteration code (without commas).                   #
#        bh = base hight of the plot                                   #
#        bw = base width of the plot                                   #
#        lp = right plot margin                                        #
#                                                                      #
########################################################################

swf <- function(dsn,fitnam,mm,bh,bw,lp){
  if (dsn==1){
    T = 13
    x0 = c(55.482, 49.669, 0, 0, 1.379, 0.175424909,0.003047264,2.86522E-06,25.628)
    t0 = 0
    ts = c(6,12,24,30,36,48,54,60,72,84,96,120,144)
    
    x = structure(c(26.802,16.311,13.706,13.248,7.54,1.278,2.694,3.08,0.757,2.465,2.673,3.183,7.631,
                    31.959,18.937,19.238,15.217,10.351,0.809,1.007,1.205,0.757,1.996,2.673,3.183,2.287,
                    0.201,2.324,2.46,1.996,3.173,4.025,6.033,5.685,3.605,1.617,0.576,0.964,1.56,
                    0.148,0.342,0.548,1.033,0.658,4.873,3.194,2.694,2.72,1.928,2.482,1.883,1.63,
                    0.963,1.216,1.467,1.996,1.211,3.333,6.033,4.323,8.685,5.658,11.935,11.354,11.511,
                    0.119149235,0.212369067,0.068721283,0.045508368,0.052250597,0.035488795,0.020468747,0.03750518,0.004509114,0.00048539,0.000190586,2.91133E-06,1.45912E-06,
                    0.006575216,0.02528774,0.37403308,0.364678377,0.36383964,0.073945204,0.019182712,0.027410062,0.001068833,0.00101137,0.023873166,0.00458047,0.000150977,
                    7.71184E-06,7.09836E-06,1.0427E-05,7.90966E-06,1.10194E-05,6.05561E-05,0.000127104,0.000700097,0.000130665,0.000127104,0.005347589,0.008514478,3.38188E-05,
                    29.1,30.554,33.702,34.631,35.397,36.486,36.97,38.787,40.764,47.71,46.214,47.22,48.631),
                  .Dim=c(13,9))
    
    lbl<-'brwb1_'
  }
  if (dsn==2){
    T = 13
    x0 = c(42.936,67.249,0.318,0,0.942,0.004754351,0.000221264,3.86508E-06,25.234)
    t0 = 0
    ts = c(6,12,24,30,36,48,54,60,72,84,96,120,144)
    
    x = structure(c(23.452,8.264,4.85,1.27,1.62,2.32,0.294,1.467,1.892,5.975,7.863,6.521,7.19,
                    21.349,14.478,10.974,2.641,0.889,4.057,1.756,1.467,3.812,8.351,9.416,9.081,9.109,
                    0.386,2.211,4.496,4.912,6.714,4.307,7.057,4.097,1.044,1.642,0.9,0.341,0.29,
                    0.165,0.679,1.75,0.393,0.823,0.992,1.161,0.968,1.609,4.123,2.476,2.565,1.642,
                    1.079,1.194,2.532,1.329,2.183,1.718,2.272,2.733,6.568,11.767,9.083,18.625,17.187,
                    0.069998907,0.243271515,0.037852212,0.022340416,0.038823194,0.015491421,0.055474227,0.000444725,0.001542024,7.90845E-05,0.00016676,0.000144574,9.48618E-06,
                    0.002753658,0.026786133,0.016669018,0.051749959,0.125866459,0.54564479,0.699697002,0.05151219,0.11800761,0.013928682,0.024149604,0.002923547,0.001451811,
                    8.37834E-06,1.41957E-05,2.74889E-05,4.53062E-05,9.06062E-06,3.61541E-05,0.002723692,0.001033137,0.001436012,6.57897E-05,0.026987213,0.00468984,9.33594E-06,
                    27.311,29.428,33.744,35.169,36.553,38.507,40.502,42.782,44.776,48.726,46.158,44.159,47.698),
                  .Dim=c(13,9))
    lbl<-'brwb2_'
  }
  if (dsn==3){
    T = 16
    x0 = c(51.963,57.741,0,0,0,0.029180401,0.007868827,3.36634E-06,28.121)
    t0 = 0
    ts = c(6,12,18,24,30,36,42,48,54,60,66,72,84,96,120,144)
    x = structure(c(43.884,24.25,26.098,11.871,4.607,1.714,4.821,5.409,3.701,3.696,1.544,4.428,1.086,2.337,2.837,3.486,
                    49.588,35.584,33.876,19.871,8.607,2.825,4.821,7.409,4.812,7.548,1.544,5.539,1.086,2.337,2.837,3.486,
                    2.97,9.081,9.62,6.067,11.211,16.213,10.215,20.106,22.492,5.637,5.636,13.85,4.782,9.3,4.267,7.448,
                    0.267,0.801,1.256,1.745,5.944,3.246,7.787,7.991,6.943,8.593,8.296,6.85,8.021,7.667,7.209,8.117,
                    0.152,1.616,0.979,4.033,5.121,2.759,3.541,4.278,4.141,6.139,3.219,5.319,4.328,3.621,4.208,5.93,
                    0.026187332,0.955193281,0.238831309,0.045299276,0.01021154,0.06608323,0.008262115,0.000554742,0.000267973,0.006013001,0.002203389,0.000375916,0.000889388,0.000948618,0.000110177,0.000105461,
                    0.067128975,0.07959944,0.248834167,0.20179482,0.416783016,0.401707567,0.566121975,0.386286929,0.372314554,0.199484893,0.167459586,0.161402409,0.14484717,0.008748025,0.006575216,0.00181514,
                    9.82105E-06,1.70278E-05,2.16351E-05,5.49741E-05,0.00010572,0.000177484,0.000259512,0.000681017,0.000364887,0.00118075,0.001914953,0.000676329,0.001127608,0.000237224,6.32642E-05,2.64337E-05,
                    27.967,31.839,33.9,35.508,37.165,38.521,39.876,40.527,41.58,42.482,43.687,43.733,42.768,42.71,42.694,44.138),
                  .Dim=c(16,9))
    lbl<-'ghhp1_'
  }
  
  x1 = rbind(x0,x)
  scl = c(max(x1[,1]),max(x1[,2]),max(x1[,3]),max(x1[,4]),max(x1[,5]),max(x1[,6]),max(x1[,7]),max(x1[,8]),max(x1[,9]))
  
  load(paste0(dpath,fitnam))
  
  pred = as.data.frame(fit, pars = "x_hat") %>%
    gather(factor_key = TRUE) %>%
    group_by(key) %>%
    summarize(lb = quantile(value, probs = 0.025,na.rm=T),
              median = quantile(value, probs = 0.5,na.rm=T),
              ub = quantile(value, probs = 0.975,na.rm=T)) 
  
  sp = split(pred,rep(1:9,each=T))
  
  for (i in 1:9){
    sp[i] = lapply(sp[i],"*",scl[i])
    sp[i] = mapply(cbind, sp[i], "obs"=list(x[,i]), SIMPLIFY=F)
    sp[i] = mapply(cbind, sp[i], "time"=list(ts), SIMPLIFY=F)
    sp[i] = mapply(rbind,list(c(NA,x0[i],x0[i],x0[i],x0[i],0)),sp[i],SIMPLIFY=F)
    sp[i] = mapply(cbind, sp[i], "met"=metb[i],SIMPLIFY = F)
    sp[i] = mapply(cbind, sp[i], "data"=lbl,SIMPLIFY= F)
    sp[i] = mapply(cbind, sp[i], "model"=mm, SIMPLIFY = F)
  }
  
  dfout <- do.call(rbind.data.frame,sp)
  
  lpt <- list()
  
  for (i in 1:9){
    lpt[[i]] <- ggplot(data=subset(dfout,met==metb[i]),aes(x=time,y=median))+
      geom_line(color="#D43F3AFF",size=0.7)+
      geom_point(aes(x=time,y=obs))+
      geom_line(aes(x = time, y = ub),color='#EEA2367F',linetype='dashed',size=1) +
      geom_line(aes(x = time, y = lb),color='#EEA2367F',linetype='dashed',size=1)+
      geom_ribbon(aes(ymin = lb, ymax = ub), fill = "#EEA2367F",alpha=0.3)+
      ggtitle(metb[i])+
      theme_bw()+
      theme(plot.title=element_text(hjust=-0.01,size=16),
            panel.grid.major = element_line(colour = "grey90",linetype='dashed'),
            panel.grid.minor = element_blank(),
            plot.margin = unit(c(0,lp,0,0),'lines'),
            axis.text=element_text(size=14),
            axis.title.y=element_text(size=14,margin=margin(0,-10,0,0)),
            axis.title.x=element_text(size=14,margin=margin(0,0,0,0))
      )+
      xlab(NULL)+
      ylab(NULL)
  }
  lpt[[1]] <- lpt[[1]] + ylab(expression(paste('mg g(pulp) '^{-1})))
  lpt[[4]] <- lpt[[4]] + ylab(expression(paste('mg g(pulp) '^{-1})))
  lpt[[7]] <- lpt[[7]] + ylab(expression(paste('mg g(pulp) '^{-1}))) + xlab('time (h)')
  lpt[[8]] <- lpt[[8]] + xlab('time (h)')
  lpt[[9]] <- lpt[[9]] + ylab('ºC') + xlab('time (h)')
  
  plt_f <- plot_grid(lpt[[1]],lpt[[2]],lpt[[3]],
                     lpt[[4]],lpt[[5]],lpt[[6]],
                     lpt[[7]],lpt[[8]],lpt[[9]],
                     ncol=3,align='hv',axis='bt')
  
  return(save_plot(paste0(outpth,lbl,mm,'.pdf'),plt_f,base_height = bh,base_width = bw))
}


#ploting Figures S32 to S34

swf(1,'brwb1_mt1.Rsave','1',12.45,14.9,1) #Figure S32 
swf(2,'brwb2_mt1.Rsave','1',12.45,14.9,1) #Figure S33
swf(3,'ghhp1_mt1.Rsave','1',12.45,14.9,1) #Figure S34
