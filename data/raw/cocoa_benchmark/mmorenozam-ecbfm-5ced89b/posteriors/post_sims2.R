#######################################################################
#                                                                     #
# R Script for computation of posterior predictions of the manuscript #
# "Exploring cocoa bean fermentation processes by kinetic modelling"  #
#                                                                     #
# The following script contains data from 4 fermentation trials co-   #
# rresponding to Camu et al. 2008 where time series for lactic acid   #
# were generated. For a more detailed description of datasets, see    #
# Supplementary Material, Section 1, and Supplementary Table S1.      #
#                                                                     #
# This code works in prompt shell of Linux OS. For its execution it   #
# needs seven arguments:                                              #
# d = dataset (integer from 1 to 4):                                  #
#     1:  ghhp2                                                       #
#     2:  ghhp3                                                       #
#     3:  ghhp4                                                       #
#     4:  ghhp5                                                       #
#                                                                     #
# mi = model iteration (integer), refers to model iteration code, for #
#      instance, MI(1,2,3), here is called as 123                     #
#                                                                     #
# The following are Stan parameters. For further understanding, see   #
# reference manual for Stan. Between parenthesis are noted the values #
# used in our paper                                                   #
# ad = adapt_delta parameter (0.995 to 0.999)                         #
#                                                                     #
# st = step_size parameter (0.1)                                      #
#                                                                     #
# mt = max_treedepth parameter (10)                                   #
#                                                                     #
# s  = seed parameter (141085)                                        #
#                                                                     #
# ni = number of iterations (3000)                                    #
#                                                                     #
# Hence, for running this script you must follow the following form:  #
#                                                                     #
# $ Rscript post_sims2.R d mi ad st mt s ni                           #
#                                                                     #
# Once the script is done running, an .Rsave file is saved with the   #
# dataset coding as name containing the Stan output.                  #
#                                                                     #
# Code written by Mauricio Moreno-Zambrano                            #
# 2021                                                                #
#                                                                     #
#######################################################################

setwd('/home') # directory name, change accordingly if needed
rm(list = ls())
args = commandArgs(TRUE)

library(rstan)

rstan_options(auto_write = TRUE)

# datasets ordered as mentioned above

if (args[1]==1){
  T = 15
  xi_r = c(47.956,51.386,0.286,0,0.048539049,0.004044921,0.000906062)
  t0 = 0
  ts = c(6,12,18,24,30,36,42,48,54,60,72,84,96,120,144)
  
  x = structure(c(34.158,24.028,12.064,7.996,7.438,4.885,1.854,1.615,2.094,2.094,3.131,2.333,3.21,5.922,0.259,
                  41.575,23.789,15.972,7.916,12.303,6.64,3.529,3.45,3.45,3.529,5.444,4.088,4.806,7.278,2.014,
                  1.883,5.389,3.098,10.491,4.833,4.799,1.501,3.896,2.751,1.883,1.467,0.564,1.189,1.501,0.148,
                  0,1.536,1.605,1.64,3.827,4.07,4.833,11.359,9.623,6.465,15.802,11.289,9.207,11.845,5.632,
                  0.278670668,0.394540199,0.508266234,0.189710452,0.013934496,0.021090713,0.013934496,0.000478731,0.000244958,0.000000015,0.000000015,0.000000015,0.000000015,0.000000015,0.000000015,
                  0.032954142,0.03287835,0.024712121,0.015809204,0.141550045,0.333357333,0.402633599,0.061646725,0.322040145,0.003740331,0.014220341,0.04177438,0.00020602,0.00018068,0.000126739,
                  0.001880001,0.00736475,0.001816176,0.002495503,0.00136194,0.000432671,0.000260111,0.000228643,0.000925035,0.000313443,0.001349454,0.000628287,0.000364887,0.000483235,0.000450981 #aab
  ),
  .Dim=c(15,7))
  lbl <- 'ghhp2_'
} 

if (args[1]==2){
  T = 15
  xi_r = c(40.246,50.851,0,0,0.059304993,0.014618742,0.001107027)
  t0 = 0
  ts = c(6,12,18,24,30,36,42,48,54,60,72,84,96,120,144)
  
  x = structure(c(29.88,35.62,31.234,18.236,8.827,8.108,1.968,4.597,2.126,13.288,6.509,11.291,4.034,1.32,8.971,
                  47.9,47.66,43.353,27.964,15.604,12.813,5.636,6.99,3.162,8.025,8.263,9.537,5.867,4.031,8.413,
                  0.254,7.115,13.419,11.193,6.668,5.068,5.835,11.2,8.625,6.677,5.67,6.892,7.105,9.863,9.242,
                  0.185,0.708,1.233,1.826,0.573,0.158,2.284,4.305,3.889,1.976,1.213,1.774,1.708,1.68,3.289,
                  0.32143359,0.0539624,0.358997363,0.229663119,0.005117894,0.001599894,0.004721622,0.007483267,0.022860791,0.002198322,0.000790845,0.000508266,0.000721259,0.000363992,0.00016147,
                  0.039166072,0.060243475,0.044556392,0.028967433,0.074115666,0.212280457,0.178200949,0.142531223,0.795994401,0.23708824,0.218227769,0.053077445,0.057532072,0.049763396,0.02937041,
                  0.002128914,0.010142804,0.006701284,0.00416062,0.000558673,0.000446846,0.001500231,0.013936638,0.001449299,0.000393693,0.000277433,0.000419912,0.001023665,0.00079279,0.000767641 #aab
  ),
  .Dim=c(15,7))
  lbl <- 'ghhp3_'
}

if (args[1]==3){
  T = 15
  xi_r = c(41.573,44.215,0.617,0,0.02761158,0.001314952,8.4949E-05)
  t0 = 0
  ts = c(6,12,18,24,30,36,42,48,54,60,72,84,96,120,144)
  
  x = structure(c(34.776,24.698,7.016,13.984,13.188,3.67,7.276,12.482,1.604,3.289,3.459,3.948,5.479,3.496,4.077,
                  37.737,30.14,12.698,18.465,17.35,6.391,13.758,17.204,4.005,6.01,6.66,6.988,7.959,7.979,7.678,
                  1.871,9.427,13.117,14.546,12.351,18.271,9.358,5.632,8.417,4.309,2.812,2.046,2.846,0.514,0.339,
                  0.339,0.653,2.429,1.001,0.827,3.16,3.96,6.119,8.522,10.367,11.064,10.89,11.029,12.073,9.218,
                  0.055602108,0.41695699,0.641344329,0.338915366,0.124191325,0.014193557,0.155629262,0.15103975,0.050244816,0.00052371,0.001130033,0.00104494,0.000354072,0.000133996,0.000000015,
                  0.071106616,0.247121205,0.831591445,0.251715531,0.734361691,0.64700854,0.560931737,0.494208275,0.485187957,0.184039063,0.15202325,0.147540079,0.005443898,0.002830805,0.003862869,
                  0.006624575,0.009123428,0.001202702,0.008168796,0.002637291,0.001858481,0.001227886,0.003200059,0.004407152,0.011485715,0.009821052,0.002265469,0.000944404,0.000781912,0.001462709 #aab
  ),
  .Dim=c(15,7))
  lbl <- 'ghhp4_'
}

if (args[1]==4){
  T = 15
  xi_r = c(49.309,51.391,0,0,0.025240111,0.007061712,0.000861307)
  t0 = 0
  ts = c(6,12,18,24,30,36,42,48,54,60,72,84,96,120,144)
  
  x = structure(c(45.631,15.927,13.69,7.45,5.934,10.343,3.783,5.55,3.314,4.761,6.694,2.861,5.596,11.384,12.127,
                  48.514,21.853,18.335,10.972,10.418,15.709,7.787,11.717,5.876,8.765,12.94,6.465,9.679,17.55,17.572,
                  2.876,2.527,13.479,17.859,20.257,13.024,11.215,10.762,13.682,5.683,11.696,10.164,8.597,6.438,2.713,
                  0,1.276,0.649,0.577,2.281,1.863,6.938,6.728,3.737,4.64,7.141,6.617,4.773,11.758,10.05,
                  0.054211479,0.227557555,0.333496484,0.074660563,0.171037468,0.026127103,0.025298295,0.019014778,0.023019255,0.003090945,0.000952996,0.000716294,0.000366515,0.000213841,0.000000015,
                  0.029642171,0.444539148,0.503396293,1.693986765,0.632280828,0.632280828,0.652995236,0.741156656,0.895179263,0.76543799,0.039619593,0.128798265,0.061788836,0.058198262,0.002055465,
                  0.003192699,0.004679054,0.005850029,0.000760603,0.005484765,0.003000254,0.003865076,0.003000254,0.00290508,0.000628287,0.003304898,0.004131978,0.009775929,0.00164497,0.000188 #aab
  ),
  .Dim=c(15,7))
  lbl <- 'ghhp5_'
}

x_d <- rbind(xi_r,x)

scl <- c()
x0 <- c()
for (i in 1:7){
  scl[i] <- max(x_d[,i])
  x0[i] <- xi_r[i]/scl[i]
  x[,i]<-x[,i]/scl[i]
}


x0 = c(x0,0)

#depending on chosen MI, initial values are randomly generated for MCMC-NUTS sampler

if (args[2]==0){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==2){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==3){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==4){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==5){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==12){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==13){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==14){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==15){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==23){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==24){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==25){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==34){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==35){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==45){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==123){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),y17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==124){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==125){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==134){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==135){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==145){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc13=abs(rnorm(1,0.5,0.2)),yc14=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==234){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==235){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==245){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1234){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1235){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}
if (args[2]==1245){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==1345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),ks1=abs(rnorm(1,0.5,0.2)),
         ks2=abs(rnorm(1,0.5,0.2)),ks3=abs(rnorm(1,0.5,0.2)),
         ks4=abs(rnorm(1,0.5,0.2)),ks5=abs(rnorm(1,0.5,0.2)),
         ks6=abs(rnorm(1,0.5,0.2)),ks7=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),ev1=abs(rnorm(1,0.5,0.2)),
         ev2=abs(rnorm(1,0.5,0.2)),ev3=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==2345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),yc21=abs(rnorm(1,0.5,0.2)),
         sigma=abs(rnorm(1,0.5,0.2)))
  }
}

if (args[2]==12345){
  ini = function(){
    list(mu1=abs(rnorm(1,0.5,0.2)),mu2=abs(rnorm(1,0.5,0.2)),
         mu3=abs(rnorm(1,0.5,0.2)),mu4=abs(rnorm(1,0.5,0.2)),
         mu5=abs(rnorm(1,0.5,0.2)),mu6=abs(rnorm(1,0.5,0.2)),
         mu7=abs(rnorm(1,0.5,0.2)),mu8=abs(rnorm(1,0.5,0.2)),
         ks1=abs(rnorm(1,0.5,0.2)),ks2=abs(rnorm(1,0.5,0.2)),
         ks3=abs(rnorm(1,0.5,0.2)),ks4=abs(rnorm(1,0.5,0.2)),
         ks5=abs(rnorm(1,0.5,0.2)),ks6=abs(rnorm(1,0.5,0.2)),
         ks7=abs(rnorm(1,0.5,0.2)),ks8=abs(rnorm(1,0.5,0.2)),
         k1=abs(rnorm(1,0.5,0.2)),k2=abs(rnorm(1,0.5,0.2)),
         k3=abs(rnorm(1,0.5,0.2)),yc1=abs(rnorm(1,0.5,0.2)),
         yc2=abs(rnorm(1,0.5,0.2)),yc3=abs(rnorm(1,0.5,0.2)),
         yc4=abs(rnorm(1,0.5,0.2)),yc5=abs(rnorm(1,0.5,0.2)),
         yc6=abs(rnorm(1,0.5,0.2)),yc7=abs(rnorm(1,0.5,0.2)),
         yc8=abs(rnorm(1,0.5,0.2)),yc9=abs(rnorm(1,0.5,0.2)),
         yc10=abs(rnorm(1,0.5,0.2)),yc11=abs(rnorm(1,0.5,0.2)),
         yc12=abs(rnorm(1,0.5,0.2)),yc13=abs(rnorm(1,0.5,0.2)),
         yc14=abs(rnorm(1,0.5,0.2)),yc15=abs(rnorm(1,0.5,0.2)),
         yc16=abs(rnorm(1,0.5,0.2)),yc17=abs(rnorm(1,0.5,0.2)),
         yc18=abs(rnorm(1,0.5,0.2)),yc19=abs(rnorm(1,0.5,0.2)),
         yc20=abs(rnorm(1,0.5,0.2)),yc21=abs(rnorm(1,0.5,0.2)),
         ev1=abs(rnorm(1,0.5,0.2)),ev2=abs(rnorm(1,0.5,0.2)),
         ev3=abs(rnorm(1,0.5,0.2)),sigma=abs(rnorm(1,0.5,0.2))
    )
  }
}

#Stan call:
fit = stan(paste0("mmc",args[2],".stan"),
           data=c("x0", "t0", "ts", "x", 'T'),
           control=list(adapt_delta=as.numeric(args[3]),
                        stepsize=as.numeric(args[4]),
                        max_treedepth=as.numeric(args[5])),
           warmup = round(as.numeric(args[7])*1/3,0),
           init = ini,
           refresh = 5,
           cores = min(4, parallel::detectCores()),
           chains=4, iter=as.numeric(args[7]), seed=as.numeric(args[6]))

save(fit, file=paste0(lbl,"mmc",args[2],".Rsave")) #saving output