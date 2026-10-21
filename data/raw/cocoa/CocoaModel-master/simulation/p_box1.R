# R Script for the reproduction of the results of "A mathematical model of cocoa bean fermentation"
# by Mauricio Moreno-Zambrano, Sergio Grimbs, Matthias S. Ullrich and Marc-Thorsten Huett
#
# Dataset: Papalexandratou, et al. (2011). Fermentation trial: Box 1
#
# 	  Mauricio Moreno-Zambrano, 2017
#     Jacobs University Bremen

# HEADER BLOCK

rm(list = ls()) #remove any object from R environment

setwd("path") #set working directory

require(rstan) #call to rstan library

rstan_options(auto_write = TRUE)

#DATA BLOCK

T = 13		#number of time points minus 1
x0 = c(55.482,  49.669,  0,  0,  1.379,  0.175424909, 0.003047264, 2.86522E-06)	#initial concentrations of Glc, Fru, EtOH, LA, Ac, Y, LAB and AAB respectively
t0 = 0	#initial time
ts = c(6, 12, 24, 30, 36, 48, 54, 60, 72, 84, 96, 120, 144)	#time points except initial time point
x = structure(c(26.802, 16.311, 13.706, 13.248, 7.54, 1.278, 2.694, 3.08, 0.757, 2.465, 2.673, 3.183, 7.631, #end of Glc time series
                31.959, 18.937, 19.238, 15.217, 10.351, 0.809, 1.007, 1.205, 0.757, 1.996, 2.673, 3.183, 2.287, #end of Fru time series
                0.201, 2.324, 2.46, 1.996, 3.173, 4.025, 6.033, 5.685, 3.605, 1.617, 0.576, 0.964, 1.56, #end of EtOH time series
                0.148, 0.342, 0.548, 1.033, 0.658, 4.873, 3.194, 2.694, 2.72, 1.928, 2.482, 1.883, 1.63, #end of LA time series
                0.963, 1.216, 1.467, 1.996, 1.211, 3.333, 6.033, 4.323, 8.685, 5.658, 11.935, 11.354, 11.511, #end of Ac time series
                0.119149235, 0.212369067, 0.068721283, 0.045508368, 0.052250597, 0.035488795, 0.020468747, 0.03750518, 0.004509114, 0.00048539, 0.000190586, 2.91133E-06, 1.45912E-06, #end of Y time series
                0.006575216, 0.02528774, 0.37403308, 0.364678377, 0.36383964, 0.073945204, 0.019182712, 0.027410062, 0.001068833, 0.00101137, 0.023873166, 0.00458047, 0.000150977, #end of LAB time series
                7.71184E-06, 7.09836E-06, 1.0427E-05, 7.90966E-06, 1.10194E-05, 6.05561E-05, 0.000127104, 0.000700097, 0.000130665, 0.000127104, 0.005347589, 0.008514478, 3.38188E-05), #end of AAB time series
              .Dim=c(13, 8)) #data stored in matrix x

#EXECUTION BLOCK

x1 = rbind(x0, x) #dummy matrix

scl = c(max(x1[,1]), max(x1[,2]), max(x1[,3]), max(x1[,4]), max(x1[,5]), max(x1[,6]), max(x1[,7]), max(x1[,8])) #scaling vector with maximum values for each time series

ini = function(){ #function to generate initial guessed values for the parameter estimation routine
  list(mu1=abs(rnorm(1,0.5,0.2)), mu2=abs(rnorm(1,0.5,0.2)), mu3=abs(rnorm(1,0.5,0.2)), mu4=abs(rnorm(1,0.5,0.2)), mu5=abs(rnorm(1,0.5,0.2)), 
       ks1=abs(rnorm(1,0.5,0.2)), ks2=abs(rnorm(1,0.5,0.2)), ks3=abs(rnorm(1,0.5,0.2)), ks4=abs(rnorm(1,0.5,0.2)), ks5=abs(rnorm(1,0.5,0.2)), 
       k1=abs(rnorm(1,0.5,0.2)), k2=abs(rnorm(1,0.5,0.2)), k3=abs(rnorm(1,0.5,0.2)), yc1=abs(rnorm(1,0.5,0.2)), yc2=abs(rnorm(1,0.5,0.2)), 
       yc3=abs(rnorm(1,0.5,0.2)), yc4=abs(rnorm(1,0.5,0.2)), yc5=abs(rnorm(1,0.5,0.2)), yc6=abs(rnorm(1,0.5,0.2)), yc7=abs(rnorm(1,0.5,0.2)), 
       yc8=abs(rnorm(1,0.5,0.2)), yc9=abs(rnorm(1,0.5,0.2)), yc10=abs(rnorm(1,0.5,0.2)), yc11=abs(rnorm(1,0.5,0.2)), sigma=abs(rnorm(1,0.5,0.2)))
}

fit = stan("model.stan", #file containing the Stan model
            data=c("x0",  "t0",  "ts",  "x",  "T", "scl"), 
            control=list(adapt_delta=0.985, 
                         stepsize=0.01, 
                         max_treedepth=12), 
            warmup = 1000, 
            init = ini, 
	          refresh = 100, 
            cores = min(4,  parallel::detectCores()), 
            chains=4,  iter=3000,  seed=100785)

save(fit,  file="output_pbox1.Rsave") #storage of the obtained output
