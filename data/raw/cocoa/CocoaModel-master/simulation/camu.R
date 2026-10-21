# R Script for the reproduction of the results of "A mathematical model of cocoa bean fermentation"
# by Mauricio Moreno-Zambrano, Sergio Grimbs, Matthias S. Ullrich and Marc-Thorsten Huett
#
# Dataset: Camu, et al. (2007)
#
# 	  Mauricio Moreno-Zambrano, 2017
#     Jacobs University Bremen

# HEADER BLOCK

rm(list = ls()) #remove any object from R environment

setwd("path") #set working directory

require(rstan) #call to rstan library

rstan_options(auto_write = TRUE)

#DATA BLOCK

T = 16		#number of time points minus 1
x0 = c(51.963, 57.741, 0, 0, 0, 0.029180401, 0.007868827, 3.36634E-06) 	#initial concentrations of Glc, Fru, EtOH, LA, Ac, Y, LAB and AAB respectively
t0 = 0 	#initial time
ts = c(6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 84, 96, 120, 144) #time points except initial time point
x = structure(c(43.884, 24.25, 26.098, 11.871, 4.607, 1.714, 4.821, 5.409, 3.701, 3.696, 1.544, 4.428, 1.086, 2.337, 2.837, 3.486, #end of Glc time series
                49.588, 35.584, 33.876, 19.871, 8.607, 2.825, 4.821, 7.409, 4.812, 7.548, 1.544, 5.539, 1.086, 2.337, 2.837, 3.486, #end of Fru time series
                2.97, 9.081, 9.62, 6.067, 11.211, 16.213, 10.215, 20.106, 22.492, 5.637, 5.636, 13.85, 4.782, 9.3, 4.267, 7.448, #end of EtOH time series
                0.267, 0.801, 1.256, 1.745, 5.944, 3.246, 7.787, 7.991, 6.943, 8.593, 8.296, 6.85, 8.021, 7.667, 7.209, 8.117, #end of LA time series
                0.152, 1.616, 0.979, 4.033, 5.121, 2.759, 3.541, 4.278, 4.141, 6.139, 3.219, 5.319, 4.328, 3.621, 4.208, 5.93, #end of Ac time series
                0.026187332, 0.955193281, 0.238831309, 0.045299276, 0.01021154, 0.06608323, 0.008262115, 0.000554742, 0.000267973, 0.006013001, 0.002203389, 0.000375916, 0.000889388, 0.000948618, 0.000110177, 0.000105461, #end of Y time series
                0.067128975, 0.07959944, 0.248834167, 0.20179482, 0.416783016, 0.401707567, 0.566121975, 0.386286929, 0.372314554, 0.199484893, 0.167459586, 0.161402409, 0.14484717, 0.008748025, 0.006575216, 0.00181514, #end of LAB time series
                9.82105E-06, 1.70278E-05, 2.16351E-05, 5.49741E-05, 0.00010572, 0.000177484, 0.000259512, 0.000681017, 0.000364887, 0.00118075, 0.001914953, 0.000676329, 0.001127608, 0.000237224, 6.32642E-05, 2.64337E-05), #end of AAB time series 
              .Dim=c(16, 8)) #data stored in matrix x 

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

save(fit,  file="output_camu.Rsave") #storage of the obtained output
