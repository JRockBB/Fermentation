# R Script for the reproduction of the results of "A mathematical model of cocoa bean fermentation"
# by Mauricio Moreno-Zambrano, Sergio Grimbs, Matthias S. Ullrich and Marc-Thorsten Huett
#
# Dataset: Papalexandratou, et al. (2011). Fermentation trial: Box 2
#
#     Mauricio Moreno-Zambrano, 2017
#     Jacobs University Bremen

# HEADER BLOCK

rm(list = ls()) #remove any object from R environment

setwd("path") #set working directory

require(rstan) #call to rstan library

rstan_options(auto_write = TRUE)

#DATA BLOCK

T = 13		#number of time points minus 1
x0 = c(42.936, 67.249, 0.318, 0, 0.942, 0.004754351, 0.000221264, 3.86508E-06)	#initial concentrations of Glc, Fru, EtOH, LA, Ac, Y, LAB and AAB respectively
t0 = 0	#initial time
ts = c(6, 12, 24, 30, 36, 48, 54, 60, 72, 84, 96, 120, 144) #time points except initial time point
x = structure(c(23.452, 8.264, 4.85, 1.27, 1.62, 2.32, 0.294, 1.467, 1.892, 5.975, 7.863, 6.521, 7.19, #end of Glc time series
                21.349, 14.478, 10.974, 2.641, 0.889, 4.057, 1.756, 1.467, 3.812, 8.351, 9.416, 9.081, 9.109, #end of Fru time series
                0.386, 2.211, 4.496, 4.912, 6.714, 4.307, 7.057, 4.097, 1.044, 1.642, 0.9, 0.341, 0.29, #end of EtOH time series
                0.165, 0.679, 1.75, 0.393, 0.823, 0.992, 1.161, 0.968, 1.609, 4.123, 2.476, 2.565, 1.642, #end of LA time series
                1.079, 1.194, 2.532, 1.329, 2.183, 1.718, 2.272, 2.733, 6.568, 11.767, 9.083, 18.625, 17.187, #end of Ac time series
                0.069998907, 0.243271515, 0.037852212, 0.022340416, 0.038823194, 0.015491421, 0.055474227, 0.000444725, 0.001542024, 7.90845E-05, 0.00016676, 0.000144574, 9.48618E-06, #end of Y time series
                0.002753658, 0.026786133, 0.016669018, 0.051749959, 0.125866459, 0.54564479, 0.699697002, 0.05151219, 0.11800761, 0.013928682, 0.024149604, 0.002923547, 0.001451811, #end of LAB time series
                8.37834E-06, 1.41957E-05, 2.74889E-05, 4.53062E-05, 9.06062E-06, 3.61541E-05, 0.002723692, 0.001033137, 0.001436012, 6.57897E-05, 0.026987213, 0.00468984, 9.33594E-06), #end of AAB time series
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

save(fit,  file="output_pbox2.Rsave") #storage of the obtained output
