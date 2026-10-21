# R Script for the reproduction of the results of "A mathematical model of cocoa bean fermentation"
# by Mauricio Moreno-Zambrano, Sergio Grimbs, Matthias S. Ullrich and Marc-Thorsten Huett
#
# This script extracts the estimated values of the parameters from the Stan outputs. It generates LaTeX code for a table containing
# the estimates for the three datasets, as well as the absolute value of the standardized mean difference, similar to Table 3 from
# our original article.
# 
# The function 'getSimsTable' was written by Charles Margossian and Bill Gillespie. Available from this source:
# https://github.com/stan-dev/stancon_talks/blob/master/2017/Contributed-Talks/05_margossian/tools/stanTools.R
#
# INSTRUCTIONS:
# 1. Set the working path to the folder where the Stan outputs were generated in your machine.
# 2. Run the script
#
#     Mauricio Moreno-Zambrano, 2017
#     Jacobs University Bremen

setwd('path') # change the working path according to your own machine

require(xtable) # R package for generating a LaTeX code with the generated estimates

getSimsTable <- function(x, ...){#function to extract parameters from Stan output. Thanks to Charles Margossian and Bill Gillespie
  require(dplyr)
  nChains <- dim(x)[2]
  nPost <- dim(x)[1]
  x %>%
    as.data.frame(...) %>%
    mutate(chain = rep(1:nChains, ea = nPost),
           iteration = rep(1:nPost, nChains))
}

load('output_camu.Rsave')
cmu = getSimsTable(fit)

# conversion factors for re-scaling parameters to real units
glu = 51.963
fru = 57.741
eto = 22.492
lac = 8.593
acc = 6.139
yst =  10^7.804*0.000000015
lab =  10^8.656*0.00000000125
aab =   10^6.835*0.00000000028
cmu$yc1 = cmu$yc1*(glu/yst)
cmu$yc2 = cmu$yc2*(glu/lab)
cmu$yc3 = cmu$yc3*(fru/yst)
cmu$yc4 = cmu$yc4*(eto/yst)
cmu$yc5 = cmu$yc5*(eto/yst)
cmu$yc6 = cmu$yc6*(eto/aab)
cmu$yc7 = cmu$yc7*(lac/lab)
cmu$yc8 = cmu$yc8*(lac/aab)
cmu$yc9 = cmu$yc9*(acc/lab)
cmu$yc10 = cmu$yc10*(acc/aab)
cmu$yc11 = cmu$yc11*(acc/aab)
cmu$ks1 = cmu$ks1*glu
cmu$ks2 = cmu$ks2*fru
cmu$ks3 = cmu$ks3*glu
cmu$ks4 = cmu$ks4*eto
cmu$ks5 = cmu$ks5*(lac/aab)
cmu$k1 = cmu$k1/eto
cmu$k2 = cmu$k2/lac
cmu$k3 = cmu$k3/(acc^2)

load('output_pbox1.Rsave')
pb1 = getSimsTable(fit)

glu = 55.482000000
fru = 49.669000000
eto = 6.033000000
lac = 4.873000000
acc = 11.935000000
yst =  10^7.151*0.000000015  
lab = 10^8.476 *0.00000000125
aab =  10^7.483*0.00000000028
pb1$yc1 = pb1$yc1*(glu/yst)
pb1$yc2 = pb1$yc2*(glu/lab)
pb1$yc3 = pb1$yc3*(fru/yst)
pb1$yc4 = pb1$yc4*(eto/yst)
pb1$yc5 = pb1$yc5*(eto/yst)
pb1$yc6 = pb1$yc6*(eto/aab)
pb1$yc7 = pb1$yc7*(lac/lab)
pb1$yc8 = pb1$yc8*(lac/aab)
pb1$yc9 = pb1$yc9*(acc/lab)
pb1$yc10 = pb1$yc10*(acc/aab)
pb1$yc11 = pb1$yc11*(acc/aab)
pb1$ks1 = pb1$ks1*glu
pb1$ks2 = pb1$ks2*fru
pb1$ks3 = pb1$ks3*glu
pb1$ks4 = pb1$ks4*eto
pb1$ks5 = pb1$ks5*(lac/aab)
pb1$k1 = pb1$k1/eto
pb1$k2 = pb1$k2/lac
pb1$k3 = pb1$k3/(acc^2)

load('output_pbox2.Rsave')
pb2 = getSimsTable(fit)

glu = 42.93600000
fru = 67.24900000
eto = 7.05700000 
lac = 4.12300000
acc = 18.62500000
yst = 10^7.210*0.000000015 
lab =  10^8.748*0.00000000125 
aab =  10^7.984*0.00000000028
pb2$yc1 = pb2$yc1*(glu/yst)
pb2$yc2 = pb2$yc2*(glu/lab)
pb2$yc3 = pb2$yc3*(fru/yst)
pb2$yc4 = pb2$yc4*(eto/yst)
pb2$yc5 = pb2$yc5*(eto/yst)
pb2$yc6 = pb2$yc6*(eto/aab)
pb2$yc7 = pb2$yc7*(lac/lab)
pb2$yc8 = pb2$yc8*(lac/aab)
pb2$yc9 = pb2$yc9*(acc/lab)
pb2$yc10 = pb2$yc10*(acc/aab)
pb2$yc11 = pb2$yc11*(acc/aab)
pb2$ks1 = pb2$ks1*glu
pb2$ks2 = pb2$ks2*fru
pb2$ks3 = pb2$ks3*glu
pb2$ks4 = pb2$ks4*eto
pb2$ks5 = pb2$ks5*(lac/aab)
pb2$k1 = pb2$k1/eto
pb2$k2 = pb2$k2/lac
pb2$k3 = pb2$k3/(acc^2)

# parameters' symbols in LaTeX
mu1 = '$\\mu_\\text{max}^{\\text{Y}_\\text{Glc}}$'
mu2 = '$\\mu_\\text{max}^{\\text{Y}_\\text{Fru}}$'
mu3 = '$\\mu_\\text{max}^\\text{LAB}$'
mu4 = '$\\mu_\\text{max}^{\\text{AAB}_\\text{EtOH}}$'
mu5 = '$\\mu_\\text{max}^{\\text{AAB}_\\text{LA}}$'
ks1 = '$K_{\\text{Glc}}^{\\text{Y}}$'
ks2 = '$K_{\\text{Fru}}^{\\text{Y}}$'
ks3 = '$K_{\\text{Glc}}^{\\text{LAB}}$'
ks4 = '$K_{\\text{EtOH}}^{\\text{AAB}}$'
ks5 = '$K_{\\text{LA}}^{\\text{AAB}}$'
k1 = '$k_{\\text{Y}}$'
k2 = '$k_{\\text{LAB}}$'
k3 = '$k_{\\text{AAB}}$'
y1 = '$Y_{\\text{Glc}|\\text{Y}}$'
y2 = '$Y_{\\text{Glc}|\\text{LAB}}$'
y3 = '$Y_{\\text{Fru}|\\text{Y}}$'
y4 = '$Y^{\\text{Glc}}_{\\text{EtOH}|\\text{Y}}$'
y5 = '$Y^{\\text{Fru}}_{\\text{EtOH}|\\text{Y}}$'
y6 = '$Y_{\\text{EtOH}|\\text{AAB}}$'
y7 = '$Y_{\\text{LA}|\\text{LAB}}$'
y8 = '$Y_{\\text{LA}|\\text{AAB}}$'
y9 = '$Y_{\\text{Ac}|\\text{LAB}}$'
y10 = '$Y^{\\text{EtOH}}_{\\text{Ac}|\\text{AAB}}$'
y11 = '$Y^{\\text{LA}}_{\\text{Ac}|\\text{AAB}}$'

par = c(mu1,mu2,mu3,mu4,mu5,ks1,ks2,ks3,ks4,ks5,k1,k2,k3,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11)

# extracting means and standard deviations from the datasets containing re-scaled parameters
mcmu = c(round(mean(cmu$mu1),6),round(mean(cmu$mu2),6),round(mean(cmu$mu3),6),round(mean(cmu$mu4),6),round(mean(cmu$mu5),6),
         round(mean(cmu$ks1),6),round(mean(cmu$ks2),6),round(mean(cmu$ks3),6),round(mean(cmu$ks4),6),round(mean(cmu$ks5),6),
         round(mean(cmu$k1),6),round(mean(cmu$k2),6),round(mean(cmu$k3),6),
         round(mean(cmu$yc1),6),round(mean(cmu$yc2),6),round(mean(cmu$yc3),6),round(mean(cmu$yc4),6),round(mean(cmu$yc5),6),
         round(mean(cmu$yc6),6),round(mean(cmu$yc7),6),round(mean(cmu$yc8),6),round(mean(cmu$yc9),6),round(mean(cmu$yc10),6),
         round(mean(cmu$yc11),6))

sdcmu = c(round(sd(cmu$mu1),6),round(sd(cmu$mu2),6),round(sd(cmu$mu3),6),round(sd(cmu$mu4),6),round(sd(cmu$mu5),6),
          round(sd(cmu$ks1),6),round(sd(cmu$ks2),6),round(sd(cmu$ks3),6),round(sd(cmu$ks4),6),round(sd(cmu$ks5),6),
          round(sd(cmu$k1),6),round(sd(cmu$k2),6),round(sd(cmu$k3),6),
          round(sd(cmu$yc1),6),round(sd(cmu$yc2),6),round(sd(cmu$yc3),6),round(sd(cmu$yc4),6),round(sd(cmu$yc5),6),
          round(sd(cmu$yc6),6),round(sd(cmu$yc7),6),round(sd(cmu$yc8),6),round(sd(cmu$yc9),6),round(sd(cmu$yc10),6),
          round(sd(cmu$yc11),6))

mpb1 = c(round(mean(pb1$mu1),6),round(mean(pb1$mu2),6),round(mean(pb1$mu3),6),round(mean(pb1$mu4),6),round(mean(pb1$mu5),6),
         round(mean(pb1$ks1),6),round(mean(pb1$ks2),6),round(mean(pb1$ks3),6),round(mean(pb1$ks4),6),round(mean(pb1$ks5),6),
         round(mean(pb1$k1),6),round(mean(pb1$k2),6),round(mean(pb1$k3),6),
         round(mean(pb1$yc1),6),round(mean(pb1$yc2),6),round(mean(pb1$yc3),6),round(mean(pb1$yc4),6),round(mean(pb1$yc5),6),
         round(mean(pb1$yc6),6),round(mean(pb1$yc7),6),round(mean(pb1$yc8),6),round(mean(pb1$yc9),6),round(mean(pb1$yc10),6),
         round(mean(pb1$yc11),6))

sdpb1 = c(round(sd(pb1$mu1),6),round(sd(pb1$mu2),6),round(sd(pb1$mu3),6),round(sd(pb1$mu4),6),round(sd(pb1$mu5),6),
          round(sd(pb1$ks1),6),round(sd(pb1$ks2),6),round(sd(pb1$ks3),6),round(sd(pb1$ks4),6),round(sd(pb1$ks5),6),
          round(sd(pb1$k1),6),round(sd(pb1$k2),6),round(sd(pb1$k3),6),
          round(sd(pb1$yc1),6),round(sd(pb1$yc2),6),round(sd(pb1$yc3),6),round(sd(pb1$yc4),6),round(sd(pb1$yc5),6),
          round(sd(pb1$yc6),6),round(sd(pb1$yc7),6),round(sd(pb1$yc8),6),round(sd(pb1$yc9),6),round(sd(pb1$yc10),6),
          round(sd(pb1$yc11),6))

mpb2 = c(round(mean(pb2$mu1),6),round(mean(pb2$mu2),6),round(mean(pb2$mu3),6),round(mean(pb2$mu4),6),round(mean(pb2$mu5),6),
         round(mean(pb2$ks1),6),round(mean(pb2$ks2),6),round(mean(pb2$ks3),6),round(mean(pb2$ks4),6),round(mean(pb2$ks5),6),
         round(mean(pb2$k1),6),round(mean(pb2$k2),6),round(mean(pb2$k3),6),
         round(mean(pb2$yc1),6),round(mean(pb2$yc2),6),round(mean(pb2$yc3),6),round(mean(pb2$yc4),6),round(mean(pb2$yc5),6),
         round(mean(pb2$yc6),6),round(mean(pb2$yc7),6),round(mean(pb2$yc8),6),round(mean(pb2$yc9),6),round(mean(pb2$yc10),6),
         round(mean(pb2$yc11),6))

sdpb2 = c(round(sd(pb2$mu1),6),round(sd(pb2$mu2),6),round(sd(pb2$mu3),6),round(sd(pb2$mu4),6),round(sd(pb2$mu5),6),
          round(sd(pb2$ks1),6),round(sd(pb2$ks2),6),round(sd(pb2$ks3),6),round(sd(pb2$ks4),6),round(sd(pb2$ks5),6),
          round(sd(pb2$k1),6),round(sd(pb2$k2),6),round(sd(pb2$k3),6),
          round(sd(pb2$yc1),6),round(sd(pb2$yc2),6),round(sd(pb2$yc3),6),round(sd(pb2$yc4),6),round(sd(pb2$yc5),6),
          round(sd(pb2$yc6),6),round(sd(pb2$yc7),6),round(sd(pb2$yc8),6),round(sd(pb2$yc9),6),round(sd(pb2$yc10),6),
          round(sd(pb2$yc11),6))

# Additional characters for LaTeX code
lp = c(' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (',' (')

rp = c(')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')',')')

# function for rounding digits and giving LaTeX format
scnt = function(x){
  y = c()
  for (i in 1:10){
    y[i] = round(x[i],digits=3)
  }
  for (i in 11:13){
    y[i] = paste('\\num{',round(x[i],digits=4),'}')
  }
  for (i in 14:24){
    y[i] = round(x[i],digits=3)
  }
  return(y)
}

camu =data.frame(scnt(mcmu),lp,scnt(sdcmu),rp)
box1=data.frame(scnt(mpb1),lp,scnt(sdpb1),rp)
box2=data.frame(scnt(mpb2),lp,scnt(sdpb2),rp)

camu$camu = paste(camu[,1],camu$lp,camu[,3],camu$rp,sep='')
camu = camu[c(-1,-2,-3,-4)]

box1$box1 = paste(box1[,1],box1$lp,box1[,3],box1$rp,sep="")
box1 = box1[c(-1,-2,-3,-4)]

box2$box2 = paste(box2[,1],box2$lp,box2[,3],box2$rp,sep="")
box2 = box2[c(-1,-2,-3,-4)]

tab = as.data.frame(cbind(par,camu,box1,box2))

# computation of |d|
tab$d = round(abs((mpb1-mpb2)/sqrt((sdpb1^2+sdpb2^2)/2)),2)

mat = xtable(tab)

# printing of LaTeX code
print(mat, sanitize.text.function = function(x) {x}, include.rownames=FALSE)
