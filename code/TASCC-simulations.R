# This code is originally from songchao1986 github:
# https://github.com/songchao1986/Nutrient-uptake

#Simulate data set. Demonstrate issues with Covino method#
library(ReacTran)
library(nlmrt)
library(minpack.lm)
library(here);i_am("code/TASCC-simulations.R")
#Advection-dispersion with transient storage and first order uptake in main channel#
NutUpTS = function (times, y, parm){
  Cs = y[1:200]
  Cts = y[201:400]
  AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
  trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
  ts.exchange = parm[5]*(Cs-Cts)
  uptake = parm[4]*Cs
  input = parm[4]*parm[1]
  dCs = trans$dC - ts.exchange - uptake + input
  dCts = parm[6]*ts.exchange
  list(c(dCs, dCts))
}

#Fit Cl and N together, assuming same residual variance#
BTCTS = function(parm){
  parm.Cl = c(parm[1],parm[3],parm[4],0,parm[5],parm[6])
  parm.N = c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6])
  yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)), rep(parm.Cl[1],200))
  yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)), rep(parm.N[1],200))
  Cl.conc = ode.1D(func=NutUpTS,y=yini.Cl, parms=parm.Cl, times=time, nspec=2)[,2*distance]
  N.conc = ode.1D(func=NutUpTS,y=yini.N, parms=parm.N, times=time, nspec=2)[,2*distance]
  return(c(Cl.conc,N.conc))
}

ResidualTS = function(parm){
  simu = BTCTS(parm)
  resid = c(data$Cl, data$N) - simu
  return(resid)
}


#Advection-dispersion with first order uptake in main channel#
NutUp = function (times, y, parm){
  Cs = y
  AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
  trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
  uptake = parm[4]*Cs
  input = parm[4]*parm[1]
  dCs = trans$dC - uptake + input
  list(c(dCs))
}


BTC = function(parm){
  parm.Cl = c(parm[1],parm[3],parm[4],0)
  parm.N = c(parm[2],parm[3],parm[4],parm[5])
  yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
  yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
  Cl.conc = ode.1D(func=NutUp,y=yini.Cl, parms=parm.Cl, times=time, nspec=1)[,2*distance]
  N.conc = ode.1D(func=NutUp,y=yini.N, parms=parm.N, times=time, nspec=1)[,2*distance]
  return(c(Cl.conc,N.conc))
}

Residual = function(parm){
  simu = BTC(parm)
  resid = c(data$Cl, data$N) - simu
  return(resid)
}

#Advection-dispersion with transient storage and M-M uptake in main channel#
NutUpMM = function (times, y, parm){
  Cs = y
  AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
  trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
  uptake = parm[4]*Cs/(Cs+parm[5])
  input = parm[4]*parm[1]/(parm[1]+parm[5])
  dCs = trans$dC - uptake + input
  list(c(dCs))
}


BTCMM = function(parm){
  parm.Cl = c(parm[1],parm[3],parm[4],0,0)
  parm.N = c(parm[2],parm[3],parm[4],parm[5],parm[6])
  yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)))
  yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)))
  Cl.conc = ode.1D(func=NutUpMM,y=yini.Cl, parms=parm.Cl, times=time, nspec=1)[,2*distance]
  N.conc = ode.1D(func=NutUpMM,y=yini.N, parms=parm.N, times=time, nspec=1)[,2*distance]
  return(c(Cl.conc,N.conc))
}

ResidualMM = function(parm){
  simu = BTCMM(parm)
  resid = c(data$Cl, data$N) - simu
  return(resid)
}


#Advection-dispersion with transient storage and M-M uptake in main channel #
NutUpTSMM = function (times, y, parm){
  Cs = y[1:200]
  Cts = y[201:400]
  AFDW = fiadeiro(v=parm[3],D=parm[2],grid=grid)
  trans = tran.1D(C=Cs, C.up=parm[1], C.down=parm[1], D=parm[2], v=parm[3], dx=grid, AFDW=AFDW)
  ts.exchange = parm[5]*(Cs-Cts)
  uptake = parm[4]*Cs/(Cs+parm[7])
  input = parm[4]*parm[1]/(parm[1]+parm[7])
  dCs = trans$dC - ts.exchange - uptake + input
  dCts = parm[6]*ts.exchange
  list(c(dCs, dCts))
}

BTCTSMM = function(parm){
  parm.Cl = c(parm[1],parm[3],parm[4],0,parm[5],parm[6],0)
  parm.N = c(parm[2],parm[3],parm[4],parm[7],parm[5],parm[6],parm[8])
  yini.Cl = c(parm.Cl[1],rep(Cl.ini,initial.length),rep(parm.Cl[1],(199-initial.length)), rep(parm.Cl[1],200))
  yini.N = c(parm.N[1],rep(N.ini,initial.length),rep(parm.N[1],(199-initial.length)), rep(parm.N[1],200))
  Cl.conc = ode.1D(func=NutUpTSMM,y=yini.Cl, parms=parm.Cl, times=time, nspec=2)[,2*distance]
  N.conc = ode.1D(func=NutUpTSMM,y=yini.N, parms=parm.N, times=time, nspec=2)[,2*distance]
  return(c(Cl.conc,N.conc))
}

ResidualTSMM = function(parm){
  simu = BTCTSMM(parm)
  resid = c(data$Cl, data$N) - simu
  return(resid)
}

#import to field data#
data.all = read.csv(here("data/test-data/LUQ13E01TPost.csv"), header=T)
width = data.all$AvgWettedWidth_m[1]
depth = data.all$AvgWettedDepth_cm[1]/100
distance = data.all$Reach.Length_meters[1]
discharge = data.all$Discharge_LitersPerSec[1]

#Calculate injected Cl in g and N in mg#
inject.Cl = data.all$Injected_NaCl_g[1]*35.5/(35.5+23) + data.all$Injected_NH4Cl_g[1]*35.5/(14+4+35.5)
inject.N = data.all$Injected_NH4Cl_g[1]*14/(14+4+35.5)*1000


#Assumping stock solution is released over one segment 0.5m#
initial.length = 1

#Injectate concentration immediately after releasing. Cl in g/m3 and N in mg/m3#
Cl.ini = inject.Cl/(width*depth*initial.length*0.5)
N.ini = inject.N/(width*depth*initial.length*0.5)

#Setting up grid, calculate time since released for each sample#
grid = setup.grid.1D(L=100, N=200)
collect.time = strptime(data.all$CollectionTime, format="%H:%M:%S")
start.time = strptime(data.all$InjectionTime[1], format="%H:%M:%S")
time = as.numeric(collect.time - start.time)
data = data.frame(Cl=data.all$ObservedCl_mgL, N=data.all$ObservedNH4N_ugL)

#Fit differet models to data, calculate AIC for each model#
mod1 = nlfb(start=c(Cl=11.7221,N=4.81,D=1.1167,U=1.107,K=0.04515), resfn=Residual, lower=c(0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf))
mod2 = nlfb(start=c(Cl=9.18,N=0.875,D=0.07537,U=1.62454,alpha=0.1423,AsA=1.69,K=0.0564), resfn=ResidualTS, lower=c(0,0,0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf,Inf,Inf))
mod3 = nlfb(start=c(Cl=11.7221,N=4.81,D=1.1167,U=1.107,Vmax=45176, Km=1000551), resfn=ResidualMM, lower=c(0,0,0,0,0,0), upper=c(12,5,Inf,Inf,Inf,Inf))
mod4 = nlfb(start=c(Cl=9.178,N=0.87,D=0.075156,U=1.6238,alpha=0.1419,AsA=1.69043,Vmax=56420,Km=1000000), resfn=ResidualTSMM, upper=c(12,5,Inf,Inf,Inf,Inf,Inf,Inf), lower=c(0,0,0,0,0,0,0,0))
