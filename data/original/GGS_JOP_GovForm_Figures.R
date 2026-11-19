library(lattice)
library(mgcv)
library(np)
library(survival)
library(sandwich)
library(rms)

#Figure 1: single run
onerun<-read.csv("GGS_JOP_GovForm_Data_Fig1.csv")
x11()
p1<-xyplot(Running.Slope.of.Cabinet.Share.with.Seat.Share.Contributed~tick,data=onerun,par.settings=list(par.main.text = list(cex = 0.8)),
xlab="Government Formations",ylab="Slope", main="Cumulative Slope Coefficient",
panel = function(x, y, subscripts, groups, ...) {
	     panel.xyplot(x,y,...)
           panel.abline(v = 2500,lty=2)
       })
p2<-xyplot(Running.Intercept.of.Cabinet.Share.with.Seat.Share.Contributed~tick,data=onerun,par.settings=list(par.main.text = list(cex = 0.8)),
xlab="Government Formations",ylab="Intercept", main="Cumulative Intercept Coefficient",
panel = function(x, y, subscripts, groups, ...) {
	     panel.xyplot(x,y,...)
           panel.abline(v = 2500,lty=2)
       })
print(p1,position=c(0,0,0.5,1),more=T)
print(p2,position=c(0.5,0,1,1))
rm(onerun)


#Loads government formation batch data and changes some names for convenience
#Data for Figures 2 and 3 and Tables 1-4
quad2d<-read.csv("GGS_JOP_GovForm_Data.csv")
names(quad2d)[74]<-"Prop_Single.Party_Majority"
names(quad2d)[75]<-"Prop_Single.Party_Minority"
names(quad2d)[76]<-"Prop_Coalition_Minority"
names(quad2d)[77]<-"Prop_Minimal_Winning"
names(quad2d)[78]<-"Prop_Surplus_Simple"
names(quad2d)[79]<-"Prop_Unity"
names(quad2d)[2]<-"beta_sim"
names(quad2d)[24]<-"enpp"
names(quad2d)[25]<-"diversity"
names(quad2d)[66]<-"delay"
attach(quad2d)

#Fix Some Variables
#New surplus variables that includes unity governments
Prop_Surplus<-quad2d$Prop_Surplus_Simple+quad2d$Prop_Unity
quad2d<-cbind(quad2d,Prop_Surplus)
psg<-Proportion.Surplus.Govt+Proportion.Unity.Govt
quad2d<-cbind(quad2d,psg)
rm(Prop_Surplus)
rm(psg)

#Calc real-world weights
real1<-read.csv("GGS_JOP_GovForm_Data_real1.csv")
bw<-npudensbw(formula=~enpp,dat=real1)
enpp_w<-npudens(bws=bw,newdata=quad2d)
enpp_weights<-enpp_w$dens
quad2d<-cbind(quad2d,enpp_weights)
rm(enpp_weights)
attach(quad2d)

#Table 1: government type
summary(Proportion.Single.Party.Majority.Govt)
summary(Proportion.Single.Party.Minority.Govt)
summary(Proportion.Coalition.Minority.Govt)
summary(Proportion.Minimal.Winning.Govt)
summary(psg)
weighted.mean(Proportion.Single.Party.Majority.Govt,enpp_weights,na.rm=TRUE)
weighted.mean(Proportion.Single.Party.Minority.Govt,enpp_weights,na.rm=TRUE)
weighted.mean(Proportion.Coalition.Minority.Govt,enpp_weights,na.rm=TRUE)
weighted.mean(Proportion.Minimal.Winning.Govt,enpp_weights,na.rm=TRUE)
weighted.mean(psg,enpp_weights,na.rm=TRUE)

#Figure 2: government type as a function of beta, parties' policy preference
x11()
op<-par(pch=46,mfrow=c(2,2))
plot(gam(Proportion.Single.Party.Minority.Govt~s(beta_sim),weights=enpp_weights),rug=FALSE,xlab=expression(beta),ylim=c(-.1,.1),
ylab="Proportion Single-Party Minority",main="(A) Single-Party Minority")
abline(h=0,col=2)
plot(gam(Proportion.Coalition.Minority.Govt~s(beta_sim),weights=enpp_weights),rug=FALSE,xlab=expression(beta),ylim=c(-.1,.1),
ylab="Proportion Coalition Minority",main="(B) Coalition Minority")
abline(h=0,col=2)
plot(gam(Proportion.Minimal.Winning.Govt~s(beta_sim),weights=enpp_weights),rug=FALSE,xlab=expression(beta),ylim=c(-.1,.1),
ylab="Proportion Minimal Winning Coalition",main="(C) Minimal Winning Coalition")
abline(h=0,col=2)
plot(gam(psg~s(beta_sim),weights=enpp_weights),rug=FALSE,xlab=expression(beta),ylim=c(-.1,.1),
ylab="Proportion Surplus Government",main="(D) Surplus Government")
abline(h=0,col=2)
par(op)

#Table 2: Gamson's Law regressions
formlist<-mat.or.vec(500000,1)
formcssfl<-mat.or.vec(500000,1)
formcsfl<-mat.or.vec(500000,1)
wghts<-mat.or.vec(500000,1)
p_spmaj<-mat.or.vec(500000,1)
p_spmin<-mat.or.vec(500000,1)
p_mwc<-mat.or.vec(500000,1)
p_mc<-mat.or.vec(500000,1)
p_surplus<-mat.or.vec(500000,1)
bfl<-mat.or.vec(500000,1)
divfl<-mat.or.vec(500000,1)
enppfl<-mat.or.vec(500000,1)
numpfl<-mat.or.vec(500000,1)
govdummy<-mat.or.vec(500000,1)
dly<-mat.or.vec(500000,1)
cnter<-0
for(i in 1:50000) {
	plc<-numParties[i]-1
	for(j in 0:plc) {
    govdummy[cnter+j+1]<-i
		formcssfl[cnter+j+1]<-quad2d[i,46+j]
		formcsfl[cnter+j+1]<-quad2d[i,56+j]
		wghts[cnter+j+1]<-enpp_weights[i]
		dly[cnter+j+1]<-delay[i]
		p_spmaj[cnter+j+1]<-Proportion.Single.Party.Majority.Govt[i]
		p_spmin[cnter+j+1]<-Proportion.Single.Party.Minority.Govt[i]
		p_mwc[cnter+j+1]<-Proportion.Minimal.Winning.Govt[i]
		p_mc[cnter+j+1]<-Proportion.Coalition.Minority.Govt[i]
		p_surplus[cnter+j+1]<-Proportion.Surplus.Govt[i]+Proportion.Unity.Govt[i]
		bfl[cnter+j+1]<-beta_sim[i]
		divfl[cnter+j+1]<-diversity[i]
		enppfl[cnter+j+1]<-enpp[i]
		numpfl[cnter+j+1]<-numParties[i]
	}
	formlist[cnter+Mean.Incumbent[i]+1]=1
	cnter<-cnter+10
}
slopecalc<-data.frame(cbind(formlist,formcssfl,formcsfl,wghts,p_spmaj,p_spmin,p_mwc,p_mc,p_surplus,bfl,divfl,enppfl,numpfl,govdummy,dly))
scalc<-subset(slopecalc,formcsfl>0.0)
with(scalc,summary(lm(formcsfl~formcssfl)))
with(scalc,summary(lm(formcsfl~formcssfl,weights=wghts)))
#With clustered standard errors
gl1uw <- ols(formcsfl~formcssfl, x=TRUE, y=TRUE, data=scalc)
gl1uwcse<-robcov(gl1uw, scalc$govdummy)
gl1uw
gl1uwcse
gl1 <- ols(formcsfl~formcssfl, x=TRUE, y=TRUE, data=scalc, weights=wghts)
gl1cse<-robcov(gl1, scalc$govdummy)
gl1
gl1cse

#Table 3: Effect of Parameters on Gamson's Law
with(scalc,summary(lm(formcsfl~formcssfl+formlist,weights=wghts)))
with(scalc,summary(lm(formcsfl~formcssfl*bfl+formcssfl*numpfl+formcssfl*enppfl+formcssfl*divfl,weights=wghts)))
#With clustered standard errors
gl2<-ols(formcsfl~formcssfl+formlist, x=TRUE, y=TRUE, data=scalc, weights=wghts)
robcov(gl2, scalc$govdummy)
gl3<-ols(formcsfl~formcssfl*bfl+formcssfl*numpfl+formcssfl*enppfl+formcssfl*divfl, x=TRUE, y=TRUE, data=scalc, weights=wghts)
robcov(gl3, scalc$govdummy)

#Summary of delay data (in text)
summary(delay)
#Weighted by effective number of parties in real-world data
weighted.mean(delay,enpp_weights,na.rm=TRUE)

#Figure 3: bargaining delay

#Weights for density plots must be normalized
enpp_w_2<-enpp_weights/sum(enpp_weights)

x11()
h1<-densityplot(~bargaining_duration,data=real1,par.settings=list(par.main.text = list(cex = 0.8)),xlab="Bargaining Delay",
main="Real World Data",plot.points=F)
h2<-densityplot(~2.08*delay,data=quad2d,par.settings=list(par.main.text = list(cex = 0.8)),weights=enpp_w_2,xlab="Bargaining Delay",
main="Weighted Simulation Data",plot.points=F)
print(h1,position=c(0,0,0.5,1),more=T)
print(h2,position=c(0.5,0,1,1))

#Table 4: bargaining delay regressions
summary(lm(delay~beta_sim+numParties+enpp+diversity,weights=enpp_weights))
delay2<-1:50000
for(i in 1:50000) {
	if(delay[i]==100) delay2[i]<-0 else delay2[i]<-1
}
quad2d<-cbind(quad2d,delay2)
srv = Surv(delay,delay2)
cx2<-coxph(srv~beta_sim+numParties+enpp+diversity,data=quad2d,weights=enpp_weights)
summary(cx2)