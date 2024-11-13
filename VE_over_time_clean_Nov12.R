lapply(c("emdbook","ggsci","zoo","ggplot2","ggthemes","ggnewscale","tidyverse","beepr","brglm2","dplyr","lme4","brms","rstan","mgcv","mvtnorm","scales","gtools","bbmle","estimateR"),require,character.only=T) #load packages
#setwd()

#Begin original fold change figure: Load data
FC=read.csv("fold_figure_all.csv") #load VE, nAbs data
FC$Variant[FC$Variant=="Omicron (initial)"]="BA.1 Dec. 2021"

FC2=aggregate(fold_red~Variant+Study,data=FC,FUN=mean) #average across estimates w/in a study
#remove studies w/ <3 estimates
x=as.data.frame.matrix(table(FC2$Study,FC2$Variant));
x$`BA.1 Dec. 2021`=3*x$`BA.1 Dec. 2021`
x$varN=rowSums(x);#rownames(x[x$varN>2,]) 
FC2=FC2[FC2$Study%in%rownames(x[x$varN>1,]),]
FC2$Variant=factor(FC2$Variant,levels=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"));FC2$Study=as.factor(FC2$Study) #reordering variants
FC2$figure="A";FC2$upper=NA;FC2$lower=NA

r1a=lmer(log(fold_red)~Variant+(1|Study),data = FC2);summary(r1a)
r1=lm(log(fold_red)~Variant,data = FC2);summary(r1);AIC(r1,r1a) #random effects model is preferred

#Means and SEs
mean_varO=data.frame(Variant=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"),logmean=NA,logse=NA,
                     logupper_bound=NA,loglower_bound=NA,fold_red=NA,figure="B",upper=NA,lower=NA,Study=NA)

coefs = fixef(r1a)
V = vcov(r1a)

DF.O = data.frame(Variant = factor("Alpha", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.O = model.matrix(~ Variant, data = DF.O)
mean_varO$logmean[1]=c(X.O %*% coefs)
mean_varO$logse[1]=sqrt(diag(X.O %*% V %*% t(X.O)))

DF.1 = data.frame(Variant = factor("Gamma", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.1 = model.matrix(~ Variant, data = DF.1)
mean_varO$logmean[2]=c(X.1 %*% coefs)
mean_varO$logse[2]=sqrt(diag(X.1 %*% V %*% t(X.1)))

DF.2 = data.frame(Variant = factor("Delta", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.2 = model.matrix(~ Variant, data = DF.2)
mean_varO$logmean[3]=c(X.2 %*% coefs)
mean_varO$logse[3]=sqrt(diag(X.2 %*% V %*% t(X.2)))

DF.3 = data.frame(Variant = factor("Beta", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.3 = model.matrix(~ Variant, data = DF.3)
mean_varO$logmean[4]=c(X.3 %*% coefs)
mean_varO$logse[4]=sqrt(diag(X.3 %*% V %*% t(X.3)))

DF.4 = data.frame(Variant = factor("BA.1", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.4 = model.matrix(~ Variant, data = DF.4)
mean_varO$logmean[5]=c(X.4 %*% coefs)
mean_varO$logse[5]=sqrt(diag(X.4 %*% V %*% t(X.4)))

DF.5 = data.frame(Variant = factor("BA.2", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.5 = model.matrix(~ Variant, data = DF.5)
mean_varO$logmean[6]=c(X.5 %*% coefs)
mean_varO$logse[6]=sqrt(diag(X.5 %*% V %*% t(X.5)))

DF.6 = data.frame(Variant = factor("BA.4/5", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.6 = model.matrix(~ Variant, data = DF.6)
mean_varO$logmean[7]=c(X.6 %*% coefs)
mean_varO$logse[7]=sqrt(diag(X.6 %*% V %*% t(X.6)))

DF.7 = data.frame(Variant = factor("BA.1 Dec. 2021", levels = c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021")))
X.7 = model.matrix(~ Variant, data = DF.7)
mean_varO$logmean[8]=c(X.7 %*% coefs)
mean_varO$logse[8]=sqrt(diag(X.7 %*% V %*% t(X.7)))

mean_varO$logupper_bound = (mean_varO$logmean + (1.96 * mean_varO$logse))
mean_varO$loglower_bound = (mean_varO$logmean- (1.96 * mean_varO$logse))
mean_varO$fold_red=exp(mean_varO$logmean)
mean_varO$upper=exp(mean_varO$logupper_bound)
mean_varO$lower=exp(mean_varO$loglower_bound)

FC3=mean_varO[c("Study","Variant","fold_red","upper","lower","figure")]
FC_fig=(rbind(FC2,FC3))

FC_fig$Variant=factor(FC_fig$Variant,levels=c("Alpha","Gamma","Delta","Beta","BA.1","BA.2","BA.4/5","BA.1 Dec. 2021"));FC_fig$Study=as.factor(FC_fig$Study)

#standard errors and CVs for text
(FC3$upper-FC3$fold_red)/1.96
((FC3$upper-FC3$fold_red)/1.96)/FC3$fold_red

FC3[9,]=c("","D614G",1,1,1,"B") #added for later merge
FC3$fold_red=as.numeric(FC3$fold_red);FC3$lower=as.numeric(FC3$lower);FC3$upper=as.numeric(FC3$upper);

#Pfizer waning ratios from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
NATR_pfizer_recent_two=2.3723404
NATR_pfizer_waned_two=1/(.71*497/83+.29*538/41)*NATR_pfizer_recent_two
NATR_pfizer_recent_three=(.71*1754/497+.29*1318/538)*NATR_pfizer_recent_two
NATR_pfizer_recent_three/NATR_pfizer_waned_two
NATR_hybrid_ratio=6300/800 #Goel et al
#Moderna boosting @28d relative to post 2nd dose from Chu et al medRxiv 0,181, 209: 1268, 150.2, 1951.7
NATR_moderna_recent_two=4.1332912
NATR_moderna_recent_three=1951.7/1268*NATR_moderna_recent_two
NATR_moderna_waned_two=150.2/1268*NATR_moderna_recent_two
NATR_moderna_recent_three/NATR_moderna_waned_two

#draws for NATRvar and NATRvac
ndsNATR=50000

#NATRvar distributions
NATRalpha=(rlnorm(ndsNATR,meanlog = mean_varO$logmean[mean_varO$Variant=="Alpha"],sd=mean_varO$logse[mean_varO$Variant=="Alpha"]))
NATRbeta=(rlnorm(ndsNATR,meanlog = mean_varO$logmean[mean_varO$Variant=="Beta"],sd=mean_varO$logse[mean_varO$Variant=="Beta"]))
NATRdelta=(rlnorm(ndsNATR,meanlog = mean_varO$logmean[mean_varO$Variant=="Delta"],sd=mean_varO$logse[mean_varO$Variant=="Delta"]))
NATRba.1=(rlnorm(ndsNATR,meanlog = mean_varO$logmean[mean_varO$Variant=="BA.1"],sd=mean_varO$logse[mean_varO$Variant=="BA.1"]))
NATRwt=rep(1,ndsNATR)

#NATRvac distributions
#From Khoury 2021
NATRpfi=(rlnorm(ndsNATR,meanlog = log(223/94),sd=log(10^sqrt(0.38^2/(24)+0.41^2/(38)))))
NATRmod=(rlnorm(ndsNATR,meanlog = log(654/158),sdlog =log(10^sqrt(0.21^2/(14)+0.33^2/(3)))))
NATRcon=(rlnorm(ndsNATR, meanlog = log(1),sdlog =log(10^sqrt(2*(0.61^2/(64))))))

#NATR_var
NATR_var_Alpha=1/FC3$fold_red[FC3$Variant=="Alpha"]
NATR_var_Beta=1/FC3$fold_red[FC3$Variant=="Beta"]
NATR_var_Delta=1/FC3$fold_red[FC3$Variant=="Delta"]
NATR_var_BA.1=1/FC3$fold_red[FC3$Variant=="BA.1"]
NATR_var_WT=1

## Must run main R script to obtain cvalues.csv
#load cvalues
cvalues=read.csv("cvalues.csv")

#Waning - load data
m2=read.csv("nAbs_waning.csv") #Relative (0-1) antibody titers over time

#Fit w/ nls; model for relative nAbs
g6b=nls(log2(nAbs)~((c0+c7*Infection)*exp(c1*Day+c3*Day*Moderna)  -c0-c7*Infection),#best model
        start=list(c0=.8,c1=-0.02,c3=0.01,c7=0.01),
        data=m2,trace = F,control = list(warnOnly=T));summary(g6b)#Table S4


#Load file case ratio unvax v. vax
vacc_irr0=read.csv("irr_unvax_vax2.csv")
vacc_irr=vacc_irr0[vacc_irr0$time<300,]

#GAM IRR v time
g0=gam(IRR_cases ~s(time),
       data = vacc_irr,method = "REML");summary(g0)

#Start date (first day with data) and end date (last day with data), first day vaccine given, and population, new end = end of relevant period
start_date = '2020-01-23'; end_date = '2022-01-06';first_vax="2020-12-14";population=332858873
new_end='2021-12-13'

#Fraction of mRNA vaccines given that are pfizer/moderna
frac_pfizer=(1.503/2.503);frac_moderna=1-frac_pfizer

#Build data frame using predict
last_date="2022-01-14";first_date="2021-04-04"
stop=length(seq(as.Date(last_date), as.Date(first_date), by = "-1 days"))
constant_IRR=length(seq(as.Date(first_date)-1, as.Date(first_vax), by = "-1 days"))
fitted_irr=data.frame(time=seq(0,stop-1,by=1),IRR=predict(g0,newdata=data.frame(time=seq(0,stop-1,by=1))))
fitted_irr=rbind(data.frame(time=seq(-constant_IRR,-1,by=1),IRR=fitted_irr$IRR[1]),fitted_irr)

#Reverse time order
fitted_irr = fitted_irr[order(-fitted_irr$time),]

#Remove time
case_ratios=data.frame(IRR=fitted_irr["IRR"])

#Start and end dates of case ratio data (leave as is if using CDC data)
case_ratio_start="2020-12-14";case_ratio_end="2022-01-14"
case_ratios$date=as.Date(NA)

#Create dataframe with case ratios of unvaccinated/vaccinated 
case_ratios$date=seq(as.Date(case_ratio_end), as.Date(case_ratio_start), by = "-1 days")
case_ratios_prevax = data.frame(IRR=1, date= seq(as.Date(case_ratio_start)-1, as.Date(start_date), by = "-1 days"))
case_ratios_all=rbind(case_ratios,case_ratios_prevax)

#Read in data set with daily cases, new fully vaccinated, and deaths w/ row 1 being most recent data
case_vax_data=read.csv("Cases_Vax_Death3.csv") #load case and vax data

#Add columns with date, number of infections, and proportion of 2nd doses given on a specific day
case_vax_data$date=seq(as.Date(end_date), as.Date(start_date), by = "-1 days")

#rolling mean and remove NAs
case_vax_data$daily_deaths_sm=rollmean(case_vax_data$Deaths,k=7,fill=NA,align="right")
case_vax_data=case_vax_data[complete.cases(case_vax_data$daily_deaths_sm),] #remove edge rows w/ NAs from rolling means

#Estimate infections from DEATHS by deconvolution #Distributions from Lewnard et al 2020 medRxiv
inf_shed=list(name="weibull",shape=5.983,scale=1.455) #delay infection to shedding
shed_sympt=list(name="weibull",shape=0.294,scale= 0.14) #shedding to symptom onset
sympt_hosp=list(name="gamma",shape=5.078,scale=0.765) #symptom onset to hospitalization
#hosp_death mean  13.7 (2.5-97.5%iles:1.7-34.6) #empirical distribution hospitalization to death
sh1=2.1;sc1=13.7/sh1
#qgamma(.025,shape=sh1,scale=sc1);qgamma(.975,shape=sh1,scale=sc1) #used to estimate sh1, sc1
hosp_death=list(name="gamma",shape=sh1,scale=sc1) #fitted distribution hospitalization to death
infec_deaths=deconvolve_incidence(incidence_data = case_vax_data$daily_deaths_sm,
                                  delay=list(inf_shed,shed_sympt,sympt_hosp,hosp_death))

#Daily infections
case_vax_data$inf=
  (1*infec_deaths$values[c(rep(NA,-infec_deaths$index_offset),
                           (1):(infec_deaths$index_offset+length(infec_deaths$values)))])/0.003

#Remove NAs
case_vax_data=case_vax_data[complete.cases(case_vax_data$inf),]

#Cumulative infections
case_vax_data$cum_inf=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_inf[i]=sum(case_vax_data$inf[i:nrow(case_vax_data)])
}

#Fraction total infected (332.9M US Pop 2021)
case_vax_data$frac_tot_inf=case_vax_data$cum_inf/population

#Proportion vaccinated each day
case_vax_data$prop_vax=case_vax_data$Vaccinated/sum(case_vax_data$Vaccinated)

#Create cumulative infected and cumulative fully vaccinated columns
case_vax_data$cum_vax=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_vax[i]=sum(case_vax_data$Vaccinated[i:nrow(case_vax_data)])
}

#Fraction total vaccinated (332.9 US Pop 2021)
case_vax_data$frac_tot_vax=case_vax_data$cum_vax/population

#Create cumulative uninfected and cumulative unvaccinated columns
case_vax_data$cum_unvax=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_unvax[i]=(population-case_vax_data$cum_vax[i])
}

case_vax_data$cum_uninf=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_uninf[i]=(population-case_vax_data$cum_inf[i])
}

#Add case ratios to main dataframe
case_vax_data=cbind(case_vax_data,IRR=case_ratios_all$IRR[case_ratios_all$date<=(new_end)])

#Create columns splitting infections between unvaccinated and vaccinated
case_vax_data$daily_inf_in_vax=NA
case_vax_data$daily_inf_in_unvax=NA


for (i in 1:nrow(case_vax_data)) {
  case_vax_data$daily_inf_in_vax=(case_vax_data$inf/(case_vax_data$IRR*(case_vax_data$cum_unvax/case_vax_data$cum_vax)))
}
for (i in 1:nrow(case_vax_data)) {
  case_vax_data$daily_inf_in_unvax=(case_vax_data$inf-case_vax_data$daily_inf_in_vax)
}


#Enter 0s in place of NAN
case_vax_data$daily_inf_in_vax[is.nan(case_vax_data$daily_inf_in_vax)]<-0
case_vax_data$daily_inf_in_unvax[is.nan(case_vax_data$daily_inf_in_unvax)]<-0

#Create columns splitting unvaccinated in previously positive and naive
case_vax_data$cum_unvax_neg=NA
case_vax_data$cum_unvax_pos=NA

#Prior to vaccines, split unvaccinated into previously positive and naive using infection data
length_prevax=(length(seq(as.Date(start_date), as.Date(first_vax), by = "1 days"))-1)

for (i in (nrow(case_vax_data)-length_prevax+1):nrow(case_vax_data)) {
  case_vax_data$cum_unvax_neg[i]=case_vax_data$cum_uninf[i]
  case_vax_data$cum_unvax_pos[i]=case_vax_data$cum_inf[i]
}

#After vaccines available, assume decision to get vaccinated independent of serostatus
#Account for cases in unvaccinated calculated for case ratios
for (i in 0:(nrow(case_vax_data)-length_prevax-1)) {
  case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax-i]=(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+1-i]-case_vax_data$daily_inf_in_unvax[nrow(case_vax_data)-length_prevax-i]-(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+1-i]/(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+1-i]+case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax+1-i])*case_vax_data$Vaccinated[nrow(case_vax_data)-length_prevax-i]))
  case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax-i]=(case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax+1-i]+case_vax_data$daily_inf_in_unvax[nrow(case_vax_data)-length_prevax-i]-(case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax+1-i]/(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+1-i]+case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax+1-i])*case_vax_data$Vaccinated[nrow(case_vax_data)-length_prevax-i]))
}

#Create column with the number of people gaining hybrid immunity each day
case_vax_data$daily_hybrid=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$daily_hybrid[nrow(case_vax_data)+1-i]=(case_vax_data$cum_unvax_pos[nrow(case_vax_data)+2-i]/(case_vax_data$cum_unvax_neg[nrow(case_vax_data)+2-i]+case_vax_data$cum_unvax_pos[nrow(case_vax_data)+2-i])*case_vax_data$Vaccinated[nrow(case_vax_data)+1-i]+case_vax_data$daily_inf_in_vax[nrow(case_vax_data)+1-i])
}

case_vax_data$daily_hybrid[is.na(case_vax_data$daily_hybrid)] = 0
case_vax_data$cum_hybrid=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_hybrid[i]=sum(case_vax_data$daily_hybrid[i:nrow(case_vax_data)])
}

#Create column with the daily number of unvaccinated testing positive
case_vax_data$daily_unvax_pos=NA

for (i in (nrow(case_vax_data)-length_prevax+1):nrow(case_vax_data)){
  case_vax_data$daily_unvax_pos[i]=case_vax_data$inf[i]
}

for (i in 1:(nrow(case_vax_data)-length_prevax)){
  case_vax_data$daily_unvax_pos[i]=case_vax_data$daily_inf_in_unvax[i]
}

#Create columns with adjusted number of seronegative people becoming vaccinated
case_vax_data$daily_vax_neg_adj=NA

for (i in 1:(nrow(case_vax_data)-length_prevax)) {
  case_vax_data$daily_vax_neg_adj[nrow(case_vax_data)-length_prevax+1-i]=(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+2-i]/(case_vax_data$cum_unvax_neg[nrow(case_vax_data)-length_prevax+2-i]+case_vax_data$cum_unvax_pos[nrow(case_vax_data)-length_prevax+2-i])*case_vax_data$Vaccinated[nrow(case_vax_data)-length_prevax+1-i]-case_vax_data$prop_vax[nrow(case_vax_data)-length_prevax+1-i]*sum(case_vax_data$daily_inf_in_vax))
}

#Replace NA with 0 thoughout dataframe
case_vax_data[is.na(case_vax_data)] = 0



VE_time_func=function(cutoff_df,NATR_var_XX){  
  #Cutoff 
  case_vax_data_dummy=case_vax_data[case_vax_data$date<=cutoff_df,]
  
  #Create columns with proportion of all individuals with hybrid immunity, vaccinated/seronegative, and
  #unvacccinated/seropositive (calculated two ways: 1. using case ratios and 2. using death data)
  case_vax_data_dummy$prop_hybrid=(case_vax_data_dummy$daily_hybrid/sum(case_vax_data_dummy$daily_hybrid))
  case_vax_data_dummy$prop_unvax_pos=(case_vax_data_dummy$daily_unvax_pos/sum(case_vax_data_dummy$daily_unvax_pos))
  case_vax_data_dummy$prop_vax_neg_adj=(case_vax_data_dummy$daily_vax_neg_adj/sum(case_vax_data_dummy$daily_vax_neg_adj))
  
  #Create dataframe to calculate VEs
  #Pull relevant columns from case_vax_data_dummy
  VE_data=data.frame(time=c(rep(0,15),1:nrow(case_vax_data_dummy)),prop_inf_unvax=NA,prop_hybrid=NA,prop_pfizer=NA,prop_moderna=NA)
  VE_data$prop_inf_unvax[8:(nrow(case_vax_data_dummy)+8-24)]=case_vax_data_dummy$prop_unvax_pos[24:nrow(case_vax_data_dummy)]
  VE_data$prop_hybrid[15:(nrow(case_vax_data_dummy)+15-1)]=case_vax_data_dummy$prop_hybrid
  VE_data$prop_pfizer[8:(nrow(case_vax_data_dummy)+8-1)]=case_vax_data_dummy$prop_vax_neg_adj
  VE_data$prop_moderna[1:(nrow(case_vax_data_dummy)+1-1)]=case_vax_data_dummy$prop_vax_neg_adj
  VE_data[is.na(VE_data)] = 0
  
  #Add columns with predicted nAbs for each class
  VE_data$nAbs_inf_unvax=predict(g6b,newdata=data.frame(Day=VE_data$time,Infection=1,Moderna=0))
  Day=VE_data$time
  VE_data$nAbs_inf_unvax_var=deltavar((c0+c7*1)*exp(c1*Day)  -c0-c7*1,
                                      meanval=coef(g6b),Sigma=vcov(g6b))
  VE_data$nAbs_inf_unvax_se=VE_data$nAbs_inf_unvax_var^0.5
  VE_data$nAbs_pfizer=predict(g6b,newdata=data.frame(Day=Day,Infection=0,Moderna=0))
  VE_data$nAbs_pfizer_var=deltavar((c0)*exp(c1*Day)  -c0,
                                   meanval=coef(g6b),Sigma=vcov(g6b))
  VE_data$nAbs_pfizer_se=VE_data$nAbs_pfizer_var^0.5
  VE_data$nAbs_moderna=predict(g6b,newdata=data.frame(Day=Day,Infection=0,Moderna=1))
  VE_data$nAbs_moderna_var=deltavar((c0)*exp(c1*Day+c3*Day)  -c0,
                                    meanval=coef(g6b),Sigma=vcov(g6b))
  VE_data$nAbs_moderna_se=VE_data$nAbs_moderna_var^0.5
  VE_data$nAbs_hybrid=VE_data$nAbs_pfizer*NATR_hybrid_ratio
  VE_data$nAbs_hybrid_se=VE_data$nAbs_pfizer_se
  
  
   #Create dataframe for prediction outputs from model
  y_susc_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_moderna=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_hybrid=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_inf=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  
  y_infect_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_moderna=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_hybrid=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_inf=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  
  #Create covariance matrices for calculations
  sigma_inf=diag(NA,nrow=2,ncol=2)
  sigma_trans=diag(NA,nrow=2,ncol=2)
  #sigma_sev=diag(NA,nrow=2,ncol=2)
  
  #Susceptibility and Infectiousness by day since most recent vaccination/infection by immune status for US
  #for (i in 1:1) {
  mm1=VE_data
  mean_inf=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma_inf[1,1]=cvalues$sigma11[cvalues$endpoint=="Susceptibility"];sigma_inf[1,2]=cvalues$sigma12[cvalues$endpoint=="Susceptibility"];
  sigma_inf[2,1]=cvalues$sigma21[cvalues$endpoint=="Susceptibility"];sigma_inf[2,2]=cvalues$sigma22[cvalues$endpoint=="Susceptibility"]
  x_susc=rmvnorm(n=nds, mean=mean_inf, sigma=sigma_inf)
  
  mean_trans=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma_trans[1,1]=cvalues$sigma11[cvalues$endpoint=="Infectiousness"];sigma_trans[1,2]=cvalues$sigma12[cvalues$endpoint=="Infectiousness"];
  sigma_trans[2,1]=cvalues$sigma21[cvalues$endpoint=="Infectiousness"];sigma_trans[2,2]=cvalues$sigma22[cvalues$endpoint=="Infectiousness"]
  x_infect=rmvnorm(n=nds, mean=mean_trans, sigma=sigma_trans)
  
  # 1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi[j]*NATR_var_Delta)))))
  # 1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_moderna[j],sd=mm1$nAbs_moderna_se[j]))*NATRmod[j]*NATR_var_Delta)))))
  # j=1
  for (j in 1:nrow(mm1)) {
    #for (v in 1:nds2) {
      #Susceptibility
      y_susc_pfizer[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi/NATR_var_XX)))))
      
      y_susc_moderna[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_moderna[j],sd=mm1$nAbs_moderna_se[j]))*NATRmod/NATR_var_XX)))))
      
      y_susc_hybrid[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi*NATR_hybrid_ratio/NATR_var_XX)))))
      
      y_susc_inf[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_inf_unvax[j],sd=mm1$nAbs_inf_unvax_se[j]))*NATRcon/NATR_var_XX)))))
      
      
      #Infectiousness
      y_infect_pfizer[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi/NATR_var_XX)))))
      
      y_infect_moderna[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_moderna[j],sd=mm1$nAbs_moderna_se[j]))*NATRmod/NATR_var_XX)))))
      
      y_infect_hybrid[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi*NATR_hybrid_ratio/NATR_var_XX)))))
      
      y_infect_inf[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds,mean= mm1$nAbs_inf_unvax[j],sd=mm1$nAbs_inf_unvax_se[j]))*NATRcon/NATR_var_XX)))))
      
    #}
  } 
  
  # #Susceptibility VE and data frame to calculate CIs for each immune status based on US population waning
  VE_susc_mRNA=frac_pfizer*(y_susc_pfizer%*%VE_data$prop_pfizer)+frac_moderna*(y_susc_moderna%*%VE_data$prop_moderna)
  VE_susc_hybrid=y_susc_hybrid%*%VE_data$prop_hybrid
  VE_susc_inf=y_susc_inf%*%VE_data$prop_inf_unvax

  # #Infectiousness VE and data frame to calculate CIs for each immune status based on US population waning
  VE_infect_mRNA=frac_pfizer*(y_infect_pfizer%*%VE_data$prop_pfizer)+frac_moderna*(y_infect_moderna%*%VE_data$prop_moderna)
  VE_infect_hybrid=y_infect_hybrid%*%VE_data$prop_hybrid
  VE_infect_inf=y_infect_inf%*%VE_data$prop_inf_unvax

  #95% CI and median VE over time values
  VE_time=data.frame(VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                     VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                     VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)
  
  VE_time[,c("VE_si_low","VE_si","VE_si_up")]=quantile(VE_susc_inf,probs = c(0.025,0.5,0.975))
  VE_time[,c("VE_ii_low","VE_ii","VE_ii_up")]=quantile(VE_infect_inf,probs = c(0.025,0.5,0.975))
  VE_time[,c("VE_sv_low","VE_sv","VE_sv_up")]=quantile(VE_susc_mRNA,probs = c(0.025,0.5,0.975))
  VE_time[,c("VE_iv_low","VE_iv","VE_iv_up")]=quantile(VE_infect_mRNA,probs = c(0.025,0.5,0.975))
  VE_time[,c("VE_ss_low","VE_ss","VE_ss_up")]=quantile(VE_susc_hybrid,probs = c(0.025,0.5,0.975))
  VE_time[,c("VE_is_low","VE_is","VE_is_up")]=quantile(VE_infect_hybrid,probs = c(0.025,0.5,0.975))
  
  return(list(VE_time,VE_susc_inf,VE_infect_inf,VE_susc_mRNA,VE_infect_mRNA,VE_susc_hybrid,VE_infect_hybrid))
}

#Create dataframe to store function outputs for each variant
VE_time_alpha=data.frame(date=NA,
                         VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                         VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                         VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)

VE_time_ba.1=data.frame(date=NA,
                        VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                        VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                        VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)

VE_time_delta=data.frame(date=NA,
                         VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                         VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                         VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)

VE_time_beta=data.frame(date=NA,
                        VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                        VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                        VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)

VE_time_WT=data.frame(date=NA,
                      VE_si_low=NA,VE_si=NA,VE_si_up=NA,VE_ii_low=NA,VE_ii=NA,VE_ii_up=NA,
                      VE_sv_low=NA,VE_sv=NA,VE_sv_up=NA,VE_iv_low=NA,VE_iv=NA,VE_iv_up=NA,
                      VE_ss_low=NA,VE_ss=NA,VE_ss_up=NA,VE_is_low=NA,VE_is=NA,VE_is_up=NA)

#Number of draws used in function
nds=nds2=50000

#Start/end of dates to be considered
start_cut=as.Date('2021-01-01');end_cut=as.Date('2021-12-13')
cutoff_df = seq(as.Date(start_cut),as.Date(end_cut),by = 'day')


#Alpha Variant
#Create matrices for function outputs
VE_susc_inf_alpha=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_inf_alpha=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_mRNA_alpha=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_mRNA_alpha=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_hybrid_alpha=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_hybrid_alpha=matrix(nrow=length(cutoff_df),ncol=nds)

#For loop to run function for each day
for (i in 1:length(cutoff_df)) {
  Day=c(rep(0,15),1:(nrow(case_vax_data)-length(cutoff_df)+i))
  VE_time_hold=VE_time_func(cutoff_df[i],NATRalpha)
  VE_time_alpha[i,2:19]=VE_time_hold[[1]]
  VE_susc_inf_alpha[i,1:nds]=VE_time_hold[[2]]
  VE_infect_inf_alpha[i,1:nds]=VE_time_hold[[3]]
  VE_susc_mRNA_alpha[i,1:nds]=VE_time_hold[[4]]
  VE_infect_mRNA_alpha[i,1:nds]=VE_time_hold[[5]]
  VE_susc_hybrid_alpha[i,1:nds]=VE_time_hold[[6]]
  VE_infect_hybrid_alpha[i,1:nds]=VE_time_hold[[7]]
};beep(3)

#Add date column
VE_time_alpha$date=cutoff_df

write.csv(VE_susc_inf_alpha,"VE_susc_inf_alpha.csv")
write.csv(VE_infect_inf_alpha,"VE_infect_inf_alpha.csv")
write.csv(VE_susc_mRNA_alpha,"VE_susc_mRNA_alpha.csv")
write.csv(VE_infect_mRNA_alpha,"VE_infect_mRNA_alpha.csv")
write.csv(VE_susc_hybrid_alpha,"VE_susc_hybrid_alpha.csv")
write.csv(VE_infect_hybrid_alpha,"VE_infect_hybrid_alpha.csv")

#Delta Variant
#Create matrices for function outputs
VE_susc_inf_delta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_inf_delta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_mRNA_delta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_mRNA_delta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_hybrid_delta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_hybrid_delta=matrix(nrow=length(cutoff_df),ncol=nds)

#For loop to run function for each day
for (i in 1:length(cutoff_df)) {
  Day=c(rep(0,15),1:(nrow(case_vax_data)-length(cutoff_df)+i))
  VE_time_hold=VE_time_func(cutoff_df[i],NATRdelta)
  VE_time_delta[i,2:19]=VE_time_hold[[1]]
  VE_susc_inf_delta[i,1:nds]=VE_time_hold[[2]]
  VE_infect_inf_delta[i,1:nds]=VE_time_hold[[3]]
  VE_susc_mRNA_delta[i,1:nds]=VE_time_hold[[4]]
  VE_infect_mRNA_delta[i,1:nds]=VE_time_hold[[5]]
  VE_susc_hybrid_delta[i,1:nds]=VE_time_hold[[6]]
  VE_infect_hybrid_delta[i,1:nds]=VE_time_hold[[7]]
};beep(3)

#Add date column
VE_time_delta$date=cutoff_df

write.csv(VE_susc_inf_delta,"VE_susc_inf_delta.csv")
write.csv(VE_infect_inf_delta,"VE_infect_inf_delta.csv")
write.csv(VE_susc_mRNA_delta,"VE_susc_mRNA_delta.csv")
write.csv(VE_infect_mRNA_delta,"VE_infect_mRNA_delta.csv")
write.csv(VE_susc_hybrid_delta,"VE_susc_hybrid_delta.csv")
write.csv(VE_infect_hybrid_delta,"VE_infect_hybrid_delta.csv")

#Omicron (BA.1) Variant
#Create matrices for function outputs
VE_susc_inf_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_inf_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_mRNA_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_mRNA_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_hybrid_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_hybrid_ba.1=matrix(nrow=length(cutoff_df),ncol=nds)

#For loop to run function for each day
for (i in 1:length(cutoff_df)) {
  Day=c(rep(0,15),1:(nrow(case_vax_data)-length(cutoff_df)+i))
  VE_time_hold=VE_time_func(cutoff_df[i],NATRba.1)
  VE_time_ba.1[i,2:19]=VE_time_hold[[1]]
  VE_susc_inf_ba.1[i,1:nds]=VE_time_hold[[2]]
  VE_infect_inf_ba.1[i,1:nds]=VE_time_hold[[3]]
  VE_susc_mRNA_ba.1[i,1:nds]=VE_time_hold[[4]]
  VE_infect_mRNA_ba.1[i,1:nds]=VE_time_hold[[5]]
  VE_susc_hybrid_ba.1[i,1:nds]=VE_time_hold[[6]]
  VE_infect_hybrid_ba.1[i,1:nds]=VE_time_hold[[7]]
};beep(3)

#Add date column
VE_time_ba.1$date=cutoff_df

write.csv(VE_susc_inf_ba.1,"VE_susc_inf_ba.1.csv")
write.csv(VE_infect_inf_ba.1,"VE_infect_inf_ba.1.csv")
write.csv(VE_susc_mRNA_ba.1,"VE_susc_mRNA_ba.1.csv")
write.csv(VE_infect_mRNA_ba.1,"VE_infect_mRNA_ba.1.csv")
write.csv(VE_susc_hybrid_ba.1,"VE_susc_hybrid_ba.1.csv")
write.csv(VE_infect_hybrid_ba.1,"VE_infect_hybrid_ba.1.csv")

#Beta Variant
#Create matrices for function outputs
VE_susc_inf_beta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_inf_beta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_mRNA_beta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_mRNA_beta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_hybrid_beta=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_hybrid_beta=matrix(nrow=length(cutoff_df),ncol=nds)

#For loop to run function for each day
for (i in 1:length(cutoff_df)) {
  Day=c(rep(0,15),1:(nrow(case_vax_data)-length(cutoff_df)+i))
  VE_time_hold=VE_time_func(cutoff_df[i],NATRbeta)
  VE_time_beta[i,2:19]=VE_time_hold[[1]]
  VE_susc_inf_beta[i,1:nds]=VE_time_hold[[2]]
  VE_infect_inf_beta[i,1:nds]=VE_time_hold[[3]]
  VE_susc_mRNA_beta[i,1:nds]=VE_time_hold[[4]]
  VE_infect_mRNA_beta[i,1:nds]=VE_time_hold[[5]]
  VE_susc_hybrid_beta[i,1:nds]=VE_time_hold[[6]]
  VE_infect_hybrid_beta[i,1:nds]=VE_time_hold[[7]]
};beep(3)

#Add date column
VE_time_beta$date=cutoff_df

write.csv(VE_susc_inf_beta,"VE_susc_inf_beta.csv")
write.csv(VE_infect_inf_beta,"VE_infect_inf_beta.csv")
write.csv(VE_susc_mRNA_beta,"VE_susc_mRNA_beta.csv")
write.csv(VE_infect_mRNA_beta,"VE_infect_mRNA_beta.csv")
write.csv(VE_susc_hybrid_beta,"VE_susc_hybrid_beta.csv")
write.csv(VE_infect_hybrid_beta,"VE_infect_hybrid_beta.csv")

#D614G (WT) Variant
#Create matrices for function outputs
VE_susc_inf_WT=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_inf_WT=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_mRNA_WT=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_mRNA_WT=matrix(nrow=length(cutoff_df),ncol=nds)
VE_susc_hybrid_WT=matrix(nrow=length(cutoff_df),ncol=nds)
VE_infect_hybrid_WT=matrix(nrow=length(cutoff_df),ncol=nds)

#For loop to run function for each day
for (i in 1:length(cutoff_df)) {
  Day=c(rep(0,15),1:(nrow(case_vax_data)-length(cutoff_df)+i))
  VE_time_hold=VE_time_func(cutoff_df[i],NATRwt)
  VE_time_WT[i,2:19]=VE_time_hold[[1]]
  VE_susc_inf_WT[i,1:nds]=VE_time_hold[[2]]
  VE_infect_inf_WT[i,1:nds]=VE_time_hold[[3]]
  VE_susc_mRNA_WT[i,1:nds]=VE_time_hold[[4]]
  VE_infect_mRNA_WT[i,1:nds]=VE_time_hold[[5]]
  VE_susc_hybrid_WT[i,1:nds]=VE_time_hold[[6]]
  VE_infect_hybrid_WT[i,1:nds]=VE_time_hold[[7]]
};beep(3)

#Add date column
VE_time_WT$date=cutoff_df

write.csv(VE_susc_inf_WT,"VE_susc_inf_WT.csv")
write.csv(VE_infect_inf_WT,"VE_infect_inf_WT.csv")
write.csv(VE_susc_mRNA_WT,"VE_susc_mRNA_WT.csv")
write.csv(VE_infect_mRNA_WT,"VE_infect_mRNA_WT.csv")
write.csv(VE_susc_hybrid_WT,"VE_susc_hybrid_WT.csv")
write.csv(VE_infect_hybrid_WT,"VE_infect_hybrid_WT.csv")

#Load variant frequency data and widen
var_freq=read.csv("var_freq_US.csv")
var_freq_wide=var_freq %>% 
  pivot_wider(id_cols=date,names_from = lineage,values_from = pred)

#Subset data to remove irrelevant dates
var_freq_wide$date=as.Date(var_freq_wide$date)
var_freq_wide=var_freq_wide[var_freq_wide$date>=start_cut&var_freq_wide$date<=end_cut,]

#replace NAs with 0
var_freq_wide[is.na(var_freq_wide)]=0

#Remove irrelevant dates from case_vax_data
case_vax_data$julian=julian(case_vax_data$date)
attach(case_vax_data)
case_vax_rev=case_vax_data[order(julian),]
detach(case_vax_data)
case_vax_rev=case_vax_rev[case_vax_rev$date>=start_cut&case_vax_rev$date<=end_cut,]

#Covariance matrices
sigma_inf=diag(NA,nrow=2,ncol=2)
sigma_trans=diag(NA,nrow=2,ncol=2)

#Draws for coefficients for VE v NATR relationships
mean_inf=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
sigma_inf[1,1]=cvalues$sigma11[cvalues$endpoint=="Susceptibility"];sigma_inf[1,2]=cvalues$sigma12[cvalues$endpoint=="Susceptibility"];
sigma_inf[2,1]=cvalues$sigma21[cvalues$endpoint=="Susceptibility"];sigma_inf[2,2]=cvalues$sigma22[cvalues$endpoint=="Susceptibility"]
x_susc=rmvnorm(n=nds, mean=mean_inf, sigma=sigma_inf)

mean_trans=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
sigma_trans[1,1]=cvalues$sigma11[cvalues$endpoint=="Infectiousness"];sigma_trans[1,2]=cvalues$sigma12[cvalues$endpoint=="Infectiousness"];
sigma_trans[2,1]=cvalues$sigma21[cvalues$endpoint=="Infectiousness"];sigma_trans[2,2]=cvalues$sigma22[cvalues$endpoint=="Infectiousness"]
x_infect=rmvnorm(n=nds, mean=mean_trans, sigma=sigma_trans)

#VE for boosted mRNA against susceptibility and infectiousness and repeat rows 
VE_susc_boost_mRNA=frac_pfizer*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))+frac_moderna*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))
# Replicate each row
VE_sb <- do.call("rbind", replicate( 
  length(cutoff_df), VE_susc_boost_mRNA, simplify = FALSE))

VE_infect_boost_mRNA=frac_pfizer*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))+frac_moderna*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))

# Replicate each row
VE_ib <- do.call("rbind", replicate( 
  length(cutoff_df), VE_infect_boost_mRNA, simplify = FALSE))


#Weighted VEs by day based on variant frequenct
VE_ih=var_freq_wide$Alpha*VE_infect_hybrid_alpha+var_freq_wide$Delta*VE_infect_hybrid_delta+
  var_freq_wide$Beta*VE_infect_hybrid_beta+var_freq_wide$Other*VE_infect_hybrid_WT+
  var_freq_wide$'Om:BA.1'*VE_infect_hybrid_ba.1

VE_ii=var_freq_wide$Alpha*VE_infect_inf_alpha+var_freq_wide$Delta*VE_infect_inf_delta+
  var_freq_wide$Beta*VE_infect_inf_beta+var_freq_wide$Other*VE_infect_inf_WT+
  var_freq_wide$'Om:BA.1'*VE_infect_inf_ba.1

VE_iv=var_freq_wide$Alpha*VE_infect_mRNA_alpha+var_freq_wide$Delta*VE_infect_mRNA_delta+
  var_freq_wide$Beta*VE_infect_mRNA_beta+var_freq_wide$Other*VE_infect_mRNA_WT+
  var_freq_wide$'Om:BA.1'*VE_infect_mRNA_ba.1

VE_sh=var_freq_wide$Alpha*VE_susc_hybrid_alpha+var_freq_wide$Delta*VE_susc_hybrid_delta+
  var_freq_wide$Beta*VE_susc_hybrid_beta+var_freq_wide$Other*VE_susc_hybrid_WT+
  var_freq_wide$'Om:BA.1'*VE_susc_hybrid_ba.1

VE_si=var_freq_wide$Alpha*VE_susc_inf_alpha+var_freq_wide$Delta*VE_susc_inf_delta+
  var_freq_wide$Beta*VE_susc_inf_beta+var_freq_wide$Other*VE_susc_inf_WT+
  var_freq_wide$'Om:BA.1'*VE_susc_inf_ba.1

VE_sv=var_freq_wide$Alpha*VE_susc_mRNA_alpha+var_freq_wide$Delta*VE_susc_mRNA_delta+
  var_freq_wide$Beta*VE_susc_mRNA_beta+var_freq_wide$Other*VE_susc_mRNA_WT+
  var_freq_wide$'Om:BA.1'*VE_susc_mRNA_ba.1

#Fraction unvaccinated; Fraction unvaccinated among previously infected
f_n=1-case_vax_rev$frac_tot_vax;f_pu=case_vax_rev$cum_unvax_pos/(case_vax_rev$cum_unvax_pos+case_vax_rev$cum_unvax_neg)
#Fraction vaccinated among previously vaccinated; Fraction twice vaccinated; Fraction boosted
f_pv=case_vax_rev$cum_hybrid/case_vax_rev$cum_vax;f_v=case_vax_rev$frac_tot_vax;f_b=case_vax_rev$Frac_boosted

#Population VE weighted by fraction in each group
VE_i_pop=f_n*f_pu*VE_ii+(1-f_pv)*(f_v-f_b)*VE_iv+(f_pv)*(f_v-f_b)*VE_ih+(f_b)*VE_ib
VE_s_pop=f_n*f_pu*VE_si+(1-f_pv)*(f_v-f_b)*VE_sv+(f_pv)*(f_v-f_b)*VE_sh+(f_b)*VE_sb

#Create dataframes to store values
frac_inf=data.frame(lower=NA,mean=rep(f_n*f_pu,2),upper=NA,VE=NA,date=rep(cutoff_df,2),plot="A: Population Immunity",Immunity="Frac. Inf.")
frac_vax=data.frame(lower=NA,mean=rep((1-f_pv)*(f_v-f_b),2),upper=NA,VE=NA,date=rep(cutoff_df,2),plot="A: Population Immunity",Immunity="Frac. Vacc.")
frac_hybrid=data.frame(lower=NA,mean=rep((f_pv)*(f_v-f_b),2),upper=NA,VE=NA,date=rep(cutoff_df,2),plot="A: Population Immunity",Immunity="Frac. Hybrid")
frac_boost=data.frame(lower=NA,mean=rep(f_b,2),upper=NA,VE=NA,date=rep(cutoff_df,2),plot="A: Population Immunity",Immunity="Frac. Boost")

VE_s_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_si_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_sv_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_sh_CI=data.frame(lower=NA,mean=NA,upper=NA)

VE_i_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_ii_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_iv_CI=data.frame(lower=NA,mean=NA,upper=NA)
VE_ih_CI=data.frame(lower=NA,mean=NA,upper=NA)

#Cycle through VE values each day and take 95% CIs and medians
for (i in 1:length(cutoff_df)) {
  VE_s_CI[i,]=quantile(VE_s_pop[i,],probs=c(0.025,0.5,0.975))
  VE_si_CI[i,]=quantile(VE_si[i,],probs=c(0.025,0.5,0.975))
  VE_sv_CI[i,]=quantile(VE_sv[i,],probs=c(0.025,0.5,0.975))
  VE_sh_CI[i,]=quantile(VE_sh[i,],probs=c(0.025,0.5,0.975))
  
  VE_i_CI[i,]=quantile(VE_i_pop[i,],probs=c(0.025,0.5,0.975))
  VE_ii_CI[i,]=quantile(VE_ii[i,],probs=c(0.025,0.5,0.975))
  VE_iv_CI[i,]=quantile(VE_iv[i,],probs=c(0.025,0.5,0.975))
  VE_ih_CI[i,]=quantile(VE_ih[i,],probs=c(0.025,0.5,0.975))
}

#Plot legend
VE_s_CI$VE=VE_i_CI$VE="Whole pop.";VE_si_CI$VE=VE_ii_CI$VE="Prev. Inf."
VE_sv_CI$VE=VE_iv_CI$VE="Vacc.";VE_sh_CI$VE=VE_ih_CI$VE="Hybrid"

#Add date columns
VE_s_CI$date=VE_i_CI$date=VE_si_CI$date=VE_ii_CI$date=VE_sv_CI$date=VE_iv_CI$date=
  VE_sh_CI$date=VE_ih_CI$date=cutoff_df

#Bind into VE_s and VE_i dataframes
VE_s_CI_all=rbind(VE_s_CI,VE_si_CI,VE_sv_CI,VE_sh_CI)
VE_i_CI_all=rbind(VE_i_CI,VE_ii_CI,VE_iv_CI,VE_ih_CI)

#Facet label
VE_s_CI_all$plot="B: Susceptibility"
VE_i_CI_all$plot="C: Infectiousness"
VE_s_CI_all$Immunity=NA;VE_i_CI_all$Immunity=NA

#Create final dataframes for Figure 3A and 3BC
FracImmune=rbind(frac_inf,frac_vax,frac_hybrid,frac_boost)
Figure3_df=rbind(VE_s_CI_all,VE_i_CI_all)

#Write csvs for use in main R script
write.csv(FracImmune,"FracImmune.csv")
write.csv(Figure3_df,"Figure3ABC2.csv")

