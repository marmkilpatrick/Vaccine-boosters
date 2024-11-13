lapply(c("emdbook","ggblend","ggpubr","ggsci","cowplot","zoo","ggplot2","ggthemes","tidyverse","brglm2","dplyr","lme4","brms","rstan","mgcv","mvtnorm","scales","gtools","bbmle","ggnewscale","estimateR","beepr","egg"),require,character.only=T) #load packages
#setwd()

mAll=read.csv("VE_nAbs_data7.csv") #load VE, nAbs data
mAll=mAll[mAll$endpoint!="Symptomatic disease",] #Remove Symptomatic disease data
mAll=mAll[mAll$Variant!="BA.2",] #Remove BA.2 data


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

#standard errors and CVs
(FC3$upper-FC3$fold_red)/1.96
((FC3$upper-FC3$fold_red)/1.96)/FC3$fold_red

FC3[9,]=c("","D614G",1,1,1,"B") #added for later merge
FC3$fold_red=as.numeric(FC3$fold_red);FC3$lower=as.numeric(FC3$lower);FC3$upper=as.numeric(FC3$upper);

#Merge fold change data with main dataframe
mAll=merge(mAll,FC3[,c("Variant","fold_red")],by="Variant") #merging VE data w/ Nabs data

#add column for NATR_tot
mAll$NATR_tot=mAll$NATR_vac*mAll$NATR_vac_istatus/mAll$fold_red

#add vaccine codes to main dataframe
vac_comp=data.frame(Vaccine=c("Astrazeneca","Convalescent","Johnson and Johnson","Moderna","Novavax","Pfizer","Sinovac","Sputnik"),
                    code=c("ChAdOx1 nCoV-19","Convalescent","Ad26.COV2.S","mRNA-1273","NVX-CoV2373","BNT162b2","CoronaVac","Sputnik V"))

mAll=merge(mAll,vac_comp,by="Vaccine")

#NATR adjustments (NATR_i) for various vaccinations states
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

#Fraction of mRNA vaccines that are Pfizer/Moderna in the US
mRNA_pfizer_ratio=(1.503/2.503)
mRNA_moderna_ratio=(1/2.503)

#Draws for NATR_var for Delta variant
nds2=10000 #number of draws (same for NATR_vac)
NATRdelta=(rlnorm(nds2,meanlog = mean_varO$logmean[mean_varO$Variant=="Delta"],sd=mean_varO$logse[mean_varO$Variant=="Delta"]))
quantile(NATRdelta,probs=c(0.025,0.975))

#Draws for NATR_vac
#From Khoury 2021
NATRpfi=(rlnorm(nds2,meanlog = log(223/94),sd=log(10^sqrt(0.38^2/(24)+0.41^2/(38)))))
quantile(NATRpfi,probs=c(0.025,0.975))
NATRmod=(rlnorm(nds2,meanlog = log(654/158),sdlog =log(10^sqrt(0.21^2/(14)+0.33^2/(3)))))
quantile(NATRmod,probs=c(0.025,0.975))
NATRcon=(rlnorm(nds2, meanlog = log(1),sdlog =log(10^sqrt(2*(0.61^2/(64))))))
quantile(NATRcon,probs=c(0.025,0.975))

#Delta NATR_var (mean)
NATR_var_Delta=1/FC3$fold_red[FC3$Variant=="Delta"]

#function to calculate effective events in each arm
VE_ci=function (par,mrow) {
  #  print(mrow)
  if(mrow$VE>=1){
    E_rat=exp(par[1])
    E1s=0;N1s=100000000;N2s=100000000;E2s=100
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1-mean(xnas1)
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(mrow$upper-VE$VE_up)^2+(mrow$lower-VE$VE_low)^2
    print(paste(c(E_rat,VE,err),sep=" "))
    err
    
    
  }
  else{
    E_rat=exp(par[1])
    E1s=100;N1s=100000000;N2s=100000000;E2s=(E1s/N1s)*N2s/(1-mrow$VE)
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    VE=1-(E1/N1)/(E2/N2)
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N);#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    VE=data.frame(VE_m=NA,VE_up=NA,VE_low=NA)
    VE$VE_m=1-mean(xnas1)
    VE$VE_up=1-quantile(xnas1,0.025)
    VE$VE_low=1-quantile(xnas1,0.975)
    err=(mrow$upper-VE$VE_up)^2+(mrow$lower-VE$VE_low)^2
    #  print(paste(c(E_rat,VE,err),sep=" "))
    err
  }
}

#calculate effective events in each arm
for (i in 1:nrow(mAll)) {
  mrow=mAll[i,c("Study","VE","upper","lower")]
  print(mrow)
  if(mrow$VE>=1){
    o1=optim(par=c(1),fn=VE_ci,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",mrow=mrow);o1
    E1s=0;N1s=100000000;N2s=100000000;E2s=100
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1
    mAll$I_v[i]=E1;mAll$I_c[i]=E2;mAll$N_v[i]=N1;mAll$N_c[i]=N2
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N,method = "brglmFit");#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    mAll$VE_m[i]=1-mean(xnas1)
    mAll$VE_up[i]=1-quantile(xnas1,0.025)
    mAll$VE_low[i]=1-quantile(xnas1,0.975)
  }else{
    o1=optim(par=c(1),fn=VE_ci,control=list(trace=T,reltol=1e-4),method="Nelder-Mead",mrow=mrow);o1
    E1s=100;N1s=100000000;N2s=100000000;E2s=(E1s/N1s)*N2s/(1-mrow$VE)
    E_rat=exp(o1$par[1])
    E1=E1s*E_rat;E2=E2s*E_rat;N1=N1s;N2=N2s
    1-(E1/N1)/(E2/N2)
    mAll$I_v[i]=E1;mAll$I_c[i]=E2;mAll$N_v[i]=N1;mAll$N_c[i]=N2
    m1=data.frame(Prev=c(E1/N1,E2/N2),N=c(N1,N2),Group=c("Vaccine","Control"))
    f2=glm(Prev~Group,data=m1,family=binomial,weights=N);#summary(f2)
    m1$fit=predict(f2,type="link",se.fit=T)$fit
    m1$se.fit=predict(f2,type="link",se.fit=T)$se.fit
    xnas1=plogis(rnorm(50000,m1$fit[m1$Group=="Vaccine"],
                       m1$se.fit[m1$Group=="Vaccine"]))/
      plogis(rnorm(50000,m1$fit[m1$Group=="Control"],
                   m1$se.fit[m1$Group=="Control"]))
    mAll$VE_m[i]=1-mean(xnas1)
    mAll$VE_up[i]=1-quantile(xnas1,0.025)
    mAll$VE_low[i]=1-quantile(xnas1,0.975)
  }
}

#Cap events and add column
cap=1000
mAll$I_v[mAll$I_c>cap]=mAll$I_v[mAll$I_c>cap]*cap/mAll$I_c[mAll$I_c>cap]
mAll$I_c[mAll$I_c>cap]=cap

#prevalence
mAll$b_s=mAll$I_c/mAll$N_c

#Likelihood for estimating coefficients for relationship between VE and nAbs
lkf=function(c0,c1) {   #This function estimates the neg log likelihood using data in data.frame m1
  llk=matrix(0,(nrow(mm1)),1)
  for (i in 1:nrow(mm1)) {
    
    llk[i]=mm1$I_c[i]*log(mm1$b_s[i])+(mm1$N_c[i]-mm1$I_c[i])*log(1-mm1$b_s[i])+(mm1$I_v[i])*log(mm1$b_s[i]*((1/(1+exp(-(c0+c1*log2(mm1$NATR_tot[i])))))))+(mm1$N_v[i]-mm1$I_v[i])*log(1-mm1$b_s[i]*((1/(1+exp(-(c0+c1*log2(mm1$NATR_tot[i])))))))
    
  }
  tkl=-sum(llk) #sum likelihood 
  return(tkl)	}

#number of draws 
nds=10000
qdraws=runif(nds)

#create dataframes for model outputs
y=matrix(data=NA,nrow =nds, ncol = 21 )
y0=matrix(data=NA,nrow =nds, ncol = 21 )
z=data.frame(endpoint=NA,NATR_tot=rep(0.0078125*2^((0:20)/2),3),mean=NA)
CI=data.frame(endpoint=NA,NATR_tot=rep(0.0078125*2^((0:20)/2),3),mean=NA)
VE_l=matrix(data=NA,nrow =nds, ncol = 21 )

#baseline prev from original UK report Andrews et al.
base_prev=115/10000

#Calculate coefficient values for model for each endpoint, fitted line/CIs
ep=unique(mAll$endpoint)

cvalues=data.frame(endpoint=c("Infectiousness","Susceptibility"),c0=NA,c1=NA,sigma11=NA,sigma12=NA,sigma21=NA,sigma22=NA)

for (i in 1:length(ep)) {
  mm1=mAll[mAll$endpoint==ep[i],]
  fp2 <- mle2(lkf,start=list(c0=0,c1=0),
              fixed=list(),control=list(trace=3))
  
  h1=fp2@details$hessian
  mean=c(coef(fp2)[1],coef(fp2)[2])
  sigma=solve(h1)
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  cvalues$endpoint[i]=ep[i]
  cvalues$c0[i]=coef(fp2)[1]
  cvalues$c1[i]=coef(fp2)[2]
  cvalues$sigma11[i]=sigma[1,1]
  cvalues$sigma12[i]=sigma[1,2]
  cvalues$sigma21[i]=sigma[2,1]
  cvalues$sigma22[i]=sigma[2,2]
  if(i<2)
  {susc_sum=summary(fp2)}  
  else
  {infect_sum=summary(fp2);
  y0=y}
  
  for (j in 0:20) {
    y[,j+1]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(0.0078125*2^(j/2))))))
  }
  for (k in 1:ncol(y)) {
    CI$endpoint[k+ncol(y)*(i-1)]=ep[i]
    CI$lower[k+ncol(y)*(i-1)]=quantile(y[,k],probs=0.025)
    CI$upper[k+ncol(y)*(i-1)]=quantile(y[,k],probs=0.975)
  }
  for(l in 1:nrow(y)){
    for(n in 1:ncol(y)){
      VE_l[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y[l,n]))/rbinom(1,10000,prob=base_prev)
    }
  }
  for (p in 1:ncol(y)) {
    z$endpoint[p+ncol(y)*(i-1)]=ep[i]
    z$mean[p+ncol(y)*(i-1)]=quantile(y[,p],probs=0.5)
    z$lower[p+ncol(y)*(i-1)]=quantile(VE_l[,p],probs=0.025)
    z$upper[p+ncol(y)*(i-1)]=quantile(VE_l[,p],probs=0.975)
  }
}

#Create cvalues.csv for use in VE_over_time R script
write.csv(cvalues,file = "cvalues.csv")

#Combine draws for VE for transmission calc. VEt = (1-VEs)*(1-VEi)
yt=1-(1-y) * (1-y0)

#Modify dataframe with fitted line and PIs
z$endpoint[43:63]="Transmission"
z$lower[z$endpoint=="Transmission"]=NA
z$upper[z$endpoint=="Transmission"]=NA

#Modify dataframe with CIs
CI$endpoint[43:63]="Transmission"
CI$lower[CI$endpoint=="Transmission"]=NA
CI$upper[CI$endpoint=="Transmission"]=NA

#Calculate fitted line and CIs for VE for transmission
for (i in 1:ncol(yt)) {
  z$mean[42+i]=quantile(yt[,i],probs=0.5)
  CI$lower[42+i]=quantile(yt[,i],probs=0.025)
  CI$upper[42+i]=quantile(yt[,i],probs=0.975)
}

#Remove unused rows for plotting
z=z[-(43:50), ]
CI=CI[-(43:50), ]



#Label for fig 1
tlab2=data.frame(endpoint=c("Susceptibility","Infectiousness","Transmission"),text1=c("A: Susceptibility","B: Infectiousness","C: Transmission"),NATR_tot=c(0.13,0.13,0.13),VE=c(0.93,0.9,0.95))

#jitter x-values
mAll$NATR_tot_jit=mAll$NATR_tot+rnorm(nrow(mAll),0,0.02) #jitter xvals

#Split data by type for ggplot2
z$data="mean"
CI$data="CI"
mAll$data="data"


#modify data to fit plot
mAll_fig2=mAll
mAll_fig2$lower[mAll_fig2$lower<=-0.1&mAll_fig2$data=="data"]=-0.1
mAll_fig2=mAll_fig2[mAll_fig2$VE>=0.1&mAll_fig2$data=="data",]

#Data into single dataframe for plot
fig2=bind_rows(mAll_fig2,CI,z)#,pred_delta)

#Reorder for plotting
fig2$endpoint=factor(fig2$endpoint,levels = c("Susceptibility", "Infectiousness","Transmission"))
tlab2$endpoint=factor(tlab2$endpoint,levels = c("Susceptibility","Infectiousness","Transmission"))

#Add "sera" to Convalescent
fig2$code[fig2$code=="Convalescent sera"]="Convalescent"

#Plot Figure 1
ggplot(data=fig2)+
  geom_ribbon(data=fig2[fig2$data=="CI",],aes(x=NATR_tot,ymin=lower, ymax=upper),fill="gray",alpha=0.5)+
  geom_line(data=fig2[fig2$data=="mean",],aes(x=NATR_tot,y=mean),color="black")+
  geom_point(data=fig2[fig2$data=="data",],aes(x=NATR_tot_jit,y=VE,col=Variant,shape=code),size=3)+
  geom_errorbar(data=fig2[fig2$data=="data",],aes(x=NATR_tot_jit,ymin=lower, ymax=upper,col=Variant), width=0)+
  #geom_ribbon(data=fig2[fig2$data=="mean",],aes(x=NATR_tot,ymin=lower, ymax=upper),fill="gray",alpha=0.5)+
  #scale_x_continuous(trans='log2',limits=c(0.0625,4.1),n.breaks=7)+
  theme_few()+
  scale_x_continuous(expand = c(0,0),trans="log2",limits=c(0.125,4),breaks=c(0.0078125*2^c(4:9)),labels = label_number(accuracy = 0.001, scale = 1, 
                                                                                                                       prefix = "", suffix = "",
                                                                                                                       big.mark = " ", decimal.mark = "."))+
  scale_y_continuous(expand = c(0, 0))+
  scale_shape_manual(values = c(18,19,17,6,15))+ 
  coord_cartesian(ylim = c(NA, 1),xlim = c(NA,NA)) +
  #scale_color_manual(values=c(  "#00BFC4", "#C77CFF","#F8766D"))+
  xlab(expression("Neutralizing antibody titer ratio ("*NATR[tot]*")"))+
  ylab("VE for susceptibility or infectiousness")+
  theme(axis.title=element_text(size=24),#legend.position = "top",#legend.position = c(.1, .80),
        axis.text=element_text(size=18),legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  #geom_text(data=tlab,aes(x=Nabs_Ratio,y=VE,label=text1))+
  geom_text(data=tlab2,aes(x=NATR_tot,y=VE,label=text1),size=6.5,hjust=0)+
  labs(col="Variant", shape="Immunity")+
  facet_wrap(.~endpoint,nrow=3, scales="free_y")+theme(strip.text.x = element_blank())

#Supplemental Figure 1
Susc_Symp_Rat=data.frame(NATR_tot=seq(0.125,4,by=0.0005))

#Relative risk susceptibility
y1_susc=matrix(data=NA,nrow =nds/10, ncol = nrow(Susc_Symp_Rat) )
for (j in 1:nrow(Susc_Symp_Rat)) {
  mean=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma[1,1]=cvalues[1,4];sigma[1,2]=cvalues[1,5];sigma[2,1]=cvalues[1,6];sigma[2,2]=cvalues[1,7]
  x=rmvnorm(n=nds/10, mean=mean, sigma=sigma)
  
  y1_susc[,j]=(1/(1+exp(-(x[,1]+x[,2]*log2(Susc_Symp_Rat$NATR_tot[j])))))
}

#Relative risk infectiousness
y1_inf=matrix(data=NA,nrow =nds/10, ncol = nrow(Susc_Symp_Rat) )
for (j in 1:nrow(Susc_Symp_Rat)) {
  mean=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma[1,1]=cvalues[2,4];sigma[1,2]=cvalues[2,5];sigma[2,1]=cvalues[2,6];sigma[2,2]=cvalues[2,7]
  x=rmvnorm(n=nds/10, mean=mean, sigma=sigma)
  
  y1_inf[,j]=(1/(1+exp(-(x[,1]+x[,2]*log2(Susc_Symp_Rat$NATR_tot[j])))))
}

##Transmission
y1_trans=y1_susc*y1_inf

#Relative risk symptomatic
y1_symp=matrix(data=NA,nrow =nds/10, ncol = nrow(Susc_Symp_Rat) )

#Relative Risk from Gardner and Kilpatrick, Viruses 2024
for (j in 1:nrow(Susc_Symp_Rat)) {
  mean=c(-1.607672,-0.5336916)
  sigma[1,1]=0.0008021563;sigma[1,2]=-4.403608e-05;sigma[2,1]=-4.403608e-05;sigma[2,2]=0.000965282
  x=rmvnorm(n=nds/10, mean=mean, sigma=sigma)
  
  y1_symp[,j]=(1/(1+exp(-(x[,1]+x[,2]*log2(Susc_Symp_Rat$NATR_tot[j])))))
}

##RR ratios 
y_rat=y1_trans/y1_symp
y_rat2=y1_susc/y1_symp


for (j in 1:nrow(Susc_Symp_Rat)) {
  Susc_Symp_Rat$ratio[j]=quantile(y_rat[,j],probs=0.5)
  Susc_Symp_Rat$lowerCI[j]=quantile(y_rat[,j],probs=0.025)
  Susc_Symp_Rat$upperCI[j]=quantile(y_rat[,j],probs=0.975)
}

for (j in 1:nrow(Susc_Symp_Rat)) {
  Susc_Symp_Rat$ratio2[j]=quantile(y_rat2[,j],probs=0.5)
  Susc_Symp_Rat$lowerCI2[j]=quantile(y_rat2[,j],probs=0.025)
  Susc_Symp_Rat$upperCI2[j]=quantile(y_rat2[,j],probs=0.975)
}

##Prep for plotting
trans_symp_rat=Susc_Symp_Rat[,c("NATR_tot","ratio","lowerCI","upperCI")]
trans_symp_rat$plot="B: Ratio of RR for transmission to RR for symptomatic disease"

susc_symp_rat=trans_symp_rat
susc_symp_rat[,c("NATR_tot","ratio","lowerCI","upperCI")]=Susc_Symp_Rat[,c("NATR_tot","ratio2","lowerCI2","upperCI2")]
susc_symp_rat$plot="A: Ratio of RR for susceptibility to RR for symptomatic disease"
Supp_Plot3=rbind(trans_symp_rat,susc_symp_rat)

#Supplemental Plot
ggplot(data=Supp_Plot3)+
  geom_ribbon(aes(x=NATR_tot,ymin=lowerCI, ymax=upperCI),fill="gray",alpha=0.5)+
  geom_line(aes(x=NATR_tot,y=ratio),color="black")+
  theme_few()+
  scale_x_continuous(expand = c(0,0),trans="log2",limits=c(0.125,4),breaks=c(0.0078125*2^c(4:9)),labels = label_number(accuracy = 0.001, scale = 1,
                                                                                                                       prefix = "", suffix = "",
                                                                                                                       big.mark = " ", decimal.mark = "."))+
  scale_y_continuous(expand = c(0, 0))+
  scale_shape_manual(values = c(18,19,17,6,15))+
  coord_cartesian(ylim = c(0, NA),xlim = c(NA,NA)) +
  xlab(expression("Neutralizing antibody titer ratio ("*NATR[tot]*")"))+
  ylab("Ratio of Relative Risk")+
  theme(strip.text = element_text(size=18),
        axis.title=element_text(size=18),#legend.position = "top",#legend.position = c(.1, .80),
        axis.text=element_text(size=18),legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        plot.margin = margin(1,1,1.5,1.2, "cm"))+
  labs(col="Variant", shape="Vaccine")+
  facet_wrap(vars(plot),ncol=1,scales="free")


#Waning figure - load data
m2=read.csv("nAbs_waning.csv") #Relative (0-1) antibody titers over time

#Fit w/ nls; model for relative nAbs
g6b=nls(log2(nAbs)~((c0+c7*Infection)*exp(c1*Day+c3*Day*Moderna)  -c0-c7*Infection),#best model
        start=list(c0=.8,c1=-0.02,c3=0.01,c7=0.01),
        data=m2,trace = F,control = list(warnOnly=T));summary(g6b)#Table S4

days=0:365
p1=data.frame(Day=rep(days,3),Infection=rep(c(0,1,0),each=length(days)),
              Moderna=rep(c(0,0,1),each=length(days)),
              Vaccine=rep(c("Pfizer","Infection","Moderna"),each=length(days)))

##SE and mean
#https://stackoverflow.com/questions/12425920/predict-cannot-display-the-standard-errors-of-the-predictions-with-se-fit-true
p1$varnt=deltavar((c0+c7*p1$Infection)*exp(c1*p1$Day+c3*p1$Day*p1$Moderna)  -c0-c7*p1$Infection,
                  meanval=coef(g6b),Sigma=vcov(g6b))
p1$se=p1$varnt^.5
p1$nAbs=predict(g6b,newdata=p1)
p1[,c("lower","upper")]=c(p1$nAbs-p1$se,p1$nAbs+p1$se)

#Making absolute neutralizing antibody predicted lines - rescaling of relative values
q1=p1;q1$response="absolute_nAbs";p1$response="relative_nAbs"
q1[q1$Vaccine=="Pfizer",c("nAbs","lower","upper")]=2^q1[q1$Vaccine=="Pfizer",c("nAbs","lower","upper")]*NATR_var_Delta*NATR_pfizer_recent_two
q1[q1$Vaccine=="Hybrid",c("nAbs","lower","upper")]=2^q1[q1$Vaccine=="Pfizer",c("nAbs","lower","upper")]*NATR_var_Delta*NATR_pfizer_recent_two*NATR_hybrid_ratio
q1[q1$Vaccine=="Infection",c("nAbs","lower","upper")]=2^q1[q1$Vaccine=="Infection",c("nAbs","lower","upper")]*NATR_var_Delta
q1[q1$Vaccine=="Moderna",c("nAbs","lower","upper")]=2^q1[q1$Vaccine=="Moderna",c("nAbs","lower","upper")]*NATR_var_Delta*NATR_moderna_recent_two

#Making absolute neutralizing antibody points - rescaling of relative values
n1=m2;
n1[n1$Vaccine=="Pfizer",c("nAbs")]=n1[n1$Vaccine=="Pfizer",c("nAbs")]*NATR_var_Delta*NATR_pfizer_recent_two
n1[n1$Vaccine=="Infection",c("nAbs")]=n1[n1$Vaccine=="Infection",c("nAbs")]*NATR_var_Delta
n1[n1$Vaccine=="Moderna",c("nAbs")]=n1[n1$Vaccine=="Moderna",c("nAbs")]*NATR_var_Delta*NATR_moderna_recent_two
m2$response="relative_nAbs";n1$response="absolute_nAbs"
n1=rbind(m2,n1)

#Making protection vs time dataframe
r1=p1;r1$response="VE_s";
names(r1)[2]="nAbs2" #renaming nAbs in r1 to nAbs2 b/c nAbs is column used for Y in plot
s1=p1;s1$response="VE_i";
names(s1)[2]="nAbs2" #renaming nAbs in s1 to nAbs2 b/c nAbs is column used for Y in plot

#Create dataframe for prediction outputs from model
y_susc_pfizer=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Pfizer",]) )
y_susc_moderna=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Moderna",]) )
y_susc_inf=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Infection",]) )

#Predictions for VE
for (i in 1:1) {
  mmPfi=r1[r1$Vaccine=="Pfizer",]
  mmMod=r1[r1$Vaccine=="Moderna",]
  mmInf=r1[r1$Vaccine=="Infection",]
  
  sigma=diag(NA,nrow=2,ncol=2)
  mean=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma[1,1]=cvalues$sigma11[cvalues$endpoint=="Susceptibility"];sigma[1,2]=cvalues$sigma12[cvalues$endpoint=="Susceptibility"];
  sigma[2,1]=cvalues$sigma21[cvalues$endpoint=="Susceptibility"];sigma[2,2]=cvalues$sigma22[cvalues$endpoint=="Susceptibility"]
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  
  sigma=diag(NA,nrow=2,ncol=2)
  mean=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  
  
  j=150
  for (j in 1:nrow(mmPfi)) {
    # for (v in 1:nds) {
    r1$nAbs[r1$Vaccine=="Pfizer"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^(mmPfi$nAbs[j])*NATR_pfizer_recent_two*NATR_var_Delta)))))
    r1$nAbs[r1$Vaccine=="Moderna"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^(mmMod$nAbs[j])*NATR_moderna_recent_two*NATR_var_Delta)))))
    r1$nAbs[r1$Vaccine=="Infection"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^(mmInf$nAbs[j])*NATR_var_Delta)))))
    
    #r1$nAbs[j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^(r1$nAbs[j])*NATR_pfizer_recent_two*NATR_var_Delta)))))
    # 
    y_susc_pfizer[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean= mmPfi$nAbs[j],sd=mmPfi$se[j])))*NATRpfi/NATRdelta)))))
    
    y_susc_moderna[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean = mmMod$nAbs[j],sd=mmMod$se[j])))*NATRmod/NATRdelta)))))
    
    y_susc_inf[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean = mmInf$nAbs[j],sd=mmInf$se[j])))*NATRcon/NATRdelta)))))
    
    #y_susc[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2(2^((rlnorm(nds,meanlog = r1$nAbs[j],sd=r1$se[j])))*NATRpfi/NATRdelta)))))
    # }
  }
  # y_susc_pfizer[is.na(y_susc_pfizer)] <- 0
  # y_susc_moderna[is.na(y_susc_moderna)] <- 0
  # y_susc_inf[is.na(y_susc_inf)] <- 0
  
  for (k in 1:nrow(mmPfi)) {
    r1$lower[r1$Vaccine=="Pfizer"][k+nrow(mmPfi)*(i-1)]=quantile(y_susc_pfizer[,k],probs=0.16)
    r1$upper[r1$Vaccine=="Pfizer"][k+nrow(mmPfi)*(i-1)]=quantile(y_susc_pfizer[,k],probs=0.84)
    r1$lower[r1$Vaccine=="Moderna"][k+nrow(mmMod)*(i-1)]=quantile(y_susc_moderna[,k],probs=0.16)
    r1$upper[r1$Vaccine=="Moderna"][k+nrow(mmMod)*(i-1)]=quantile(y_susc_moderna[,k],probs=0.84)
    r1$lower[r1$Vaccine=="Infection"][k+nrow(mmInf)*(i-1)]=quantile(y_susc_inf[,k],probs=0.16)
    r1$upper[r1$Vaccine=="Infection"][k+nrow(mmInf)*(i-1)]=quantile(y_susc_inf[,k],probs=0.84)
  }
}

y_infect_pfizer=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Pfizer",]) )
y_infect_moderna=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Moderna",]) )
y_infect_inf=matrix(data=NA,nrow =nds, ncol = nrow(r1[r1$Vaccine=="Infection",]) )

for (i in 1:1) {
  mmPfi=s1[s1$Vaccine=="Pfizer",]
  mmMod=s1[s1$Vaccine=="Moderna",]
  mmInf=s1[s1$Vaccine=="Infection",]
  
  sigma=diag(NA,nrow=2,ncol=2)
  mean=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma[1,1]=cvalues$sigma11[cvalues$endpoint=="Infectiousness"];sigma[1,2]=cvalues$sigma12[cvalues$endpoint=="Infectiousness"];
  sigma[2,1]=cvalues$sigma21[cvalues$endpoint=="Infectiousness"];sigma[2,2]=cvalues$sigma22[cvalues$endpoint=="Infectiousness"]
  x=rmvnorm(n=nds, mean=mean, sigma=sigma)
  for (j in 1:nrow(mmPfi)) {
    s1$nAbs[r1$Vaccine=="Pfizer"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^(mmPfi$nAbs[j])*NATR_pfizer_recent_two*NATR_var_Delta)))))
    s1$nAbs[r1$Vaccine=="Moderna"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^(mmMod$nAbs[j])*NATR_moderna_recent_two*NATR_var_Delta)))))
    s1$nAbs[r1$Vaccine=="Infection"][j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^(mmInf$nAbs[j])*NATR_var_Delta)))))
    
    #r1$nAbs[j]=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^(r1$nAbs[j])*NATR_pfizer_recent_two*NATR_var_Delta)))))
    # 
    y_infect_pfizer[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean= mmPfi$nAbs[j],sd=mmPfi$se[j])))*NATRpfi/NATRdelta)))))
    
    y_infect_moderna[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean = mmMod$nAbs[j],sd=mmMod$se[j])))*NATRmod/NATRdelta)))))
    
    y_infect_inf[,j]=1-(1/(1+exp(-(x[,1]+x[,2]*log2((2^(rnorm(nds,mean = mmInf$nAbs[j],sd=mmInf$se[j])))*NATRcon/NATRdelta)))))
  }
  # y_infect_pfizer[is.na(y_infect_pfizer)] <- 0
  # y_infect_moderna[is.na(y_infect_moderna)] <- 0
  # y_infect_inf[is.na(y_infect_inf)] <- 0
  for (k in 1:nrow(mmPfi)) {
    s1$lower[s1$Vaccine=="Pfizer"][k+nrow(mmPfi)*(i-1)]=quantile(y_infect_pfizer[,k],probs=0.16)
    s1$upper[s1$Vaccine=="Pfizer"][k+nrow(mmPfi)*(i-1)]=quantile(y_infect_pfizer[,k],probs=0.84)
    s1$lower[s1$Vaccine=="Moderna"][k+nrow(mmMod)*(i-1)]=quantile(y_infect_moderna[,k],probs=0.16)
    s1$upper[s1$Vaccine=="Moderna"][k+nrow(mmMod)*(i-1)]=quantile(y_infect_moderna[,k],probs=0.84)
    s1$lower[s1$Vaccine=="Infection"][k+nrow(mmInf)*(i-1)]=quantile(y_infect_inf[,k],probs=0.16)
    s1$upper[s1$Vaccine=="Infection"][k+nrow(mmInf)*(i-1)]=quantile(y_infect_inf[,k],probs=0.84)
  }
}

#eliminating working columns to bind to other dataframes
r1=r1[,c("Day","nAbs","lower","upper","Vaccine","response")]
s1=s1[,c("Day","nAbs","lower","upper","Vaccine","response")]
p1=p1[,c("Day","nAbs","lower","upper","Vaccine","response")]
q1=q1[,c("Day","nAbs","lower","upper","Vaccine","response")]

#unlog p1
p1[,c("nAbs","lower","upper")]=2^p1[,c("nAbs","lower","upper")]


#combining absolute and relative antibodies and protection into 1 dataframe for plotting
p1=rbind(p1,q1,r1,s1) 

#Panel labels
tlab3=data.frame(response=c("absolute_nAbs","relative_nAbs","VE_s","VE_i"),
                 text1=c("A: Neutralizing antibody titer relative to convalescent sera with WT virus","B: Neutralizing antibody titer relative to peak from panel A",
                         "C: VE for susceptibility", "D: VE for infectiousness"),Day=c(50),nAbs=c(1.93,1.005,.90,0.78))

pred_delta_boost=data.frame(Variant="Delta",vaccine=rep(c("Pfizer","Moderna"),2),
                            endpoint=rep(c(rep(c("Susceptibility","Infectiousness" ),each=2)),1),
                            nAbs=1,NATR_tot=NA,x_lower=NA,x_upper=NA, mean=NA,
                            lowerPI=NA,upperPI=NA,lowerCI=NA,upperCI=NA)


#Create dataframe for prediction outputs from model
y1PfiS=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModS=matrix(data=NA,nrow =nds, ncol = 1 )
y1PfiI=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModI=matrix(data=NA,nrow =nds, ncol = 1 )



# #baseline nabs for Pfizer/immune status
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Pfizer"&pred_delta$status=="waned"]=NATR_pfizer_waned_two
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Pfizer"&pred_delta$status=="full"]=NATR_pfizer_recent_two
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Pfizer"&pred_delta$status=="boosted"]=NATR_pfizer_recent_three
# 
# #baseline nabs for Moderna/immune status
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Moderna"&pred_delta$status=="waned"]=NATR_moderna_waned_two
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Moderna"&pred_delta$status=="full"]=NATR_moderna_recent_two
# pred_delta$NATR_vac_istatus[pred_delta$vaccine=="Moderna"&pred_delta$status=="boosted"]=NATR_moderna_recent_three
# 
# #baseline nabs for D614G prediction including upper/lower CIs
# pred_delta$NATR_tot[pred_delta$Prediction=="D614G"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="D614G"]/FC3$fold_red[FC3$Variant=="D614G"]
# pred_delta$NATR_tot_lower[pred_delta$Prediction=="D614G"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="D614G"]/FC3$upper[FC3$Variant=="D614G"]
# pred_delta$NATR_tot_upper[pred_delta$Prediction=="D614G"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="D614G"]/FC3$lower[FC3$Variant=="D614G"]
# 
# #baseline nabs for Delta prediction including upper/lower CIs
# pred_delta$NATR_tot[pred_delta$Prediction=="Delta"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="Delta"]/FC3$fold_red[FC3$Variant=="Delta"]
# pred_delta$NATR_tot_lower[pred_delta$Prediction=="Delta"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="Delta"]/FC3$upper[FC3$Variant=="Delta"]
# pred_delta$NATR_tot_upper[pred_delta$Prediction=="Delta"]=pred_delta$NATR_vac_istatus[pred_delta$Prediction=="Delta"]/FC3$lower[FC3$Variant=="Delta"]


#Predictions for pfizer/moderna boosted
for (i in 1:1) {
  mm1=pred_delta_boost
  mean_susc=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma_susc=diag(NA,nrow=2,ncol=2)
  sigma_susc[1,1]=cvalues[1,4];sigma_susc[1,2]=cvalues[1,5];sigma_susc[2,1]=cvalues[1,6];sigma_susc[2,2]=cvalues[1,7]
  x_susc=rmvnorm(n=nds, mean=mean_susc, sigma=sigma_susc)
  
  mean_infect=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma_infect=diag(NA,nrow=2,ncol=2)
  sigma_infect[1,1]=cvalues[1,4];sigma_infect[1,2]=cvalues[1,5];sigma_infect[2,1]=cvalues[1,6];sigma_infect[2,2]=cvalues[1,7]
  x_infect=rmvnorm(n=nds, mean=mean_infect, sigma=sigma_infect)
  
  y1PfiS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility",]$nAbs*NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta))))))
  y1PfiI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness",]$nAbs*NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta))))))
  y1ModS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility",]$nAbs*NATRmod*NATR_moderna_recent_three/(NATR_moderna_recent_two*NATRdelta))))))
  y1ModI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness",]$nAbs*NATRmod*NATR_moderna_recent_three/(NATR_moderna_recent_two*NATRdelta))))))
  
  
  # for(l in 1:nrow(y1)){
  #   for(n in 1:ncol(y1)){
  #     VE_l1[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y1[l,n]))/rbinom(1,10000,prob=base_prev)
  #   }
  # }
  pred_delta_boost$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))
  pred_delta_boost$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.025)
  pred_delta_boost$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.975)
  
  pred_delta_boost$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))
  pred_delta_boost$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.025)
  pred_delta_boost$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.975)
  
  pred_delta_boost$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))
  pred_delta_boost$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.025)
  pred_delta_boost$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.975)
  
  pred_delta_boost$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))
  pred_delta_boost$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.025)
  pred_delta_boost$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.975)
  
}

pred_delta_two=data.frame(Variant="Delta",vaccine=rep(c("Pfizer","Moderna"),2),
                          endpoint=rep(c(rep(c("Susceptibility","Infectiousness" ),each=2)),1),
                          nAbs=1,NATR_tot=NA,x_lower=NA,x_upper=NA, mean=NA,
                          lowerPI=NA,upperPI=NA,lowerCI=NA,upperCI=NA)


#Create dataframe for prediction outputs from model
y1PfiS=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModS=matrix(data=NA,nrow =nds, ncol = 1 )
y1PfiI=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModI=matrix(data=NA,nrow =nds, ncol = 1 )

#Predictions for pfizer/moderna recent two-dose
for (i in 1:1) {
  mm1=pred_delta_two
  mean_susc=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma_susc=diag(NA,nrow=2,ncol=2)
  sigma_susc[1,1]=cvalues[1,4];sigma_susc[1,2]=cvalues[1,5];sigma_susc[2,1]=cvalues[1,6];sigma_susc[2,2]=cvalues[1,7]
  x_susc=rmvnorm(n=nds, mean=mean_susc, sigma=sigma_susc)
  
  mean_infect=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma_infect=diag(NA,nrow=2,ncol=2)
  sigma_infect[1,1]=cvalues[1,4];sigma_infect[1,2]=cvalues[1,5];sigma_infect[2,1]=cvalues[1,6];sigma_infect[2,2]=cvalues[1,7]
  x_infect=rmvnorm(n=nds, mean=mean_infect, sigma=sigma_infect)
  
  y1PfiS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility",]$nAbs*NATRpfi/(NATRdelta))))))
  y1PfiI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness",]$nAbs*NATRpfi/(NATRdelta))))))
  y1ModS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility",]$nAbs*NATRmod/(NATRdelta))))))
  y1ModI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness",]$nAbs*NATRmod/(NATRdelta))))))
  
  
  # for(l in 1:nrow(y1)){
  #   for(n in 1:ncol(y1)){
  #     VE_l1[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y1[l,n]))/rbinom(1,10000,prob=base_prev)
  #   }
  # }
  pred_delta_two$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_two*NATR_var_Delta))))))
  pred_delta_two$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.025)
  pred_delta_two$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.975)
  
  pred_delta_two$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_two*NATR_var_Delta))))))
  pred_delta_two$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.025)
  pred_delta_two$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.975)
  
  pred_delta_two$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_moderna_recent_two*NATR_var_Delta))))))
  pred_delta_two$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.025)
  pred_delta_two$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.975)
  
  pred_delta_two$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_moderna_recent_two*NATR_var_Delta))))))
  pred_delta_two$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.025)
  pred_delta_two$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.975)
}

#Predictions for moderna/pfizer waned
pred_delta_waned=data.frame(Variant="Delta",vaccine=rep(c("Pfizer","Moderna"),2),
                            endpoint=rep(c(rep(c("Susceptibility","Infectiousness" ),each=2)),1),
                            nAbs=1,NATR_tot=NA,x_lower=NA,x_upper=NA, mean=NA,
                            lowerPI=NA,upperPI=NA,lowerCI=NA,upperCI=NA)

pred_delta_waned$nAbs_pfizer=predict(g6b,newdata=data.frame(Day=240,Infection=0,Moderna=0))
pred_delta_waned$nAbs_pfizer_var=deltavar((c0)*exp(c1*240)  -c0,
                                          meanval=coef(g6b),Sigma=vcov(g6b))
pred_delta_waned$nAbs_pfizer_se=pred_delta_waned$nAbs_pfizer_var^0.5

pred_delta_waned$nAbs_moderna=predict(g6b,newdata=data.frame(Day=240,Infection=0,Moderna=1))
pred_delta_waned$nAbs_moderna_var=deltavar((c0)*exp(c1*240+c3*240)  -c0,
                                           meanval=coef(g6b),Sigma=vcov(g6b))
pred_delta_waned$nAbs_moderna_se=pred_delta_waned$nAbs_moderna_var^0.5

#Create dataframe for prediction outputs from model
y1PfiS=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModS=matrix(data=NA,nrow =nds, ncol = 1 )
y1PfiI=matrix(data=NA,nrow =nds, ncol = 1 )
y1ModI=matrix(data=NA,nrow =nds, ncol = 1 )

#Predictions
for (i in 1:1) {
  mm1=pred_delta_waned
  mean_susc=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
  sigma_susc=diag(NA,nrow=2,ncol=2)
  sigma_susc[1,1]=cvalues[1,4];sigma_susc[1,2]=cvalues[1,5];sigma_susc[2,1]=cvalues[1,6];sigma_susc[2,2]=cvalues[1,7]
  x_susc=rmvnorm(n=nds, mean=mean_susc, sigma=sigma_susc)
  
  mean_infect=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
  sigma_infect=diag(NA,nrow=2,ncol=2)
  sigma_infect[1,1]=cvalues[1,4];sigma_infect[1,2]=cvalues[1,5];sigma_infect[2,1]=cvalues[1,6];sigma_infect[2,2]=cvalues[1,7]
  x_infect=rmvnorm(n=nds, mean=mean_infect, sigma=sigma_infect)
  
  y1PfiS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility",]$nAbs_pfizer,sd=mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility",]$nAbs_pfizer_se))*NATRpfi/(NATRdelta))))))
  y1PfiI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness",]$nAbs_pfizer,sd=mm1[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness",]$nAbs_pfizer_se))*NATRpfi/(NATRdelta))))))
  y1ModS[,i]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility",]$nAbs_moderna,sd=mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility",]$nAbs_moderna_se))*NATRmod/(NATRdelta))))))
  y1ModI[,i]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness",]$nAbs_moderna,sd=mm1[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness",]$nAbs_moderna_se))*NATRmod/(NATRdelta))))))
  
  
  # for(l in 1:nrow(y1)){
  #   for(n in 1:ncol(y1)){
  #     VE_l1[l,n]=1-rbinom(1,10000,prob=base_prev*(1-y1[l,n]))/rbinom(1,10000,prob=base_prev)
  #   }
  # }
  pred_delta_waned$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_waned_two*NATR_var_Delta))))))
  pred_delta_waned$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.025)
  pred_delta_waned$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Susceptibility"]=quantile(y1PfiS,probs=0.975)
  
  pred_delta_waned$mean[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_waned_two*NATR_var_Delta))))))
  pred_delta_waned$lowerCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.025)
  pred_delta_waned$upperCI[mm1$vaccine=="Pfizer"&mm1$endpoint=="Infectiousness"]=quantile(y1PfiI,probs=0.975)
  
  pred_delta_waned$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_moderna_waned_two*NATR_var_Delta))))))
  pred_delta_waned$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.025)
  pred_delta_waned$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Susceptibility"]=quantile(y1ModS,probs=0.975)
  
  pred_delta_waned$mean[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_moderna_waned_two*NATR_var_Delta))))))
  pred_delta_waned$lowerCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.025)
  pred_delta_waned$upperCI[mm1$vaccine=="Moderna"&mm1$endpoint=="Infectiousness"]=quantile(y1ModI,probs=0.975)
}

#make point for boosted on Figure 2
br1=pred_delta_boost[pred_delta_boost$endpoint=="Susceptibility",]
bp=data.frame(Day=c(365,364),nAbs=br1$mean,lower=br1$lowerCI,upper=br1$upperCI,response="VE_s")
br2=pred_delta_boost[pred_delta_boost$endpoint=="Infectiousness",]
bp2=data.frame(Day=c(365,364),nAbs=br2$mean,lower=br2$lowerCI,upper=br2$upperCI,response="VE_i")

#Reorder for plotting
p1$response=factor(p1$response,levels = c("absolute_nAbs","relative_nAbs","VE_s","VE_i"))
n1$response=factor(n1$response,levels = c("absolute_nAbs","relative_nAbs","VE_s","VE_i"))
tlab3$response=factor(tlab3$response,levels = c("absolute_nAbs","relative_nAbs","VE_s","VE_i"))
bp$response=factor(bp$response,levels = c("absolute_nAbs","relative_nAbs","VE_s","VE_i"))
bp2$response=factor(bp2$response,levels = c("absolute_nAbs","relative_nAbs","VE_s","VE_i"))

n1$Vaccine[n1$Vaccine=="Pfizer"]="BNT162b2"
n1$Vaccine[n1$Vaccine=="Moderna"]="mRNA-1273"
p1$Vaccine[p1$Vaccine=="Pfizer"]="BNT162b2"
p1$Vaccine[p1$Vaccine=="Moderna"]="mRNA-1273"

mylabeller=as_labeller(c(absolute_nAbs="Antibody~titer",
                         relative_nAbs="Relative~antibody~titer",
                         VE_s="VE[S]",
                         VE_i="VE[I]"),default = label_parsed)

#Figure 2
ggplot()+
  geom_point(data=n1,aes(x=Day,y=nAbs,col=Vaccine),size=3)+
  geom_ribbon(data=p1[p1$Vaccine!="Hybrid",],aes(x=Day,y=nAbs,ymin=lower, ymax=upper,fill=Vaccine),colour=NA,alpha=0.1)+
  geom_line(data=p1[p1$Vaccine!="Hybrid",],aes(x=Day,y=nAbs,col=Vaccine))+theme_few()+
  #ylab(expression("Neutralizing antibody titres, "*VE[S]*", or "*VE[I]))+
  xlab("Days since infection/vaccination")+#ylab("Neutralizing antibody titres or VE against susceptibility")+
  labs(col="Immunity", fill="Immunity")+
  theme(axis.title=element_text(size=18),#legend.position = c(.1, .80),
        axis.text=element_text(size=18),legend.text=element_text(size=15),
        legend.title=element_text(size=20),strip.text = element_text(size=17))+
  theme(strip.text.x = element_blank())+
  scale_color_manual(values=c("#0066FF","#00BA38","#C77CFF"))+
  scale_fill_manual(values=c("#0066FF","#00BA38","#C77CFF"))+
  geom_text(data=tlab3,aes(x=Day,y=nAbs,label=text1),size=6,hjust=0)+
  geom_point(data=bp2,aes(x=Day,y=nAbs),size=3,col=c("#0066FF", "#C77CFF"))+
  geom_text(data=bp2,aes(x=Day,y=nAbs+0.02*c(1,-1),label=c("BNT162b2 Boosted","mRNA-1273 Boosted")),size=5,hjust=1.1)+
  geom_errorbar(data=bp2,aes(x=Day,ymin=lower, ymax=upper), width=0,col=c("#0066FF", "#C77CFF"))+
  geom_point(data=bp,aes(x=Day,y=nAbs),size=3,col=c("#0066FF", "#C77CFF"))+
  geom_text(data=bp,aes(x=Day,y=nAbs+0.02*c(1,-1),label=c("BNT162b2 Boosted","mRNA-1273 Boosted")),size=5,hjust=1.1)+
  geom_errorbar(data=bp,aes(x=Day,ymin=lower, ymax=upper), width=0,col=c("#0066FF", "#C77CFF"))+
  theme(axis.title.y = element_blank(),legend.position="top",strip.background=element_blank(),strip.placement="outside")+  
  facet_wrap(.~response,nrow=4, scales="free_y",strip.position = "left",
             labeller = mylabeller)

#Load file
vacc_irr0=read.csv("irr_unvax_vax2.csv") #Infection ratio unvax v. vax
vacc_irr=vacc_irr0[vacc_irr0$time<275,]

#GAM
g0=gam(IRR_cases ~s(time),
       data = vacc_irr,method = "REML");summary(g0)

#Plot data with model
ggplot(data = vacc_irr,
       aes(x=time, y=IRR_cases)) + geom_point()+
  geom_smooth(method="gam",formula = y ~ s(x))+
  xlab("")+
  ylab("Ratio of COVID-19 cases in unvaccinated v. vaccinated")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position="bottom")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=25),
        legend.text=element_text(size=25),legend.title=element_blank())+
  scale_x_continuous(expand = c(0,0), breaks=c(15+30*(0:8)),
                     labels = c("Apr 21","May 21","Jun 21",
                                "Jul 21","Aug 21","Sep 21","Oct 21",
                                "Nov 21","Dec 21"))+
  scale_y_continuous(
    expand = c(0, 0), limits=c(0,14))

#GAM
g1=gam(IRR_deaths ~s(time),
       data = vacc_irr,method = "REML");summary(g0)

#Plot data with model
ggplot(data = vacc_irr,
       aes(x=time, y=IRR_deaths)) + geom_point()+
  geom_smooth(method="gam",formula = y ~ s(x))+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position="bottom")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=25),
        legend.text=element_text(size=25),legend.title=element_blank())+
  scale_x_continuous(expand = c(0,0), n.breaks=10)+
  scale_y_continuous(
    expand = c(0, 0), limits=c(0,23))

#Start date (first day with data) and end date (last day with data), first day vaccine given, and population
start_date = '2020-01-23'; end_date = '2022-01-06';first_vax="2020-12-14";population=332858873
new_end='2021-12-13'
#Fraction of mRNA vaccines given that are pfizer
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

#Read in data set with USA daily cases, new fully vaccinated, and deaths w/ row 1 being most recent data
case_vax_data=read.csv("Cases_Vax_Death3.csv") #load case and vax data

#Add columns with date, number of infections, and proportion of 2nd doses given on a specific day
case_vax_data$date=seq(as.Date(end_date), as.Date(start_date), by = "-1 days")

#Rolling mean and remove NAs
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
infec_deaths_US=deconvolve_incidence(incidence_data = case_vax_data$daily_deaths_sm,
                                     delay=list(inf_shed,shed_sympt,sympt_hosp,hosp_death))
case_vax_data$inf=
  (1*infec_deaths_US$values[c(rep(NA,-infec_deaths_US$index_offset),
                              (1):(infec_deaths_US$index_offset+length(infec_deaths_US$values)))])/0.003

#remove rows with NAs in infections
case_vax_data=case_vax_data[complete.cases(case_vax_data$inf),]

#cumulative infections
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

#Store dataframe for Figure 3BC creation later
case_vax_data_hold=case_vax_data

#Function that takes data and dataframe with population infection/vaccination data and gives VE values for that time
VE_at_time=function(cutoff_date,daily_cases_vaxes){
  #Cutoff
  mm2=daily_cases_vaxes[daily_cases_vaxes$date<=cutoff_date,]
  mm2$prop_hybrid=(mm2$daily_hybrid/sum(mm2$daily_hybrid))
  mm2$prop_unvax_pos=(mm2$daily_unvax_pos/sum(mm2$daily_unvax_pos))
  mm2$prop_vax_neg_adj=(mm2$daily_vax_neg_adj/sum(mm2$daily_vax_neg_adj))
  
  #Create dataframe to calculate VEs
  #Pull relevant columns from mm2
  VE_data=data.frame(time=c(rep(0,15),1:nrow(mm2)),prop_inf_unvax=NA,prop_hybrid=NA,prop_pfizer=NA,prop_moderna=NA)
  VE_data$prop_inf_unvax[8:(nrow(mm2)+8-24)]=mm2$prop_unvax_pos[24:nrow(mm2)]
  VE_data$prop_hybrid[15:(nrow(mm2)+15-1)]=mm2$prop_hybrid
  VE_data$prop_pfizer[8:(nrow(mm2)+8-1)]=mm2$prop_vax_neg_adj
  VE_data$prop_moderna[1:(nrow(mm2)+1-1)]=mm2$prop_vax_neg_adj
  VE_data[is.na(VE_data)] = 0
  
  #Add columns with predicted nAbs for each class
  VE_data$nAbs_inf_unvax=predict(g6b,newdata=data.frame(Day=VE_data$time,Infection=1,Moderna=0))
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
  VE_data$nAbs_hybrid=VE_data$nAbs_pfizer
  VE_data$nAbs_hybrid_se=VE_data$nAbs_pfizer_se
  
  #Susceptibility
  VE_data$V_s_pfi=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^VE_data$nAbs_pfizer*NATR_pfizer_recent_two*NATR_var_Delta)))))
  
  VE_data$V_s_mod=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^VE_data$nAbs_moderna*NATR_moderna_recent_two*NATR_var_Delta)))))
  
  VE_data$V_s_hybrid=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^VE_data$nAbs_pfizer*NATR_pfizer_recent_two*NATR_hybrid_ratio*NATR_var_Delta)))))
  
  VE_data$V_s_inf=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(2^VE_data$nAbs_inf_unvax*NATR_var_Delta)))))
  
  VE_data$V_s_boost=mRNA_pfizer_ratio*(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))+
    mRNA_moderna_ratio*(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))
  
  VE_data$V_s_boost_pfi=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))
  
  
  #Infectiousness
  VE_data$V_i_pfi=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^VE_data$nAbs_pfizer*NATR_pfizer_recent_two*NATR_var_Delta)))))
  
  VE_data$V_i_mod=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^VE_data$nAbs_moderna*NATR_moderna_recent_two*NATR_var_Delta)))))
  
  VE_data$V_i_hybrid=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^VE_data$nAbs_pfizer*NATR_pfizer_recent_two*NATR_hybrid_ratio*NATR_var_Delta)))))
  
  VE_data$V_i_inf=1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(2^VE_data$nAbs_inf_unvax*NATR_var_Delta)))))
  
  VE_data$V_i_boost=mRNA_pfizer_ratio*(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))+
    mRNA_moderna_ratio*(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_moderna_recent_three*NATR_var_Delta))))))
  
  VE_data$V_i_boost_pfi=(1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_three*NATR_var_Delta))))))
  
  #Create dataframe for prediction outputs from model
  y_susc_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_moderna=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_hybrid=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_inf=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_boost=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_susc_boost_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  
  y_infect_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_moderna=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_hybrid=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_inf=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_boost=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  y_infect_boost_pfizer=matrix(data=NA,nrow =nds2, ncol = nrow(VE_data) )
  
  #Create covariance matrices for calculations
  sigma_inf=diag(NA,nrow=2,ncol=2)
  sigma_trans=diag(NA,nrow=2,ncol=2)
  
  
  #Susceptibilitydata.frame()#Susceptibility and Infectiousness by day since most recent vaccination/infection by immune status for US
  for (i in 1:1) {
    mm1=VE_data
    mean_inf=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
    sigma_inf[1,1]=cvalues$sigma11[cvalues$endpoint=="Susceptibility"];sigma_inf[1,2]=cvalues$sigma12[cvalues$endpoint=="Susceptibility"];
    sigma_inf[2,1]=cvalues$sigma21[cvalues$endpoint=="Susceptibility"];sigma_inf[2,2]=cvalues$sigma22[cvalues$endpoint=="Susceptibility"]
    x_susc=rmvnorm(n=nds2, mean=mean_inf, sigma=sigma_inf)
    
    mean_trans=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
    sigma_trans[1,1]=cvalues$sigma11[cvalues$endpoint=="Infectiousness"];sigma_trans[1,2]=cvalues$sigma12[cvalues$endpoint=="Infectiousness"];
    sigma_trans[2,1]=cvalues$sigma21[cvalues$endpoint=="Infectiousness"];sigma_trans[2,2]=cvalues$sigma22[cvalues$endpoint=="Infectiousness"]
    x_infect=rmvnorm(n=nds2, mean=mean_trans, sigma=sigma_trans)
    
    for (j in 1:nrow(mm1)) {
      # for (v in 1:nds2) {
      #Susceptibility
      y_susc_pfizer[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi/NATRdelta)))))
      
      y_susc_moderna[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_moderna[j],sd=mm1$nAbs_moderna_se[j]))*NATRmod/NATRdelta)))))
      
      y_susc_hybrid[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi*NATR_hybrid_ratio/NATRdelta)))))
      
      y_susc_inf[,j]=1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_inf_unvax[j],sd=mm1$nAbs_inf_unvax_se[j]))*NATRcon/NATRdelta)))))
      
      y_susc_boost[,j]=mRNA_pfizer_ratio*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta)))))))+
        mRNA_moderna_ratio*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATR_moderna_recent_three*NATRmod/(NATR_moderna_recent_two*NATRdelta)))))))
      
      y_susc_boost_pfizer[,j]=(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta)))))))
      
      
      #Infectiousness
      y_infect_pfizer[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi/NATRdelta)))))
      
      y_infect_moderna[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_moderna[j],sd=mm1$nAbs_moderna_se[j]))*NATRmod/NATRdelta)))))
      
      y_infect_hybrid[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_pfizer[j],sd=mm1$nAbs_pfizer_se[j]))*NATRpfi*NATR_hybrid_ratio/NATRdelta)))))
      
      y_infect_inf[,j]=1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(2^(rnorm(nds2,mean= mm1$nAbs_inf_unvax[j],sd=mm1$nAbs_inf_unvax_se[j]))*NATRcon/NATRdelta)))))
      
      y_infect_boost[,j]=mRNA_pfizer_ratio*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta)))))))+
        mRNA_moderna_ratio*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATR_moderna_recent_three*NATRmod/(NATR_moderna_recent_two*NATRdelta)))))))
      
      y_infect_boost_pfizer[,j]=(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATR_pfizer_recent_three*NATRpfi/(NATR_pfizer_recent_two*NATRdelta)))))))    
      
      
      # }
    }
    
  }
  
  #Susceptibility VE and data frame to calculate CIs for each immune status based on US population waning
  CI_susc_pfizer=as.data.frame((y_susc_pfizer%*%VE_data$prop_pfizer))
  CI_susc_moderna=as.data.frame((y_susc_moderna%*%VE_data$prop_moderna))
  CI_susc_mRNA=as.data.frame((mRNA_pfizer_ratio*(y_susc_pfizer%*%VE_data$prop_pfizer)+mRNA_moderna_ratio*(y_susc_moderna%*%VE_data$prop_moderna)))
  CI_susc_hybrid=as.data.frame((y_susc_hybrid%*%VE_data$prop_hybrid))
  CI_susc_inf=as.data.frame((y_susc_inf%*%VE_data$prop_inf_unvax))
  CI_susc_boost_mRNA=as.data.frame(y_susc_boost[,1])
  CI_susc_boost_pfizer=as.data.frame(y_susc_boost_pfizer[,1])
  
  VE_susc_pfizer=y_susc_pfizer%*%VE_data$prop_pfizer
  VE_susc_moderna=y_susc_moderna%*%VE_data$prop_moderna
  VE_susc_mRNA=mRNA_pfizer_ratio*(y_susc_pfizer%*%VE_data$prop_pfizer)+mRNA_moderna_ratio*(y_susc_moderna%*%VE_data$prop_moderna)
  VE_susc_hybrid=y_susc_hybrid%*%VE_data$prop_hybrid
  VE_susc_inf=y_susc_inf%*%VE_data$prop_inf_unvax
  VE_susc_boost_mRNA=y_susc_boost[,1]
  VE_susc_boost_pfizer=y_susc_boost_pfizer[,1]
  
  #Infectiousness VE and data frame to calculate CIs for each immune status based on US population waning
  CI_infect_pfizer=as.data.frame((y_infect_pfizer%*%VE_data$prop_pfizer))
  CI_infect_moderna=as.data.frame((y_infect_moderna%*%VE_data$prop_moderna))
  CI_infect_mRNA=as.data.frame((mRNA_pfizer_ratio*(y_infect_pfizer%*%VE_data$prop_pfizer)+mRNA_moderna_ratio*(y_infect_moderna%*%VE_data$prop_moderna)))
  CI_infect_hybrid=as.data.frame((y_infect_hybrid%*%VE_data$prop_hybrid))
  CI_infect_inf=as.data.frame((y_infect_inf%*%VE_data$prop_inf_unvax))
  CI_infect_boost_mRNA=as.data.frame(y_infect_boost[,1])
  CI_infect_boost_pfizer=as.data.frame(y_infect_boost_pfizer[,1])
  
  
  VE_infect_pfizer=y_infect_pfizer%*%VE_data$prop_pfizer
  VE_infect_moderna=y_infect_moderna%*%VE_data$prop_moderna
  VE_infect_mRNA=mRNA_pfizer_ratio*(y_infect_pfizer%*%VE_data$prop_pfizer)+mRNA_moderna_ratio*(y_infect_moderna%*%VE_data$prop_moderna)
  VE_infect_hybrid=y_infect_hybrid%*%VE_data$prop_hybrid
  VE_infect_inf=y_infect_inf%*%VE_data$prop_inf_unvax
  VE_infect_boost_mRNA=y_infect_boost[,1]
  VE_infect_boost_pfizer=y_infect_boost_pfizer[,1]
  
  #Deterministic
  V_s_pfi=VE_data$V_s_pfi%*%VE_data$prop_pfizer
  V_s_mod=VE_data$V_s_mod%*%VE_data$prop_moderna
  V_s_hybrid=VE_data$V_s_hybrid%*%VE_data$prop_hybrid
  V_s_inf=VE_data$V_s_inf%*%VE_data$prop_inf_unvax
  V_s_boost=VE_data$V_s_boost[1]
  V_s_boost_pfi=VE_data$V_s_boost_pfi[1]
  
  V_i_pfi=VE_data$V_i_pfi%*%VE_data$prop_pfizer
  V_i_mod=VE_data$V_i_mod%*%VE_data$prop_moderna
  V_i_hybrid=VE_data$V_i_hybrid%*%VE_data$prop_hybrid
  V_i_inf=VE_data$V_i_inf%*%VE_data$prop_inf_unvax
  V_i_boost=VE_data$V_i_boost[1]
  V_i_boost_pfi=VE_data$V_i_boost_pfi[1]
  
  return(list(CI_susc_pfizer,CI_susc_moderna,CI_susc_mRNA,CI_susc_hybrid,
              CI_susc_inf,CI_susc_boost_mRNA,CI_susc_boost_pfizer,
              VE_susc_pfizer,VE_susc_moderna,VE_susc_mRNA,VE_susc_hybrid,
              VE_susc_inf,VE_susc_boost_mRNA,VE_susc_boost_pfizer,
              CI_infect_pfizer,CI_infect_moderna,CI_infect_mRNA,CI_infect_hybrid,
              CI_infect_inf,CI_infect_boost_mRNA,CI_infect_boost_pfizer,
              VE_infect_pfizer,VE_infect_moderna,VE_infect_mRNA,VE_infect_hybrid,
              VE_infect_inf,VE_infect_boost_mRNA,VE_infect_boost_pfizer,mm2,
              V_s_pfi,V_s_mod,V_s_hybrid,V_s_inf,V_s_boost,V_s_boost_pfi,
              V_i_pfi,V_i_mod,V_i_hybrid,V_i_inf,V_i_boost,V_i_boost_pfi))
}

#United States
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data$date<=as.Date("2021-09-01"),]))
US_S=VE_at_time(as.Date("2021-09-01"),case_vax_data)

#Store VE and CI outputs and dataframe for US September
CI_susc_pfizer_sept=US_S[[1]];CI_susc_moderna_sept=US_S[[2]];CI_susc_mRNA_sept=US_S[[3]]
CI_susc_hybrid_sept=US_S[[4]];CI_susc_inf_sept=US_S[[5]];CI_susc_boost_mRNA_sept=US_S[[6]]
CI_susc_boost_pfizer_sept=US_S[[7]];VE_susc_pfizer_sept=US_S[[8]];VE_susc_moderna_sept=US_S[[9]]
VE_susc_mRNA_sept=US_S[[10]];VE_susc_hybrid_sept=US_S[[11]];VE_susc_inf_sept=US_S[[12]]
VE_susc_boost_mRNA_sept=US_S[[13]];VE_susc_boost_pfizer_sept=US_S[[14]];CI_infect_pfizer_sept=US_S[[15]]
CI_infect_moderna_sept=US_S[[16]];CI_infect_mRNA_sept=US_S[[17]];CI_infect_hybrid_sept=US_S[[18]]
CI_infect_inf_sept=US_S[[19]];CI_infect_boost_mRNA_sept=US_S[[20]];CI_infect_boost_pfizer_sept=US_S[[21]]
VE_infect_pfizer_sept=US_S[[22]];VE_infect_moderna_sept=US_S[[23]];VE_infect_mRNA_sept=US_S[[24]]
VE_infect_hybrid_sept=US_S[[25]];VE_infect_inf_sept=US_S[[26]];VE_infect_boost_mRNA_sept=US_S[[27]]
VE_infect_boost_pfizer_sept=US_S[[28]]

case_vax_data_sept=US_S[[29]]

V_s_pfi_sept=US_S[[30]]
V_s_mod_sept=US_S[[31]]
V_s_hybrid_sept=US_S[[32]]
V_s_inf_sept=US_S[[33]]
V_s_boost_sept=US_S[[34]]
V_s_boost_pfi_sept=US_S[[35]]

V_i_pfi_sept=US_S[[36]]
V_i_mod_sept=US_S[[37]]
V_i_hybrid_sept=US_S[[38]]
V_i_inf_sept=US_S[[39]]
V_i_boost_sept=US_S[[40]]
V_i_boost_pfi_sept=US_S[[41]]

#Store VE and CI outputs and dataframe for US October
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data$date<=as.Date("2021-10-01"),]))
US_O=VE_at_time(as.Date("2021-10-01"),case_vax_data)

CI_susc_pfizer_oct=US_O[[1]];CI_susc_moderna_oct=US_O[[2]];CI_susc_mRNA_oct=US_O[[3]]
CI_susc_hybrid_oct=US_O[[4]];CI_susc_inf_oct=US_O[[5]];CI_susc_boost_mRNA_oct=US_O[[6]]
CI_susc_boost_pfizer_oct=US_O[[7]];VE_susc_pfizer_oct=US_O[[8]];VE_susc_moderna_oct=US_O[[9]]
VE_susc_mRNA_oct=US_O[[10]];VE_susc_hybrid_oct=US_O[[11]];VE_susc_inf_oct=US_O[[12]]
VE_susc_boost_mRNA_oct=US_O[[13]];VE_susc_boost_pfizer_oct=US_O[[14]];CI_infect_pfizer_oct=US_O[[15]]
CI_infect_moderna_oct=US_O[[16]];CI_infect_mRNA_oct=US_O[[17]];CI_infect_hybrid_oct=US_O[[18]]
CI_infect_inf_oct=US_O[[19]];CI_infect_boost_mRNA_oct=US_O[[20]];CI_infect_boost_pfizer_oct=US_O[[21]]
VE_infect_pfizer_oct=US_O[[22]];VE_infect_moderna_oct=US_O[[23]];VE_infect_mRNA_oct=US_O[[24]]
VE_infect_hybrid_oct=US_O[[25]];VE_infect_inf_oct=US_O[[26]];VE_infect_boost_mRNA_oct=US_O[[27]]
VE_infect_boost_pfizer_oct=US_O[[28]]

case_vax_data_oct=US_O[[29]]

V_s_pfi_oct=US_O[[30]]
V_s_mod_oct=US_O[[31]]
V_s_hybrid_oct=US_O[[32]]
V_s_inf_oct=US_O[[33]]
V_s_boost_oct=US_O[[34]]
V_s_boost_pfi_oct=US_O[[35]]

V_i_pfi_oct=US_O[[36]]
V_i_mod_oct=US_O[[37]]
V_i_hybrid_oct=US_O[[38]]
V_i_inf_oct=US_O[[39]]
V_i_boost_oct=US_O[[40]]
V_i_boost_pfi_oct=US_O[[41]]

#Store VE and CI outputs and dataframe for US November
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data$date<=as.Date("2021-11-01"),]))
US_N=VE_at_time(as.Date("2021-11-01"),case_vax_data)

CI_susc_pfizer_nov=US_N[[1]];CI_susc_moderna_nov=US_N[[2]];CI_susc_mRNA_nov=US_N[[3]]
CI_susc_hybrid_nov=US_N[[4]];CI_susc_inf_nov=US_N[[5]];CI_susc_boost_mRNA_nov=US_N[[6]]
CI_susc_boost_pfizer_nov=US_N[[7]];VE_susc_pfizer_nov=US_N[[8]];VE_susc_moderna_nov=US_N[[9]]
VE_susc_mRNA_nov=US_N[[10]];VE_susc_hybrid_nov=US_N[[11]];VE_susc_inf_nov=US_N[[12]]
VE_susc_boost_mRNA_nov=US_N[[13]];VE_susc_boost_pfizer_nov=US_N[[14]];CI_infect_pfizer_nov=US_N[[15]]
CI_infect_moderna_nov=US_N[[16]];CI_infect_mRNA_nov=US_N[[17]];CI_infect_hybrid_nov=US_N[[18]]
CI_infect_inf_nov=US_N[[19]];CI_infect_boost_mRNA_nov=US_N[[20]];CI_infect_boost_pfizer_nov=US_N[[21]]
VE_infect_pfizer_nov=US_N[[22]];VE_infect_moderna_nov=US_N[[23]];VE_infect_mRNA_nov=US_N[[24]]
VE_infect_hybrid_nov=US_N[[25]];VE_infect_inf_nov=US_N[[26]];VE_infect_boost_mRNA_nov=US_N[[27]]
VE_infect_boost_pfizer_nov=US_N[[28]]

case_vax_data_nov=US_N[[29]]

V_s_pfi_nov=US_N[[30]]
V_s_mod_nov=US_N[[31]]
V_s_hybrid_nov=US_N[[32]]
V_s_inf_nov=US_N[[33]]
V_s_boost_nov=US_N[[34]]
V_s_boost_pfi_nov=US_N[[35]]

V_i_pfi_nov=US_N[[36]]
V_i_mod_nov=US_N[[37]]
V_i_hybrid_nov=US_N[[38]]
V_i_inf_nov=US_N[[39]]
V_i_boost_nov=US_N[[40]]
V_i_boost_pfi_nov=US_N[[41]]

#Store VE and CI outputs and dataframe for US December
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data$date<=as.Date("2021-12-01"),]))
US_D=VE_at_time(as.Date("2021-12-01"),case_vax_data)

CI_susc_pfizer_dec=US_D[[1]];CI_susc_moderna_dec=US_D[[2]];CI_susc_mRNA_dec=US_D[[3]]
CI_susc_hybrid_dec=US_D[[4]];CI_susc_inf_dec=US_D[[5]];CI_susc_boost_mRNA_dec=US_D[[6]]
CI_susc_boost_pfizer_dec=US_D[[7]];VE_susc_pfizer_dec=US_D[[8]];VE_susc_moderna_dec=US_D[[9]]
VE_susc_mRNA_dec=US_D[[10]];VE_susc_hybrid_dec=US_D[[11]];VE_susc_inf_dec=US_D[[12]]
VE_susc_boost_mRNA_dec=US_D[[13]];VE_susc_boost_pfizer_dec=US_D[[14]];CI_infect_pfizer_dec=US_D[[15]]
CI_infect_moderna_dec=US_D[[16]];CI_infect_mRNA_dec=US_D[[17]];CI_infect_hybrid_dec=US_D[[18]]
CI_infect_inf_dec=US_D[[19]];CI_infect_boost_mRNA_dec=US_D[[20]];CI_infect_boost_pfizer_dec=US_D[[21]]
VE_infect_pfizer_dec=US_D[[22]];VE_infect_moderna_dec=US_D[[23]];VE_infect_mRNA_dec=US_D[[24]]
VE_infect_hybrid_dec=US_D[[25]];VE_infect_inf_dec=US_D[[26]];VE_infect_boost_mRNA_dec=US_D[[27]]
VE_infect_boost_pfizer_dec=US_D[[28]]

case_vax_data_dec=US_D[[29]]

V_s_pfi_dec=US_D[[30]]
V_s_mod_dec=US_D[[31]]
V_s_hybrid_dec=US_D[[32]]
V_s_inf_dec=US_D[[33]]
V_s_boost_dec=US_D[[34]]
V_s_boost_pfi_dec=US_D[[35]]

V_i_pfi_dec=US_D[[36]]
V_i_mod_dec=US_D[[37]]
V_i_hybrid_dec=US_D[[38]]
V_i_inf_dec=US_D[[39]]
V_i_boost_dec=US_D[[40]]
V_i_boost_pfi_dec=US_D[[41]]

#California
#Read in data set with daily cases, new fully vaccinated, and deaths w/ row 1 being most recent data
case_vax_data_CA=read.csv("Cases_Vax_CA2.csv") #load case and vax data

#Start date (first day with data) and end date (last day with data), first day vaccine given, and population_CA
start_date_CA = '2020-01-22'; end_date_CA = '2022-01-06';first_vax_CA="2020-12-14";population_CA=39240000
new_end_CA='2021-12-13'

#Add columns with date, number of infections, and proportion of 2nd doses given on a specific day
case_vax_data_CA$date=seq(as.Date(end_date_CA), as.Date(start_date_CA), by = "-1 days")

#Rolling mean and remove NAs
case_vax_data_CA$daily_deaths_sm=rollmean(case_vax_data_CA$Deaths,k=7,fill=NA,align="right")
case_vax_data_CA=case_vax_data_CA[complete.cases(case_vax_data_CA$daily_deaths_sm),] #remove edge rows w/ NAs from rolling means

#Estimate infections from DEATHS by deconvolution #Distributions from Lewnard et al 2020 medRxiv
infec_deaths_CA=deconvolve_incidence(incidence_data = case_vax_data_CA$daily_deaths_sm,
                                     delay=list(inf_shed,shed_sympt,sympt_hosp,hosp_death))

#Daily and cumulative infections
case_vax_data_CA$inf=
  (1*infec_deaths_CA$values[c(rep(NA,-infec_deaths_CA$index_offset),
                              (1):(infec_deaths_CA$index_offset+length(infec_deaths_CA$values)))])/0.003

case_vax_data_CA$cum_inf=NA

case_vax_data_CA=case_vax_data_CA[complete.cases(case_vax_data_CA$inf),]

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_inf[i]=sum(case_vax_data_CA$inf[i:nrow(case_vax_data_CA)])
}

#Fraction total infected (332.9M US Pop 2021)
case_vax_data_CA$frac_tot_inf=case_vax_data_CA$cum_inf/population_CA

#Proportion vaccinated each day
case_vax_data_CA$prop_vax=case_vax_data_CA$Vaccinated/sum(case_vax_data_CA$Vaccinated)

#Create cumulative infected and cumulative fully vaccinated columns
case_vax_data_CA$cum_vax=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_vax[i]=sum(case_vax_data_CA$Vaccinated[i:nrow(case_vax_data_CA)])
}

#Fraction total vaccinated (332.9 US Pop 2021)
case_vax_data_CA$frac_tot_vax=case_vax_data_CA$cum_vax/population_CA

#Create cumulative uninfected and cumulative unvaccinated columns
case_vax_data_CA$cum_unvax=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_unvax[i]=(population_CA-case_vax_data_CA$cum_vax[i])
}

case_vax_data_CA$cum_uninf=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_uninf[i]=(population_CA-case_vax_data_CA$cum_inf[i])
}

#Add case ratios to main dataframe
case_vax_data_CA=cbind(case_vax_data_CA,IRR=c(case_ratios_all$IRR[case_ratios_all$date<=new_end_CA],1))

#Create columns splitting infections between unvaccinated and vaccinated
case_vax_data_CA$daily_inf_in_vax=NA
case_vax_data_CA$daily_inf_in_unvax=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$daily_inf_in_vax=(case_vax_data_CA$inf/(case_vax_data_CA$IRR*(case_vax_data_CA$cum_unvax/case_vax_data_CA$cum_vax)))
}
for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$daily_inf_in_unvax=(case_vax_data_CA$inf-case_vax_data_CA$daily_inf_in_vax)
}

#Enter 0s in place of NAN
case_vax_data_CA$daily_inf_in_vax[is.nan(case_vax_data_CA$daily_inf_in_vax)]<-0
case_vax_data_CA$daily_inf_in_unvax[is.nan(case_vax_data_CA$daily_inf_in_unvax)]<-0

#Create columns splitting unvaccinated in previously positive and naive
case_vax_data_CA$cum_unvax_neg=NA
case_vax_data_CA$cum_unvax_pos=NA

#Prior to vaccines, split unvaccinated into previously positive and naive using infection data
length_prevax_CA=(length(seq(as.Date(start_date_CA), as.Date(first_vax_CA), by = "1 days"))-1)

for (i in (nrow(case_vax_data_CA)-length_prevax_CA+1):nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_unvax_neg[i]=case_vax_data_CA$cum_uninf[i]
  case_vax_data_CA$cum_unvax_pos[i]=case_vax_data_CA$cum_inf[i]
}

#After vaccines available, assume decision to get vaccinated independent of serostatus
#Account for cases in unvaccinated calculated for case ratios
for (i in 0:(nrow(case_vax_data_CA)-length_prevax_CA-1)) {
  case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA-i]=(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+1-i]-case_vax_data_CA$daily_inf_in_unvax[nrow(case_vax_data_CA)-length_prevax_CA-i]-(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+1-i]/(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+1-i]+case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA+1-i])*case_vax_data_CA$Vaccinated[nrow(case_vax_data_CA)-length_prevax_CA-i]))
  case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA-i]=(case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA+1-i]+case_vax_data_CA$daily_inf_in_unvax[nrow(case_vax_data_CA)-length_prevax_CA-i]-(case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA+1-i]/(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+1-i]+case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA+1-i])*case_vax_data_CA$Vaccinated[nrow(case_vax_data_CA)-length_prevax_CA-i]))
}

#Create column with the number of people gaining hybrid immunity each day
case_vax_data_CA$daily_hybrid=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$daily_hybrid[nrow(case_vax_data_CA)+1-i]=(case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)+2-i]/(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)+2-i]+case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)+2-i])*case_vax_data_CA$Vaccinated[nrow(case_vax_data_CA)+1-i]+case_vax_data_CA$daily_inf_in_vax[nrow(case_vax_data_CA)+1-i])
}

case_vax_data_CA$daily_hybrid[is.na(case_vax_data_CA$daily_hybrid)] = 0
case_vax_data_CA$cum_hybrid=NA

for (i in 1:nrow(case_vax_data_CA)) {
  case_vax_data_CA$cum_hybrid[i]=sum(case_vax_data_CA$daily_hybrid[i:nrow(case_vax_data_CA)])
}

#Create column with the daily number of unvaccinated testing positive
case_vax_data_CA$daily_unvax_pos=NA

for (i in (nrow(case_vax_data_CA)-length_prevax_CA+1):nrow(case_vax_data_CA)){
  case_vax_data_CA$daily_unvax_pos[i]=case_vax_data_CA$inf[i]
}

for (i in 1:(nrow(case_vax_data_CA)-length_prevax_CA)){
  case_vax_data_CA$daily_unvax_pos[i]=case_vax_data_CA$daily_inf_in_unvax[i]
}

#Create columns with adjusted number of seronegative people becoming vaccinated
case_vax_data_CA$daily_vax_neg_adj=NA

for (i in 1:(nrow(case_vax_data_CA)-length_prevax_CA)) {
  case_vax_data_CA$daily_vax_neg_adj[nrow(case_vax_data_CA)-length_prevax_CA+1-i]=(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+2-i]/(case_vax_data_CA$cum_unvax_neg[nrow(case_vax_data_CA)-length_prevax_CA+2-i]+case_vax_data_CA$cum_unvax_pos[nrow(case_vax_data_CA)-length_prevax_CA+2-i])*case_vax_data_CA$Vaccinated[nrow(case_vax_data_CA)-length_prevax_CA+1-i]-case_vax_data_CA$prop_vax[nrow(case_vax_data_CA)-length_prevax_CA+1-i]*sum(case_vax_data_CA$daily_inf_in_vax))
}

#Replace NA with 0 thoughout dataframe
case_vax_data_CA[is.na(case_vax_data_CA)] = 0

#Store VE and CI outputs and dataframe for CA September
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_CA$date<=as.Date("2021-09-01"),]))
CA_S=VE_at_time(as.Date("2021-09-01"),case_vax_data_CA)

CI_susc_pfizer_CA_sept=CA_S[[1]];CI_susc_moderna_CA_sept=CA_S[[2]];CI_susc_mRNA_CA_sept=CA_S[[3]]
CI_susc_hybrid_CA_sept=CA_S[[4]];CI_susc_inf_CA_sept=CA_S[[5]];CI_susc_boost_mRNA_CA_sept=CA_S[[6]]
CI_susc_boost_pfizer_CA_sept=CA_S[[7]];VE_susc_pfizer_CA_sept=CA_S[[8]];VE_susc_moderna_CA_sept=CA_S[[9]]
VE_susc_mRNA_CA_sept=CA_S[[10]];VE_susc_hybrid_CA_sept=CA_S[[11]];VE_susc_inf_CA_sept=CA_S[[12]]
VE_susc_boost_mRNA_CA_sept=CA_S[[13]];VE_susc_boost_pfizer_CA_sept=CA_S[[14]];CI_infect_pfizer_CA_sept=CA_S[[15]]
CI_infect_moderna_CA_sept=CA_S[[16]];CI_infect_mRNA_CA_sept=CA_S[[17]];CI_infect_hybrid_CA_sept=CA_S[[18]]
CI_infect_inf_CA_sept=CA_S[[19]];CI_infect_boost_mRNA_CA_sept=CA_S[[20]];CI_infect_boost_pfizer_CA_sept=CA_S[[21]]
VE_infect_pfizer_CA_sept=CA_S[[22]];VE_infect_moderna_CA_sept=CA_S[[23]];VE_infect_mRNA_CA_sept=CA_S[[24]]
VE_infect_hybrid_CA_sept=CA_S[[25]];VE_infect_inf_CA_sept=CA_S[[26]];VE_infect_boost_mRNA_CA_sept=CA_S[[27]]
VE_infect_boost_pfizer_CA_sept=CA_S[[28]]

case_vax_data_CA_sept=CA_S[[29]]

V_s_pfi_CA_sept=CA_S[[30]]
V_s_mod_CA_sept=CA_S[[31]]
V_s_hybrid_CA_sept=CA_S[[32]]
V_s_inf_CA_sept=CA_S[[33]]
V_s_boost_CA_sept=CA_S[[34]]
V_s_boost_pfi_CA_sept=CA_S[[35]]

V_i_pfi_CA_sept=CA_S[[36]]
V_i_mod_CA_sept=CA_S[[37]]
V_i_hybrid_CA_sept=CA_S[[38]]
V_i_inf_CA_sept=CA_S[[39]]
V_i_boost_CA_sept=CA_S[[40]]
V_i_boost_pfi_CA_sept=CA_S[[41]]

#Store VE and CI outputs and dataframe for CA October
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_CA$date<=as.Date("2021-10-01"),]))
CA_O=VE_at_time(as.Date("2021-10-01"),case_vax_data_CA)

CI_susc_pfizer_CA_oct=CA_O[[1]];CI_susc_moderna_CA_oct=CA_O[[2]];CI_susc_mRNA_CA_oct=CA_O[[3]]
CI_susc_hybrid_CA_oct=CA_O[[4]];CI_susc_inf_CA_oct=CA_O[[5]];CI_susc_boost_mRNA_CA_oct=CA_O[[6]]
CI_susc_boost_pfizer_CA_oct=CA_O[[7]];VE_susc_pfizer_CA_oct=CA_O[[8]];VE_susc_moderna_CA_oct=CA_O[[9]]
VE_susc_mRNA_CA_oct=CA_O[[10]];VE_susc_hybrid_CA_oct=CA_O[[11]];VE_susc_inf_CA_oct=CA_O[[12]]
VE_susc_boost_mRNA_CA_oct=CA_O[[13]];VE_susc_boost_pfizer_CA_oct=CA_O[[14]];CI_infect_pfizer_CA_oct=CA_O[[15]]
CI_infect_moderna_CA_oct=CA_O[[16]];CI_infect_mRNA_CA_oct=CA_O[[17]];CI_infect_hybrid_CA_oct=CA_O[[18]]
CI_infect_inf_CA_oct=CA_O[[19]];CI_infect_boost_mRNA_CA_oct=CA_O[[20]];CI_infect_boost_pfizer_CA_oct=CA_O[[21]]
VE_infect_pfizer_CA_oct=CA_O[[22]];VE_infect_moderna_CA_oct=CA_O[[23]];VE_infect_mRNA_CA_oct=CA_O[[24]]
VE_infect_hybrid_CA_oct=CA_O[[25]];VE_infect_inf_CA_oct=CA_O[[26]];VE_infect_boost_mRNA_CA_oct=CA_O[[27]]
VE_infect_boost_pfizer_CA_oct=CA_O[[28]]

case_vax_data_CA_oct=CA_O[[29]]

V_s_pfi_CA_oct=CA_O[[30]]
V_s_mod_CA_oct=CA_O[[31]]
V_s_hybrid_CA_oct=CA_O[[32]]
V_s_inf_CA_oct=CA_O[[33]]
V_s_boost_CA_oct=CA_O[[34]]
V_s_boost_pfi_CA_oct=CA_O[[35]]

V_i_pfi_CA_oct=CA_O[[36]]
V_i_mod_CA_oct=CA_O[[37]]
V_i_hybrid_CA_oct=CA_O[[38]]
V_i_inf_CA_oct=CA_O[[39]]
V_i_boost_CA_oct=CA_O[[40]]
V_i_boost_pfi_CA_oct=CA_O[[41]]

#Store VE and CI outputs and dataframe for CA November
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_CA$date<=as.Date("2021-11-01"),]))
CA_N=VE_at_time(as.Date("2021-11-01"),case_vax_data_CA)

CI_susc_pfizer_CA_nov=CA_N[[1]];CI_susc_moderna_CA_nov=CA_N[[2]];CI_susc_mRNA_CA_nov=CA_N[[3]]
CI_susc_hybrid_CA_nov=CA_N[[4]];CI_susc_inf_CA_nov=CA_N[[5]];CI_susc_boost_mRNA_CA_nov=CA_N[[6]]
CI_susc_boost_pfizer_CA_nov=CA_N[[7]];VE_susc_pfizer_CA_nov=CA_N[[8]];VE_susc_moderna_CA_nov=CA_N[[9]]
VE_susc_mRNA_CA_nov=CA_N[[10]];VE_susc_hybrid_CA_nov=CA_N[[11]];VE_susc_inf_CA_nov=CA_N[[12]]
VE_susc_boost_mRNA_CA_nov=CA_N[[13]];VE_susc_boost_pfizer_CA_nov=CA_N[[14]];CI_infect_pfizer_CA_nov=CA_N[[15]]
CI_infect_moderna_CA_nov=CA_N[[16]];CI_infect_mRNA_CA_nov=CA_N[[17]];CI_infect_hybrid_CA_nov=CA_N[[18]]
CI_infect_inf_CA_nov=CA_N[[19]];CI_infect_boost_mRNA_CA_nov=CA_N[[20]];CI_infect_boost_pfizer_CA_nov=CA_N[[21]]
VE_infect_pfizer_CA_nov=CA_N[[22]];VE_infect_moderna_CA_nov=CA_N[[23]];VE_infect_mRNA_CA_nov=CA_N[[24]]
VE_infect_hybrid_CA_nov=CA_N[[25]];VE_infect_inf_CA_nov=CA_N[[26]];VE_infect_boost_mRNA_CA_nov=CA_N[[27]]
VE_infect_boost_pfizer_CA_nov=CA_N[[28]]

case_vax_data_CA_nov=CA_N[[29]]

V_s_pfi_CA_nov=CA_N[[30]]
V_s_mod_CA_nov=CA_N[[31]]
V_s_hybrid_CA_nov=CA_N[[32]]
V_s_inf_CA_nov=CA_N[[33]]
V_s_boost_CA_nov=CA_N[[34]]
V_s_boost_pfi_CA_nov=CA_N[[35]]

V_i_pfi_CA_nov=CA_N[[36]]
V_i_mod_CA_nov=CA_N[[37]]
V_i_hybrid_CA_nov=CA_N[[38]]
V_i_inf_CA_nov=CA_N[[39]]
V_i_boost_CA_nov=CA_N[[40]]
V_i_boost_pfi_CA_nov=CA_N[[41]]

#Store VE and CI outputs and dataframe for CA December
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_CA$date<=as.Date("2021-12-01"),]))
CA_D=VE_at_time(as.Date("2021-12-01"),case_vax_data_CA)

CI_susc_pfizer_CA_dec=CA_D[[1]];CI_susc_moderna_CA_dec=CA_D[[2]];CI_susc_mRNA_CA_dec=CA_D[[3]]
CI_susc_hybrid_CA_dec=CA_D[[4]];CI_susc_inf_CA_dec=CA_D[[5]];CI_susc_boost_mRNA_CA_dec=CA_D[[6]]
CI_susc_boost_pfizer_CA_dec=CA_D[[7]];VE_susc_pfizer_CA_dec=CA_D[[8]];VE_susc_moderna_CA_dec=CA_D[[9]]
VE_susc_mRNA_CA_dec=CA_D[[10]];VE_susc_hybrid_CA_dec=CA_D[[11]];VE_susc_inf_CA_dec=CA_D[[12]]
VE_susc_boost_mRNA_CA_dec=CA_D[[13]];VE_susc_boost_pfizer_CA_dec=CA_D[[14]];CI_infect_pfizer_CA_dec=CA_D[[15]]
CI_infect_moderna_CA_dec=CA_D[[16]];CI_infect_mRNA_CA_dec=CA_D[[17]];CI_infect_hybrid_CA_dec=CA_D[[18]]
CI_infect_inf_CA_dec=CA_D[[19]];CI_infect_boost_mRNA_CA_dec=CA_D[[20]];CI_infect_boost_pfizer_CA_dec=CA_D[[21]]
VE_infect_pfizer_CA_dec=CA_D[[22]];VE_infect_moderna_CA_dec=CA_D[[23]];VE_infect_mRNA_CA_dec=CA_D[[24]]
VE_infect_hybrid_CA_dec=CA_D[[25]];VE_infect_inf_CA_dec=CA_D[[26]];VE_infect_boost_mRNA_CA_dec=CA_D[[27]]
VE_infect_boost_pfizer_CA_dec=CA_D[[28]]

case_vax_data_CA_dec=CA_D[[29]]

V_s_pfi_CA_dec=CA_D[[30]]
V_s_mod_CA_dec=CA_D[[31]]
V_s_hybrid_CA_dec=CA_D[[32]]
V_s_inf_CA_dec=CA_D[[33]]
V_s_boost_CA_dec=CA_D[[34]]
V_s_boost_pfi_CA_dec=CA_D[[35]]

V_i_pfi_CA_dec=CA_D[[36]]
V_i_mod_CA_dec=CA_D[[37]]
V_i_hybrid_CA_dec=CA_D[[38]]
V_i_inf_CA_dec=CA_D[[39]]
V_i_boost_CA_dec=CA_D[[40]]
V_i_boost_pfi_CA_dec=CA_D[[41]]

#New Zealand
#Read in data set with daily cases, new fully vaccinated, and deaths w/ row 1 being most recent data
case_vax_data_NZ=read.csv("Cases_Vax_NZ2.csv") #load case and vax data

#Start date (first day with data) and end date (last day with data), first day vaccine given, and population_NZ
start_date_NZ = '2020-01-22'; end_date_NZ = '2022-01-06';first_vax_NZ="2021-02-21";population_NZ=5129727 
new_end_NZ='2021-12-13'

#Add columns with date, number of infections, and proportion of 2nd doses given on a specific day
case_vax_data_NZ$date=seq(as.Date(end_date_NZ), as.Date(start_date_NZ), by = "-1 days")

#Estimate infections from DEATHS by deconvolution #Distributions from Lewnard et al 2020 medRxiv
infec_deaths_NZ=deconvolve_incidence(incidence_data = case_vax_data_NZ$Deaths_rolling,
                                     delay=list(inf_shed,shed_sympt,sympt_hosp,hosp_death))
#Daily and cumulative infections
case_vax_data_NZ$inf=
  (1*infec_deaths_NZ$values[c(rep(NA,-infec_deaths_NZ$index_offset),
                              (1):(infec_deaths_NZ$index_offset+length(infec_deaths_NZ$values)))])/0.003


case_vax_data_NZ$cum_inf=NA
case_vax_data_NZ=case_vax_data_NZ[complete.cases(case_vax_data_NZ$Cases),]

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_inf[i]=sum(case_vax_data_NZ$inf[i:nrow(case_vax_data_NZ)])
}

#Fraction total infected (332.9M US Pop 2021)
case_vax_data_NZ$frac_tot_inf=case_vax_data_NZ$cum_inf/population_NZ

#Proportion vaccinated each day
case_vax_data_NZ$prop_vax=case_vax_data_NZ$Vaccinated/sum(case_vax_data_NZ$Vaccinated)

#Create cumulative infected and cumulative fully vaccinated columns
case_vax_data_NZ$cum_vax=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_vax[i]=sum(case_vax_data_NZ$Vaccinated[i:nrow(case_vax_data_NZ)])
}

#Fraction total vaccinated (332.9 US Pop 2021)
case_vax_data_NZ$frac_tot_vax=case_vax_data_NZ$cum_vax/population_NZ

#Create cumulative uninfected and cumulative unvaccinated columns
case_vax_data_NZ$cum_unvax=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_unvax[i]=(population_NZ-case_vax_data_NZ$cum_vax[i])
}

case_vax_data_NZ$cum_uninf=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_uninf[i]=(population_NZ-case_vax_data_NZ$cum_inf[i])
}

#Add case ratios to main dataframe
case_vax_data_NZ=cbind(case_vax_data_NZ,IRR=c(case_ratios_all$IRR[case_ratios_all$date<=new_end_NZ],1))

#Create columns splitting infections between unvaccinated and vaccinated
case_vax_data_NZ$daily_inf_in_vax=NA
case_vax_data_NZ$daily_inf_in_unvax=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$daily_inf_in_vax=(case_vax_data_NZ$inf/(case_vax_data_NZ$IRR*(case_vax_data_NZ$cum_unvax/case_vax_data_NZ$cum_vax)))
}
for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$daily_inf_in_unvax=(case_vax_data_NZ$inf-case_vax_data_NZ$daily_inf_in_vax)
}

#Enter 0s in place of NAN
case_vax_data_NZ$daily_inf_in_vax[is.nan(case_vax_data_NZ$daily_inf_in_vax)]<-0
case_vax_data_NZ$daily_inf_in_unvax[is.nan(case_vax_data_NZ$daily_inf_in_unvax)]<-0

#Create columns splitting unvaccinated in previously positive and naive
case_vax_data_NZ$cum_unvax_neg=NA
case_vax_data_NZ$cum_unvax_pos=NA

#Prior to vaccines, split unvaccinated into previously positive and naive using infection data
length_prevax_NZ=(length(seq(as.Date(start_date_NZ), as.Date(first_vax_NZ), by = "1 days"))-1)

for (i in (nrow(case_vax_data_NZ)-length_prevax_NZ+1):nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_unvax_neg[i]=case_vax_data_NZ$cum_uninf[i]
  case_vax_data_NZ$cum_unvax_pos[i]=case_vax_data_NZ$cum_inf[i]
}

#After vaccines available, assume decision to get vaccinated independent of serostatus
#Account for cases in unvaccinated calculated for case ratios
for (i in 0:(nrow(case_vax_data_NZ)-length_prevax_NZ-1)) {
  case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ-i]=(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]-case_vax_data_NZ$daily_inf_in_unvax[nrow(case_vax_data_NZ)-length_prevax_NZ-i]-(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]/(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]+case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i])*case_vax_data_NZ$Vaccinated[nrow(case_vax_data_NZ)-length_prevax_NZ-i]))
  case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ-i]=(case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]+case_vax_data_NZ$daily_inf_in_unvax[nrow(case_vax_data_NZ)-length_prevax_NZ-i]-(case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]/(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]+case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i])*case_vax_data_NZ$Vaccinated[nrow(case_vax_data_NZ)-length_prevax_NZ-i]))
}

#Create column with the number of people gaining hybrid immunity each day
case_vax_data_NZ$daily_hybrid=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$daily_hybrid[nrow(case_vax_data_NZ)+1-i]=(case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)+2-i]/(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)+2-i]+case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)+2-i])*case_vax_data_NZ$Vaccinated[nrow(case_vax_data_NZ)+1-i]+case_vax_data_NZ$daily_inf_in_vax[nrow(case_vax_data_NZ)+1-i])
}

case_vax_data_NZ$daily_hybrid[is.na(case_vax_data_NZ$daily_hybrid)] = 0
case_vax_data_NZ$cum_hybrid=NA

for (i in 1:nrow(case_vax_data_NZ)) {
  case_vax_data_NZ$cum_hybrid[i]=sum(case_vax_data_NZ$daily_hybrid[i:nrow(case_vax_data_NZ)])
}

#Create column with the daily number of unvaccinated testing positive
case_vax_data_NZ$daily_unvax_pos=NA

for (i in (nrow(case_vax_data_NZ)-length_prevax_NZ+1):nrow(case_vax_data_NZ)){
  case_vax_data_NZ$daily_unvax_pos[i]=case_vax_data_NZ$inf[i]
}

for (i in 1:(nrow(case_vax_data_NZ)-length_prevax_NZ)){
  case_vax_data_NZ$daily_unvax_pos[i]=case_vax_data_NZ$daily_inf_in_unvax[i]
}

#Create columns with adjusted number of seronegative people becoming vaccinated
case_vax_data_NZ$daily_vax_neg_adj=NA

for (i in 1:(nrow(case_vax_data_NZ)-length_prevax_NZ)) {
  case_vax_data_NZ$daily_vax_neg_adj[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]=(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+2-i]/(case_vax_data_NZ$cum_unvax_neg[nrow(case_vax_data_NZ)-length_prevax_NZ+2-i]+case_vax_data_NZ$cum_unvax_pos[nrow(case_vax_data_NZ)-length_prevax_NZ+2-i])*case_vax_data_NZ$Vaccinated[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]-case_vax_data_NZ$prop_vax[nrow(case_vax_data_NZ)-length_prevax_NZ+1-i]*sum(case_vax_data_NZ$daily_inf_in_vax))
}

#Replace NA with 0 thoughout dataframe
case_vax_data_NZ[is.na(case_vax_data_NZ)] = 0

#Store VE and CI outputs and dataframe for NZ September
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_NZ$date<=as.Date("2021-09-01"),]))
NZ_S=VE_at_time(as.Date("2021-09-01"),case_vax_data_NZ)

CI_susc_pfizer_NZ_sept=NZ_S[[1]];CI_susc_moderna_NZ_sept=NZ_S[[2]];CI_susc_mRNA_NZ_sept=NZ_S[[3]]
CI_susc_hybrid_NZ_sept=NZ_S[[4]];CI_susc_inf_NZ_sept=NZ_S[[5]];CI_susc_boost_mRNA_NZ_sept=NZ_S[[6]]
CI_susc_boost_pfizer_NZ_sept=NZ_S[[7]];VE_susc_pfizer_NZ_sept=NZ_S[[8]];VE_susc_moderna_NZ_sept=NZ_S[[9]]
VE_susc_mRNA_NZ_sept=NZ_S[[10]];VE_susc_hybrid_NZ_sept=NZ_S[[11]];VE_susc_inf_NZ_sept=NZ_S[[12]]
VE_susc_boost_mRNA_NZ_sept=NZ_S[[13]];VE_susc_boost_pfizer_NZ_sept=NZ_S[[14]];CI_infect_pfizer_NZ_sept=NZ_S[[15]]
CI_infect_moderna_NZ_sept=NZ_S[[16]];CI_infect_mRNA_NZ_sept=NZ_S[[17]];CI_infect_hybrid_NZ_sept=NZ_S[[18]]
CI_infect_inf_NZ_sept=NZ_S[[19]];CI_infect_boost_mRNA_NZ_sept=NZ_S[[20]];CI_infect_boost_pfizer_NZ_sept=NZ_S[[21]]
VE_infect_pfizer_NZ_sept=NZ_S[[22]];VE_infect_moderna_NZ_sept=NZ_S[[23]];VE_infect_mRNA_NZ_sept=NZ_S[[24]]
VE_infect_hybrid_NZ_sept=NZ_S[[25]];VE_infect_inf_NZ_sept=NZ_S[[26]];VE_infect_boost_mRNA_NZ_sept=NZ_S[[27]]
VE_infect_boost_pfizer_NZ_sept=NZ_S[[28]]

case_vax_data_NZ_sept=NZ_S[[29]]

V_s_pfi_NZ_sept=NZ_S[[30]]
V_s_mod_NZ_sept=NZ_S[[31]]
V_s_hybrid_NZ_sept=NZ_S[[32]]
V_s_inf_NZ_sept=NZ_S[[33]]
V_s_boost_NZ_sept=NZ_S[[34]]
V_s_boost_pfi_NZ_sept=NZ_S[[35]]

V_i_pfi_NZ_sept=NZ_S[[36]]
V_i_mod_NZ_sept=NZ_S[[37]]
V_i_hybrid_NZ_sept=NZ_S[[38]]
V_i_inf_NZ_sept=NZ_S[[39]]
V_i_boost_NZ_sept=NZ_S[[40]]
V_i_boost_pfi_NZ_sept=NZ_S[[41]]

#Store VE and CI outputs and dataframe for NZ October
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_NZ$date<=as.Date("2021-10-01"),]))
NZ_O=VE_at_time(as.Date("2021-10-01"),case_vax_data_NZ)

CI_susc_pfizer_NZ_oct=NZ_O[[1]];CI_susc_moderna_NZ_oct=NZ_O[[2]];CI_susc_mRNA_NZ_oct=NZ_O[[3]]
CI_susc_hybrid_NZ_oct=NZ_O[[4]];CI_susc_inf_NZ_oct=NZ_O[[5]];CI_susc_boost_mRNA_NZ_oct=NZ_O[[6]]
CI_susc_boost_pfizer_NZ_oct=NZ_O[[7]];VE_susc_pfizer_NZ_oct=NZ_O[[8]];VE_susc_moderna_NZ_oct=NZ_O[[9]]
VE_susc_mRNA_NZ_oct=NZ_O[[10]];VE_susc_hybrid_NZ_oct=NZ_O[[11]];VE_susc_inf_NZ_oct=NZ_O[[12]]
VE_susc_boost_mRNA_NZ_oct=NZ_O[[13]];VE_susc_boost_pfizer_NZ_oct=NZ_O[[14]];CI_infect_pfizer_NZ_oct=NZ_O[[15]]
CI_infect_moderna_NZ_oct=NZ_O[[16]];CI_infect_mRNA_NZ_oct=NZ_O[[17]];CI_infect_hybrid_NZ_oct=NZ_O[[18]]
CI_infect_inf_NZ_oct=NZ_O[[19]];CI_infect_boost_mRNA_NZ_oct=NZ_O[[20]];CI_infect_boost_pfizer_NZ_oct=NZ_O[[21]]
VE_infect_pfizer_NZ_oct=NZ_O[[22]];VE_infect_moderna_NZ_oct=NZ_O[[23]];VE_infect_mRNA_NZ_oct=NZ_O[[24]]
VE_infect_hybrid_NZ_oct=NZ_O[[25]];VE_infect_inf_NZ_oct=NZ_O[[26]];VE_infect_boost_mRNA_NZ_oct=NZ_O[[27]]
VE_infect_boost_pfizer_NZ_oct=NZ_O[[28]]

case_vax_data_NZ_oct=NZ_O[[29]]

V_s_pfi_NZ_oct=NZ_O[[30]]
V_s_mod_NZ_oct=NZ_O[[31]]
V_s_hybrid_NZ_oct=NZ_O[[32]]
V_s_inf_NZ_oct=NZ_O[[33]]
V_s_boost_NZ_oct=NZ_O[[34]]
V_s_boost_pfi_NZ_oct=NZ_O[[35]]

V_i_pfi_NZ_oct=NZ_O[[36]]
V_i_mod_NZ_oct=NZ_O[[37]]
V_i_hybrid_NZ_oct=NZ_O[[38]]
V_i_inf_NZ_oct=NZ_O[[39]]
V_i_boost_NZ_oct=NZ_O[[40]]
V_i_boost_pfi_NZ_oct=NZ_O[[41]]

#Store VE and CI outputs and dataframe for NZ November
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_NZ$date<=as.Date("2021-11-01"),]))
NZ_N=VE_at_time(as.Date("2021-11-01"),case_vax_data_NZ)

CI_susc_pfizer_NZ_nov=NZ_N[[1]];CI_susc_moderna_NZ_nov=NZ_N[[2]];CI_susc_mRNA_NZ_nov=NZ_N[[3]]
CI_susc_hybrid_NZ_nov=NZ_N[[4]];CI_susc_inf_NZ_nov=NZ_N[[5]];CI_susc_boost_mRNA_NZ_nov=NZ_N[[6]]
CI_susc_boost_pfizer_NZ_nov=NZ_N[[7]];VE_susc_pfizer_NZ_nov=NZ_N[[8]];VE_susc_moderna_NZ_nov=NZ_N[[9]]
VE_susc_mRNA_NZ_nov=NZ_N[[10]];VE_susc_hybrid_NZ_nov=NZ_N[[11]];VE_susc_inf_NZ_nov=NZ_N[[12]]
VE_susc_boost_mRNA_NZ_nov=NZ_N[[13]];VE_susc_boost_pfizer_NZ_nov=NZ_N[[14]];CI_infect_pfizer_NZ_nov=NZ_N[[15]]
CI_infect_moderna_NZ_nov=NZ_N[[16]];CI_infect_mRNA_NZ_nov=NZ_N[[17]];CI_infect_hybrid_NZ_nov=NZ_N[[18]]
CI_infect_inf_NZ_nov=NZ_N[[19]];CI_infect_boost_mRNA_NZ_nov=NZ_N[[20]];CI_infect_boost_pfizer_NZ_nov=NZ_N[[21]]
VE_infect_pfizer_NZ_nov=NZ_N[[22]];VE_infect_moderna_NZ_nov=NZ_N[[23]];VE_infect_mRNA_NZ_nov=NZ_N[[24]]
VE_infect_hybrid_NZ_nov=NZ_N[[25]];VE_infect_inf_NZ_nov=NZ_N[[26]];VE_infect_boost_mRNA_NZ_nov=NZ_N[[27]]
VE_infect_boost_pfizer_NZ_nov=NZ_N[[28]]

case_vax_data_NZ_nov=NZ_N[[29]]

V_s_pfi_NZ_nov=NZ_N[[30]]
V_s_mod_NZ_nov=NZ_N[[31]]
V_s_hybrid_NZ_nov=NZ_N[[32]]
V_s_inf_NZ_nov=NZ_N[[33]]
V_s_boost_NZ_nov=NZ_N[[34]]
V_s_boost_pfi_NZ_nov=NZ_N[[35]]

V_i_pfi_NZ_nov=NZ_N[[36]]
V_i_mod_NZ_nov=NZ_N[[37]]
V_i_hybrid_NZ_nov=NZ_N[[38]]
V_i_inf_NZ_nov=NZ_N[[39]]
V_i_boost_NZ_nov=NZ_N[[40]]
V_i_boost_pfi_NZ_nov=NZ_N[[41]]

#Store VE and CI outputs and dataframe for NZ December
Day=c(rep(0,15),1:nrow(case_vax_data[case_vax_data_NZ$date<=as.Date("2021-12-01"),]))
NZ_D=VE_at_time(as.Date("2021-12-01"),case_vax_data_NZ)

CI_susc_pfizer_NZ_dec=NZ_D[[1]];CI_susc_moderna_NZ_dec=NZ_D[[2]];CI_susc_mRNA_NZ_dec=NZ_D[[3]]
CI_susc_hybrid_NZ_dec=NZ_D[[4]];CI_susc_inf_NZ_dec=NZ_D[[5]];CI_susc_boost_mRNA_NZ_dec=NZ_D[[6]]
CI_susc_boost_pfizer_NZ_dec=NZ_D[[7]];VE_susc_pfizer_NZ_dec=NZ_D[[8]];VE_susc_moderna_NZ_dec=NZ_D[[9]]
VE_susc_mRNA_NZ_dec=NZ_D[[10]];VE_susc_hybrid_NZ_dec=NZ_D[[11]];VE_susc_inf_NZ_dec=NZ_D[[12]]
VE_susc_boost_mRNA_NZ_dec=NZ_D[[13]];VE_susc_boost_pfizer_NZ_dec=NZ_D[[14]];CI_infect_pfizer_NZ_dec=NZ_D[[15]]
CI_infect_moderna_NZ_dec=NZ_D[[16]];CI_infect_mRNA_NZ_dec=NZ_D[[17]];CI_infect_hybrid_NZ_dec=NZ_D[[18]]
CI_infect_inf_NZ_dec=NZ_D[[19]];CI_infect_boost_mRNA_NZ_dec=NZ_D[[20]];CI_infect_boost_pfizer_NZ_dec=NZ_D[[21]]
VE_infect_pfizer_NZ_dec=NZ_D[[22]];VE_infect_moderna_NZ_dec=NZ_D[[23]];VE_infect_mRNA_NZ_dec=NZ_D[[24]]
VE_infect_hybrid_NZ_dec=NZ_D[[25]];VE_infect_inf_NZ_dec=NZ_D[[26]];VE_infect_boost_mRNA_NZ_dec=NZ_D[[27]]
VE_infect_boost_pfizer_NZ_dec=NZ_D[[28]]

case_vax_data_NZ_dec=NZ_D[[29]]

V_s_pfi_NZ_dec=NZ_D[[30]]
V_s_mod_NZ_dec=NZ_D[[31]]
V_s_hybrid_NZ_dec=NZ_D[[32]]
V_s_inf_NZ_dec=NZ_D[[33]]
V_s_boost_NZ_dec=NZ_D[[34]]
V_s_boost_pfi_NZ_dec=NZ_D[[35]]

V_i_pfi_NZ_dec=NZ_D[[36]]
V_i_mod_NZ_dec=NZ_D[[37]]
V_i_hybrid_NZ_dec=NZ_D[[38]]
V_i_inf_NZ_dec=NZ_D[[39]]
V_i_boost_NZ_dec=NZ_D[[40]]
V_i_boost_pfi_NZ_dec=NZ_D[[41]]


#Rt function 
Rt=function(parsR){
  with(parsR,
       R0*((f_n*(1-f_pu))+(f_n*f_pu)*(1-VE_si)*(1-VE_ii)+(1-f_pv)*(f_v-f_b)*(1-VE_iv)*(1-VE_iv)+(f_pv)*(f_v-f_b)*(1-VE_is)*(1-VE_is)+(f_b)*(1-VE_ib)*(1-VE_ib)))}

#parameter values for scenario A - R0=3.2;  59% vacc., 81% prev. inf.
parsA_dec=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_dec$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_dec$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_dec$cum_unvax_pos[1]/(case_vax_data_dec$cum_unvax_pos[1]+case_vax_data_dec$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_dec$cum_hybrid)/max(case_vax_data_dec$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_dec, #prev infected waned Inf
  VE_ii=V_i_inf_dec, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_dec)+mRNA_moderna_ratio*(V_s_mod_dec), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_dec)+mRNA_moderna_ratio*(V_i_mod_dec), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_dec, #hybrid VEi
  VE_is=V_i_hybrid_dec, #hybrid VEt
  VE_sb=V_s_boost_dec, #boost inf
  VE_ib=V_i_boost_dec #boost trans
  
)
pars1A_dec=data.frame(f_b=seq(0.00,round(max(case_vax_data_dec$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsA_dec, f_b=seq(0.00,round(max(case_vax_data_dec$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=3.2;  59% vacc., 81% prev. inf.")

#Scenario B - 100% vacc., R0 = 7, 81% prev. inf.
parsB_dec=parsA_dec; parsB_dec[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,parsA_dec$f_n*parsA_dec$f_pu+parsA_dec$f_v*parsA_dec$f_pv)
pars1B_dec=data.frame(f_b=seq(0,1,by=0.001),Rtb=Rt(cbind(parsB_dec, f_b=seq(0,1,by=0.001))),Scenario="R0=7.0;  100% vacc., 81% prev. inf.")

#Scenario C - R0=7.0;  59% vacc., 81% prev. inf.
parsC_dec=parsA_dec;parsC_dec["R0"]=7
pars1C_dec=data.frame(f_b=seq(0.00,round(max(case_vax_data_dec$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsC_dec, f_b=seq(0.00,round(max(case_vax_data_dec$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  59% vacc., 81% prev. inf.")

# #Scenario D - 65% vacc., R0 3.7, Prev inf = 66%
parsD_dec=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_CA_dec$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_CA_dec$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_CA_dec$cum_unvax_pos[1]/(case_vax_data_CA_dec$cum_unvax_pos[1]+case_vax_data_CA_dec$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_CA_dec$cum_hybrid)/max(case_vax_data_CA_dec$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_CA_dec, #prev infected waned Inf
  VE_ii=V_i_inf_CA_dec, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_CA_dec)+mRNA_moderna_ratio*(V_s_mod_CA_dec), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_CA_dec)+mRNA_moderna_ratio*(V_i_mod_CA_dec), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_CA_dec, #hybrid VEi
  VE_is=V_i_hybrid_CA_dec, #hybrid VEt
  VE_sb=V_s_boost_CA_dec, #boost inf
  VE_ib=V_i_boost_CA_dec #boost trans
)
pars1D_dec=data.frame(f_b=seq(0.00,round(max(case_vax_data_CA_dec$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsD_dec, f_b=seq(0.00,round(max(case_vax_data_CA_dec$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=3.2;  65% vacc., 66% prev. inf.")

parsD_dec$f_n*parsD_dec$f_pu+parsD_dec$f_v*parsD_dec$f_pv

# #Scenario E - 71% vacc., R0 3.7, Prev inf = 0.3%
parsE_dec=data.frame(
  R0=7.0,
  f_n=1-max(case_vax_data_NZ_dec$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_NZ_dec$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_NZ_dec$cum_unvax_pos[1]/(case_vax_data_NZ_dec$cum_unvax_pos[1]+case_vax_data_NZ_dec$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_NZ_dec$cum_hybrid)/max(case_vax_data_NZ_dec$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_NZ_dec, #prev infected waned Inf
  VE_ii=V_i_inf_NZ_dec, #prev infected waned trans
  VE_sv=V_s_pfi_NZ_dec, #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=V_i_pfi_NZ_dec, #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_NZ_dec, #hybrid VEi
  VE_is=V_i_hybrid_NZ_dec, #hybrid VEt
  VE_sb=V_s_boost_NZ_dec, #boost inf
  VE_ib=V_i_boost_NZ_dec #boost trans
)
pars1E_dec=data.frame(f_b=seq(0.00,round(max(case_vax_data_NZ_dec$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsE_dec, f_b=seq(0.00,round(max(case_vax_data_NZ_dec$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  71% vacc., 0.3% prev. inf.")


# #Generate CIs for each scenario fitted Rt line
nd1=10000
draws=data.frame() #generate 10,000 draws from each estimate on logit scale, then transform
draws_dec=cbind("VE_si"=unlist(CI_susc_inf_dec),
                "VE_ii"=unlist(CI_infect_inf_dec),
                "VE_sv"=unlist(CI_susc_mRNA_dec),
                "VE_iv"=unlist(CI_infect_mRNA_dec),
                "VE_ss"=unlist(CI_susc_hybrid_dec),
                "VE_is"=unlist(CI_infect_hybrid_dec),
                "VE_sb"=unlist(CI_susc_boost_mRNA_dec),
                "VE_ib"=unlist(CI_infect_boost_mRNA_dec))

draws_CA_dec=cbind("VE_si"=unlist(CI_susc_inf_CA_dec),
                   "VE_ii"=unlist(CI_infect_inf_CA_dec),
                   "VE_sv"=unlist(CI_susc_mRNA_CA_dec),
                   "VE_iv"=unlist(CI_infect_mRNA_CA_dec),
                   "VE_ss"=unlist(CI_susc_hybrid_CA_dec),
                   "VE_is"=unlist(CI_infect_hybrid_CA_dec),
                   "VE_sb"=unlist(CI_susc_boost_mRNA_CA_dec),
                   "VE_ib"=unlist(CI_infect_boost_mRNA_CA_dec))

draws_NZ_dec=cbind("VE_si"=unlist(CI_susc_inf_NZ_dec),
                   "VE_ii"=unlist(CI_infect_inf_NZ_dec),
                   "VE_sv"=unlist(CI_susc_pfizer_NZ_dec),
                   "VE_iv"=unlist(CI_infect_pfizer_NZ_dec),
                   "VE_ss"=unlist(CI_susc_hybrid_NZ_dec),
                   "VE_is"=unlist(CI_infect_hybrid_NZ_dec),
                   "VE_sb"=unlist(CI_susc_boost_pfizer_NZ_dec),
                   "VE_ib"=unlist(CI_infect_boost_pfizer_NZ_dec))


for (i in 1:1001) { #calculate CIs for each scenario
  if (i<596) { #2 scenarios that have max boosters = 59.4%
    pars1A_dec[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsA_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_dec,f_b=pars1A_dec$f_b[i])))
    pars1C_dec[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsC_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_dec,f_b=pars1C_dec$f_b[i])))
  }
  if(i<653){ #scenario that has max booster = 71%
    pars1D_dec[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsD_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_CA_dec,f_b=pars1D_dec$f_b[i])))
  }
  if (i <711) {
    pars1E_dec[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsE_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_NZ_dec,f_b=pars1E_dec$f_b[i])))
    
  }
  pars1B_dec[i,c("lower","upper")]=
    quantile(probs = c(0.16,0.84),x=Rt(cbind(parsB_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_dec,f_b=pars1B_dec$f_b[i])))
}


# #group data into single dataframe
Rt_dat_dec=data.frame(rbind(pars1A_dec,pars1B_dec,pars1C_dec,pars1D_dec,pars1E_dec))#,pars1E,pars1F)) 

rt_ratCI_dec=data.frame(month="Dec",lower=NA,upper=NA,
                        rt_rat=pars1A_dec$Rtb[pars1A_dec$f_b==max(pars1A_dec$f_b)]/pars1A_dec$Rtb,f_b=pars1A_dec$f_b,f_bmax=max(pars1A_dec$f_b))
for (i in 1:nrow(rt_ratCI_dec)) {
  rt_ratCI_dec[i,c("lower","upper")]=quantile(probs = c(0.025,0.975),
                                              x=Rt(cbind(parsA_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_dec,f_b=max(pars1A_dec$f_b)))
                                              /Rt(cbind(parsA_dec[c("R0","f_n","f_v","f_pu","f_pv")],draws_dec,f_b=(i-1)/1000)))
}


##November
#parameter values for scenario A - R0=3.2;  58% vacc., 76% prev. inf.
#pars A
parsA_nov=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_nov$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_nov$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_nov$cum_unvax_pos[1]/(case_vax_data_nov$cum_unvax_pos[1]+case_vax_data_nov$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_nov$cum_hybrid)/max(case_vax_data_nov$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_nov, #prev infected waned Inf
  VE_ii=V_i_inf_nov, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_nov)+mRNA_moderna_ratio*(V_s_mod_nov), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_nov)+mRNA_moderna_ratio*(V_i_mod_nov), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_nov, #hybrid VEi
  VE_is=V_i_hybrid_nov, #hybrid VEt
  VE_sb=V_s_boost_nov, #boost inf
  VE_ib=V_i_boost_nov #boost trans
)
pars1A_nov=data.frame(f_b=seq(0.00,round(max(case_vax_data_nov$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsA_nov, f_b=seq(0.00,round(max(case_vax_data_nov$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=3.2;  58% vacc., 76% prev. inf.")

#Scenario B - 100% vacc., R0 = 7, Prev inf = 76%
parsB_nov=parsA_nov; parsB_nov[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,parsA_nov$f_n*parsA_nov$f_pu+parsA_nov$f_v*parsA_nov$f_pv)
pars1B_nov=data.frame(f_b=seq(0,1,by=0.001),Rtb=Rt(cbind(parsB_nov, f_b=seq(0,1,by=0.001))),Scenario="R0=7.0;  100% vacc., 54.9% prev. inf.")

#Scenario C - 58% vacc., R0 = 7, Prev inf = 76%
parsC_nov=parsA_nov;parsC_nov["R0"]=7
pars1C_nov=data.frame(f_b=seq(0.00,round(max(case_vax_data_nov$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsC_nov, f_b=seq(0.00,round(max(case_vax_data_nov$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  58% vacc., 76% prev. inf.")

# #Scenario D - 63% vacc., R0 3.2, Prev inf = 64%
# #Scenario D - 63% vacc., R0 3.2, Prev inf = 64%
parsD_nov=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_CA_nov$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_CA_nov$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_CA_nov$cum_unvax_pos[1]/(case_vax_data_CA_nov$cum_unvax_pos[1]+case_vax_data_CA_nov$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_CA_nov$cum_hybrid)/max(case_vax_data_CA_nov$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_CA_nov, #prev infected waned Inf
  VE_ii=V_i_inf_CA_nov, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_CA_nov)+mRNA_moderna_ratio*(V_s_mod_CA_nov), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_CA_nov)+mRNA_moderna_ratio*(V_i_mod_CA_nov), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_CA_nov, #hybrid VEi
  VE_is=V_i_hybrid_CA_nov, #hybrid VEt
  VE_sb=V_s_boost_CA_nov, #boost inf
  VE_ib=V_i_boost_CA_nov #boost trans
)
pars1D_nov=data.frame(f_b=seq(0.00,round(max(case_vax_data_CA_nov$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsD_nov, f_b=seq(0.00,round(max(case_vax_data_CA_nov$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  63% vacc., 64% prev. inf.")

parsD_nov$f_n*parsD_nov$f_pu+parsD_nov$f_v*parsD_nov$f_pv

# #Scenario E - 62% vacc., R0 7.0, Prev inf = 0.2%
parsE_nov=data.frame(
  R0=7.0,
  f_n=1-max(case_vax_data_NZ_nov$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_NZ_nov$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_NZ_nov$cum_unvax_pos[1]/(case_vax_data_NZ_nov$cum_unvax_pos[1]+case_vax_data_NZ_nov$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_NZ_nov$cum_hybrid)/max(case_vax_data_NZ_nov$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_NZ_nov, #prev infected waned Inf
  VE_ii=V_i_inf_NZ_nov, #prev infected waned trans
  VE_sv=V_s_pfi_NZ_nov, #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=V_i_pfi_NZ_nov, #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_NZ_nov, #hybrid VEi
  VE_is=V_i_hybrid_NZ_nov, #hybrid VEt
  VE_sb=V_s_boost_NZ_nov, #boost inf
  VE_ib=V_i_boost_NZ_nov #boost trans
)
pars1E_nov=data.frame(f_b=seq(0.00,round(max(case_vax_data_NZ_nov$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsE_nov, f_b=seq(0.00,round(max(case_vax_data_NZ_nov$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  62% vacc., 0.2% prev. inf.")


# #Generate CIs for each scenario fitted Rt line
nd1=10000
draws=data.frame() #generate 10,000 draws from each estimate on logit scale, then transform
draws_nov=cbind("VE_si"=unlist(CI_susc_inf_nov),
                "VE_ii"=unlist(CI_infect_inf_nov),
                "VE_sv"=unlist(CI_susc_mRNA_nov),
                "VE_iv"=unlist(CI_infect_mRNA_nov),
                "VE_ss"=unlist(CI_susc_hybrid_nov),
                "VE_is"=unlist(CI_infect_hybrid_nov),
                "VE_sb"=unlist(CI_susc_boost_mRNA_nov),
                "VE_ib"=unlist(CI_infect_boost_mRNA_nov))

draws_CA_nov=cbind("VE_si"=unlist(CI_susc_inf_CA_nov),
                   "VE_ii"=unlist(CI_infect_inf_CA_nov),
                   "VE_sv"=unlist(CI_susc_mRNA_CA_nov),
                   "VE_iv"=unlist(CI_infect_mRNA_CA_nov),
                   "VE_ss"=unlist(CI_susc_hybrid_CA_nov),
                   "VE_is"=unlist(CI_infect_hybrid_CA_nov),
                   "VE_sb"=unlist(CI_susc_boost_mRNA_CA_nov),
                   "VE_ib"=unlist(CI_infect_boost_mRNA_CA_nov))

draws_NZ_nov=cbind("VE_si"=unlist(CI_susc_inf_NZ_nov),
                   "VE_ii"=unlist(CI_infect_inf_NZ_nov),
                   "VE_sv"=unlist(CI_susc_pfizer_NZ_nov),
                   "VE_iv"=unlist(CI_infect_pfizer_NZ_nov),
                   "VE_ss"=unlist(CI_susc_hybrid_NZ_nov),
                   "VE_is"=unlist(CI_infect_hybrid_NZ_nov),
                   "VE_sb"=unlist(CI_susc_boost_pfizer_NZ_nov),
                   "VE_ib"=unlist(CI_infect_boost_pfizer_NZ_nov))


for (i in 1:1001) { #calculate CIs for each scenario
  if (i<581) { #2 scenarios that have max boosters = 60.4%
    pars1A_nov[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsA_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_nov,f_b=pars1A_nov$f_b[i])))
    pars1C_nov[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsC_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_nov,f_b=pars1C_nov$f_b[i])))
  }
  if (i <623) {
    pars1E_nov[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsE_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_NZ_nov,f_b=pars1E_nov$f_b[i])))
    
  }
  if(i<635){ #scenario that has max booster = 75%
    pars1D_nov[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsD_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_CA_nov,f_b=pars1D_nov$f_b[i])))
  }
  
  pars1B_nov[i,c("lower","upper")]=
    quantile(probs = c(0.16,0.84),x=Rt(cbind(parsB_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_nov,f_b=pars1B_nov$f_b[i])))
}

# #group data into single dataframe
Rt_dat_nov=data.frame(rbind(pars1A_nov,pars1B_nov,pars1C_nov,pars1D_nov,pars1E_nov))#,pars1E,pars1F)) 

rt_ratCI_nov=data.frame(month="Nov",lower=NA,upper=NA,
                        rt_rat=pars1A_nov$Rtb[pars1A_nov$f_b==max(pars1A_nov$f_b)]/pars1A_nov$Rtb,f_b=pars1A_nov$f_b,f_bmax=max(pars1A_nov$f_b))


for (i in 1:nrow(rt_ratCI_nov)) {
  rt_ratCI_nov[i,c("lower","upper")]=quantile(probs = c(0.025,0.975),
                                              x=Rt(cbind(parsA_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_nov,f_b=max(pars1A_nov$f_b)))
                                              /Rt(cbind(parsA_nov[c("R0","f_n","f_v","f_pu","f_pv")],draws_nov,f_b=(i-1)/1000)))
}

##October
#parameter values for scenario A - R0=3.2;  56% vacc., 73% prev. inf.
#pars A
parsA_oct=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_oct$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_oct$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_oct$cum_unvax_pos[1]/(case_vax_data_oct$cum_unvax_pos[1]+case_vax_data_oct$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_oct$cum_hybrid)/max(case_vax_data_oct$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_oct, #prev infected waned Inf
  VE_ii=V_i_inf_oct, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_oct)+mRNA_moderna_ratio*(V_s_mod_oct), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_oct)+mRNA_moderna_ratio*(V_i_mod_oct), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_oct, #hybrid VEi
  VE_is=V_i_hybrid_oct, #hybrid VEt
  VE_sb=V_s_boost_oct, #boost inf
  VE_ib=V_i_boost_oct #boost trans
)
pars1A_oct=data.frame(f_b=seq(0.00,round(max(case_vax_data_oct$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsA_oct, f_b=seq(0.00,round(max(case_vax_data_oct$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=3.2;  56% vacc., 73% prev. inf.")

#Scenario B - 100% vacc., R0 = 7, Prev inf = 73%
parsB_oct=parsA_oct; parsB_oct[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,parsA_oct$f_n*parsA_oct$f_pu+parsA_oct$f_v*parsA_oct$f_pv)
pars1B_oct=data.frame(f_b=seq(0,1,by=0.001),Rtb=Rt(cbind(parsB_oct, f_b=seq(0,1,by=0.001))),Scenario="R0=7.0;  100% vacc., 52.0% prev. inf.")

#Scenario C - 56% vacc., R0 = 7, Prev inf = 73%
parsC_oct=parsA_oct;parsC_oct["R0"]=7
pars1C_oct=data.frame(f_b=seq(0.00,round(max(case_vax_data_oct$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsC_oct, f_b=seq(0.00,round(max(case_vax_data_oct$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  56% vacc., 73% prev. inf.")

# #Scenario D - 62% vacc., R0 3.2, Prev inf = 62%
parsD_oct=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_CA_oct$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_CA_oct$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_CA_oct$cum_unvax_pos[1]/(case_vax_data_CA_oct$cum_unvax_pos[1]+case_vax_data_CA_oct$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_CA_oct$cum_hybrid)/max(case_vax_data_CA_oct$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_CA_oct, #prev infected waned Inf
  VE_ii=V_i_inf_CA_oct, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_CA_oct)+mRNA_moderna_ratio*(V_s_mod_CA_oct), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_CA_oct)+mRNA_moderna_ratio*(V_i_mod_CA_oct), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_CA_oct, #hybrid VEi
  VE_is=V_i_hybrid_CA_oct, #hybrid VEt
  VE_sb=V_s_boost_CA_oct, #boost inf
  VE_ib=V_i_boost_CA_oct #boost trans
)
pars1D_oct=data.frame(f_b=seq(0.00,round(max(case_vax_data_CA_oct$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsD_oct, f_b=seq(0.00,round(max(case_vax_data_CA_oct$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=3.2;  62% vacc., 62% prev. inf.")

parsD_oct$f_n*parsD_oct$f_pu+parsD_oct$f_v*parsD_oct$f_pv

parsE_oct=data.frame(
  R0=7.0,
  f_n=1-max(case_vax_data_NZ_oct$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_NZ_oct$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_NZ_oct$cum_unvax_pos[1]/(case_vax_data_NZ_oct$cum_unvax_pos[1]+case_vax_data_NZ_oct$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_NZ_oct$cum_hybrid)/max(case_vax_data_NZ_oct$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_NZ_oct, #prev infected waned Inf
  VE_ii=V_i_inf_NZ_oct, #prev infected waned trans
  VE_sv=V_s_pfi_NZ_oct, #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=V_i_pfi_NZ_oct, #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_NZ_oct, #hybrid VEi
  VE_is=V_i_hybrid_NZ_oct, #hybrid VEt
  VE_sb=V_s_boost_NZ_oct, #boost inf
  VE_ib=V_i_boost_NZ_oct #boost trans
)
pars1E_oct=data.frame(f_b=seq(0.00,round(max(case_vax_data_NZ_oct$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsE_oct, f_b=seq(0.00,round(max(case_vax_data_NZ_oct$frac_tot_vax), digits = 3),by=0.001))),Scenario="R0=7.0;  38% vacc., 0.2% prev. inf.")


# #Generate CIs for each scenario fitted Rt line
nd1=10000
draws=data.frame() #generate 10,000 draws from each estimate on logit scale, then transform
draws_oct=cbind("VE_si"=unlist(CI_susc_inf_oct),
                "VE_ii"=unlist(CI_infect_inf_oct),
                "VE_sv"=unlist(CI_susc_mRNA_oct),
                "VE_iv"=unlist(CI_infect_mRNA_oct),
                "VE_ss"=unlist(CI_susc_hybrid_oct),
                "VE_is"=unlist(CI_infect_hybrid_oct),
                "VE_sb"=unlist(CI_susc_boost_mRNA_oct),
                "VE_ib"=unlist(CI_infect_boost_mRNA_oct))


draws_CA_oct=cbind("VE_si"=unlist(CI_susc_inf_CA_oct),
                   "VE_ii"=unlist(CI_infect_inf_CA_oct),
                   "VE_sv"=unlist(CI_susc_mRNA_CA_oct),
                   "VE_iv"=unlist(CI_infect_mRNA_CA_oct),
                   "VE_ss"=unlist(CI_susc_hybrid_CA_oct),
                   "VE_is"=unlist(CI_infect_hybrid_CA_oct),
                   "VE_sb"=unlist(CI_susc_boost_mRNA_CA_oct),
                   "VE_ib"=unlist(CI_infect_boost_mRNA_CA_oct))

draws_NZ_oct=cbind("VE_si"=unlist(CI_susc_inf_NZ_oct),
                   "VE_ii"=unlist(CI_infect_inf_NZ_oct),
                   "VE_sv"=unlist(CI_susc_pfizer_NZ_oct),
                   "VE_iv"=unlist(CI_infect_pfizer_NZ_oct),
                   "VE_ss"=unlist(CI_susc_hybrid_NZ_oct),
                   "VE_is"=unlist(CI_infect_hybrid_NZ_oct),
                   "VE_sb"=unlist(CI_susc_boost_pfizer_NZ_oct),
                   "VE_ib"=unlist(CI_infect_boost_pfizer_NZ_oct))

for (i in 1:1001) { #calculate CIs for each scenario
  if (i <385) {
    pars1E_oct[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsE_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_NZ_oct,f_b=pars1E_oct$f_b[i])))
  }
  if (i<563) { #2 scenarios that have max boosters = 60.4%
    pars1A_oct[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsA_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_oct,f_b=pars1A_oct$f_b[i])))
    pars1C_oct[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsC_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_oct,f_b=pars1C_oct$f_b[i])))
  }
  if(i<616){ #scenario that has max booster = 75%
    pars1D_oct[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsD_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_CA_oct,f_b=pars1D_oct$f_b[i])))
  }
  pars1B_oct[i,c("lower","upper")]=
    quantile(probs = c(0.16,0.84),x=Rt(cbind(parsB_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_oct,f_b=pars1B_oct$f_b[i])))
}

# #group data into single dataframe
Rt_dat_oct=data.frame(rbind(pars1A_oct,pars1B_oct,pars1C_oct,pars1D_oct,pars1E_oct))#,pars1E,pars1F)) 

rt_ratCI_oct=data.frame(month="Oct",lower=NA,upper=NA,
                        rt_rat=pars1A_oct$Rtb[pars1A_oct$f_b==max(pars1A_oct$f_b)]/pars1A_oct$Rtb,f_b=pars1A_oct$f_b,f_bmax=max(pars1A_oct$f_b))


for (i in 1:nrow(rt_ratCI_oct)) {
  rt_ratCI_oct[i,c("lower","upper")]=quantile(probs = c(0.025,0.975),
                                              x=Rt(cbind(parsA_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_oct,f_b=max(pars1A_oct$f_b)))
                                              /Rt(cbind(parsA_oct[c("R0","f_n","f_v","f_pu","f_pv")],draws_oct,f_b=(i-1)/1000)))
}

##September
#Rt function 
Rt=function(parsR){
  with(parsR,
       R0*((f_n*(1-f_pu))+(f_n*f_pu)*(1-VE_si)*(1-VE_ii)+(1-f_pv)*(f_v-f_b)*(1-VE_sv)*(1-VE_iv)+(f_pv)*(f_v-f_b)*(1-VE_ss)*(1-VE_is)+(f_b)*(1-VE_sb)*(1-VE_ib)))}

#parameter values for scenario A - R0=3.2;  53% vacc., 67% prev. inf.
#pars A
parsA_sept=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_sept$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_sept$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_sept$cum_unvax_pos[1]/(case_vax_data_sept$cum_unvax_pos[1]+case_vax_data_sept$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_sept$cum_hybrid)/max(case_vax_data_sept$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_sept, #prev infected waned Inf
  VE_ii=V_i_inf_sept, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_sept)+mRNA_moderna_ratio*(V_s_mod_sept), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_sept)+mRNA_moderna_ratio*(V_i_mod_sept), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_sept, #hybrid VEi
  VE_is=V_i_hybrid_sept, #hybrid VEt
  VE_sb=V_s_boost_sept, #boost inf
  VE_ib=V_i_boost_sept #boost trans
)
pars1A_sept=data.frame(f_b=seq(0.00,round(max(case_vax_data_sept$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsA_sept, f_b=seq(0.00,round(max(case_vax_data_sept$frac_tot_vax), digits = 3),by=0.001))),Scenario="US/D:           R0=3.2, V: 53%, I: 67%")

#Scenario B - 100% vacc., R0 = 7, Prev inf = 67%
parsB_sept=parsA_sept; parsB_sept[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,parsA_sept$f_n*parsA_sept$f_pu+parsA_sept$f_v*parsA_sept$f_pv)
pars1B_sept=data.frame(f_b=seq(0,1,by=0.001),Rtb=Rt(cbind(parsB_sept, f_b=seq(0,1,by=0.001))),Scenario="US/ND-100:  R0=7.0, V: 100%,  I: 67%")

#Scenario C - 53% vacc., R0 = 7, Prev inf = 67%
parsC_sept=parsA_sept;parsC_sept["R0"]=7
pars1C_sept=data.frame(f_b=seq(0.00,round(max(case_vax_data_sept$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsC_sept, f_b=seq(0.00,round(max(case_vax_data_sept$frac_tot_vax), digits = 3),by=0.001))),Scenario="US/ND:         R0=7.0, V: 53%, I: 67%")

# #Scenario D - 59% vacc., R0 3.2, Prev inf = 60%
parsD_sept=data.frame(
  R0=3.2,
  f_n=1-max(case_vax_data_CA_sept$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_CA_sept$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_CA_sept$cum_unvax_pos[1]/(case_vax_data_CA_sept$cum_unvax_pos[1]+case_vax_data_CA_sept$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_CA_sept$cum_hybrid)/max(case_vax_data_CA_sept$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_CA_sept, #prev infected waned Inf
  VE_ii=V_i_inf_CA_sept, #prev infected waned trans
  VE_sv=mRNA_pfizer_ratio*(V_s_pfi_CA_sept)+mRNA_moderna_ratio*(V_s_mod_CA_sept), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=mRNA_pfizer_ratio*(V_i_pfi_CA_sept)+mRNA_moderna_ratio*(V_i_mod_CA_sept), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_CA_sept, #hybrid VEi
  VE_is=V_i_hybrid_CA_sept, #hybrid VEt
  VE_sb=V_s_boost_CA_sept, #boost inf
  VE_ib=V_i_boost_CA_sept #boost trans
)
pars1D_sept=data.frame(f_b=seq(0.00,round(max(case_vax_data_CA_sept$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsD_sept, f_b=seq(0.00,round(max(case_vax_data_CA_sept$frac_tot_vax), digits = 3),by=0.001))),Scenario="CA/D:           R0=3.2, V: 59%, I: 60%")

parsD_sept$f_n*parsD_sept$f_pu+parsD_sept$f_v*parsD_sept$f_pv


# #Scenario E - 24% vacc., R0 3.2, Prev inf = 0.2%
parsE_sept=data.frame(
  R0=7.0,
  f_n=1-max(case_vax_data_NZ_sept$frac_tot_vax), #fraction not vaccinated
  f_v=max(case_vax_data_NZ_sept$frac_tot_vax), #fraction vaccinated
  f_pu=case_vax_data_NZ_sept$cum_unvax_pos[1]/(case_vax_data_NZ_sept$cum_unvax_pos[1]+case_vax_data_NZ_sept$cum_unvax_neg[1]), #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=max(case_vax_data_NZ_sept$cum_hybrid)/max(case_vax_data_NZ_sept$cum_vax), #Used cdc data on RR unvaccinated v. vaccinated
  VE_si=V_s_inf_NZ_sept, #prev infected waned Inf
  VE_ii=V_i_inf_NZ_sept, #prev infected waned trans
  VE_sv=V_s_pfi_NZ_sept, #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_iv=V_i_pfi_NZ_sept, #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_ss=V_s_hybrid_NZ_sept, #hybrid VEi
  VE_is=V_i_hybrid_NZ_sept, #hybrid VEt
  VE_sb=V_s_boost_NZ_sept, #boost inf
  VE_ib=V_i_boost_NZ_sept #boost trans
)
pars1E_sept=data.frame(f_b=seq(0.00,round(max(case_vax_data_NZ_sept$frac_tot_vax), digits = 3),by=0.001),Rtb=Rt(cbind(parsE_sept, f_b=seq(0.00,round(max(case_vax_data_NZ_sept$frac_tot_vax), digits = 3),by=0.001))),Scenario="NZ/ND:         R0=7.0, V: 24%, I: 0.2%")

# #Generate CIs for each scenario fitted Rt line
nd1=10000
draws=data.frame() #generate 10,000 draws from each estimate on logit scale, then transform
draws_sept=cbind("VE_si"=unlist(CI_susc_inf_sept),
                 "VE_ii"=unlist(CI_infect_inf_sept),
                 "VE_sv"=unlist(CI_susc_mRNA_sept),
                 "VE_iv"=unlist(CI_infect_mRNA_sept),
                 "VE_ss"=unlist(CI_susc_hybrid_sept),
                 "VE_is"=unlist(CI_infect_hybrid_sept),
                 "VE_sb"=unlist(CI_susc_boost_mRNA_sept),
                 "VE_ib"=unlist(CI_infect_boost_mRNA_sept))

draws_CA_sept=cbind("VE_si"=unlist(CI_susc_inf_CA_sept),
                    "VE_ii"=unlist(CI_infect_inf_CA_sept),
                    "VE_sv"=unlist(CI_susc_mRNA_CA_sept),
                    "VE_iv"=unlist(CI_infect_mRNA_CA_sept),
                    "VE_ss"=unlist(CI_susc_hybrid_CA_sept),
                    "VE_is"=unlist(CI_infect_hybrid_CA_sept),
                    "VE_sb"=unlist(CI_susc_boost_mRNA_CA_sept),
                    "VE_ib"=unlist(CI_infect_boost_mRNA_CA_sept))

draws_NZ_sept=cbind("VE_si"=unlist(CI_susc_inf_NZ_sept),
                    "VE_ii"=unlist(CI_infect_inf_NZ_sept),
                    "VE_sv"=unlist(CI_susc_pfizer_NZ_sept),
                    "VE_iv"=unlist(CI_infect_pfizer_NZ_sept),
                    "VE_ss"=unlist(CI_susc_hybrid_NZ_sept),
                    "VE_is"=unlist(CI_infect_hybrid_NZ_sept),
                    "VE_sb"=unlist(CI_susc_boost_pfizer_NZ_sept),
                    "VE_ib"=unlist(CI_infect_boost_pfizer_NZ_sept))

#Draws for point with vaccines used to double vax instead of boosting
sigma_inf=diag(NA,nrow=2,ncol=2)
sigma_trans=diag(NA,nrow=2,ncol=2)

mean_inf=c(cvalues$c0[cvalues$endpoint=="Susceptibility"],cvalues$c1[cvalues$endpoint=="Susceptibility"])
sigma_inf[1,1]=cvalues$sigma11[cvalues$endpoint=="Susceptibility"];sigma_inf[1,2]=cvalues$sigma12[cvalues$endpoint=="Susceptibility"];
sigma_inf[2,1]=cvalues$sigma21[cvalues$endpoint=="Susceptibility"];sigma_inf[2,2]=cvalues$sigma22[cvalues$endpoint=="Susceptibility"]
x_susc=rmvnorm(n=nds, mean=mean_inf, sigma=sigma_inf)

mean_trans=c(cvalues$c0[cvalues$endpoint=="Infectiousness"],cvalues$c1[cvalues$endpoint=="Infectiousness"])
sigma_trans[1,1]=cvalues$sigma11[cvalues$endpoint=="Infectiousness"];sigma_trans[1,2]=cvalues$sigma12[cvalues$endpoint=="Infectiousness"];
sigma_trans[2,1]=cvalues$sigma21[cvalues$endpoint=="Infectiousness"];sigma_trans[2,2]=cvalues$sigma22[cvalues$endpoint=="Infectiousness"]
x_infect=rmvnorm(n=nds, mean=mean_trans, sigma=sigma_trans)


draws_sept_star=cbind("VE_si"=unlist(CI_susc_inf_sept),
                      "VE_ii"=unlist(CI_infect_inf_sept),
                      "VE_sv"=unlist(CI_susc_mRNA_sept),
                      "VE_iv"=unlist(CI_infect_mRNA_sept),
                      "VE_sv_star"=unlist(frac_pfizer*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATRpfi/NATRdelta))))))+frac_moderna*(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATRmod/NATRdelta))))))),
                      "VE_iv_star"=unlist(frac_pfizer*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATRpfi/NATRdelta))))))+frac_moderna*(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATRmod/NATRdelta))))))),
                      "VE_ss"=unlist(CI_susc_hybrid_sept),
                      "VE_is"=unlist(CI_infect_hybrid_sept),
                      "VE_ss_star"=unlist(1-(1/(1+exp(-(x_susc[,1]+x_susc[,2]*log2(NATRpfi*NATR_hybrid_ratio/NATRdelta)))))),
                      "VE_is_star"=unlist(1-(1/(1+exp(-(x_infect[,1]+x_infect[,2]*log2(NATRpfi*NATR_hybrid_ratio/NATRdelta)))))),
                      "VE_sb"=unlist(CI_susc_boost_mRNA_sept),
                      "VE_ib"=unlist(CI_infect_boost_mRNA_sept))

for (i in 1:1001) { #calculate CIs for each scenario
  if (i<245) {
    pars1E_sept[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsE_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_NZ_sept,f_b=pars1E_sept$f_b[i])))
  }
  if (i<534) { #2 scenarios that have max boosters = 60.4%
    pars1A_sept[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsA_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_sept,f_b=pars1A_sept$f_b[i])))
    pars1C_sept[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsC_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_sept,f_b=pars1C_sept$f_b[i])))
  }
  if(i<589){ #scenario that has max booster = 75%
    pars1D_sept[i,c("lower","upper")]=
      quantile(probs = c(0.16,0.84),x=Rt(cbind(parsD_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_CA_sept,f_b=pars1D_sept$f_b[i])))
  }
  pars1B_sept[i,c("lower","upper")]=
    quantile(probs = c(0.16,0.84),x=Rt(cbind(parsB_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_sept,f_b=pars1B_sept$f_b[i])))
}

#Add point for shots going to unvaxxed
Rt_star=function(parsR){
  with(parsR,
       R0*((f_n*(1-f_pu))+(f_n*f_pu)*(1-VE_si)*(1-VE_ii)+(1-f_pv)*(f_v-f_b)*(1-VE_sv)*(1-VE_iv)+(1-f_pv_star)*(f_v_star)*(1-VE_sv_star)*(1-VE_iv_star)+(f_pv)*(f_v-f_b)*(1-VE_ss)*(1-VE_is)+(f_pv_star)*(f_v_star)*(1-VE_ss_star)*(1-VE_is_star)+(f_b)*(1-VE_sb)*(1-VE_ib)))}


parsF_sept=parsA_sept
parsF_sept[,c("R0","f_n","f_v_star","f_pv_star","f_pu","VE_sv_star","VE_iv_star","VE_ss_star","VE_is_star")]=c(3.2,(1-0.5322761*1.5),0.2661381,1-0.6093848-(0.7478375*(0.202/0.4677239)),0.7478375*(0.202/0.4677239),1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_two*NATR_var_Delta))))),
                                                                                                               1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_two*NATR_var_Delta))))),
                                                                                                               1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Susceptibility"]+cvalues$c1[cvalues$endpoint=="Susceptibility"]*log2(NATR_pfizer_recent_two*NATR_hybrid_ratio*NATR_var_Delta))))),
                                                                                                               1-(1/(1+exp(-(cvalues$c0[cvalues$endpoint=="Infectiousness"]+cvalues$c1[cvalues$endpoint=="Infectiousness"]*log2(NATR_pfizer_recent_two*NATR_hybrid_ratio*NATR_var_Delta))))));parsF_sept$f_b=0
pars1F_sept=data.frame(f_b=0,Rtb=Rt_star(parsF_sept),Scenario="R0=3.2; f_vacc=79.8%, no third doses")
pars1F_sept[c("lower","upper")]=
  quantile(probs = c(0.16,0.84),x=Rt_star(cbind(parsF_sept[c("R0","f_n","f_v","f_pu","f_pv","f_pv_star","f_v_star")],draws_sept_star,f_b=pars1F_sept$f_b)))

# #group data into single dataframe
Rt_dat_sept=data.frame(rbind(pars1A_sept,pars1B_sept,pars1C_sept,pars1D_sept,pars1E_sept))#,pars1E,pars1F)) 

Rt_dat_sept$plot=pars1F_sept$plot="D: The impact of third doses on transmission"

Rt_dat_sept$Scenario[Rt_dat_sept$Scenario=="US/D:           R0=3.2, V: 53%, I: 67%"]="US/D"
Rt_dat_sept$Scenario[Rt_dat_sept$Scenario=="US/ND-100:  R0=7.0, V: 100%,  I: 67%"]="US/ND-100"
Rt_dat_sept$Scenario[Rt_dat_sept$Scenario=="US/ND:         R0=7.0, V: 53%, I: 67%"]="US/ND"
Rt_dat_sept$Scenario[Rt_dat_sept$Scenario=="CA/D:           R0=3.2, V: 59%, I: 60%"]="CA/D"
Rt_dat_sept$Scenario[Rt_dat_sept$Scenario=="NZ/ND:         R0=7.0, V: 24%, I: 0.2%"]="NZ/ND"
Rt_dat_sept$Scenario=factor(Rt_dat_sept$Scenario,levels=c("US/D","CA/D","NZ/ND","US/ND","US/ND-100"))

Rt_dat_sept2=Rt_dat_sept[Rt_dat_sept$Scenario=="US/D",]

#Fig 3D
Figure3D_plot=ggplot(data=Rt_dat_sept2,aes(x=f_b,y=Rtb),)+
  geom_line(colour= "#FF410DFF")+
  geom_ribbon(aes(ymin=lower,ymax=upper),fill="#FF410DFF",alpha=0.3)+
  geom_point(data=pars1F_sept,aes(x=f_b,y=Rtb),size=5, colour= "#FF410DFF")+
  geom_errorbar(data=pars1F_sept,aes(x=f_b,ymin=lower,ymax=upper),width=0, colour= "#FF410DFF")+
  theme_few()+scale_x_continuous(limits=c(0, 1.01), expand = c(.01, 0)) +
  scale_y_continuous(limits=c(0.8, 1.5),expand = c(0,0),labels = scales::number_format(accuracy = 0.01))+
  geom_hline(yintercept=1, linetype="longdash")+
  xlab("Fraction of population receiving 3rd dose")+
  ylab(expression("Effective Reproductive Number ("*R[t]*")"))+
  labs(fill="")+#,title="D: The impact of third doses on transmission")+
  scale_color_tron()+
  scale_fill_tron()+
  geom_text(data=pars1F_sept,aes(x=f_b,y=Rtb-0.07,label="R0 = 3.2; 79.8% vacc., no 3rd doses"),size=5,hjust=-0.05)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        #strip.text.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"),
        axis.title=element_text(size=15),legend.position = "right",
        axis.text=element_text(size=15),legend.text=element_text(size=15),
        legend.title=element_blank(),title = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=15))+
  facet_wrap(vars(plot),ncol=1)

#Load data from VE_over_time R script
Figure3_df=read.csv("Figure3ABC2.csv")
FracImmune=read.csv("FracImmune.csv")

#+/-SE instead of 95%
Figure3_df$lower2=Figure3_df$mean-(Figure3_df$mean-Figure3_df$lower)/1.96
Figure3_df$upper2=Figure3_df$mean+(Figure3_df$mean-Figure3_df$lower)/1.96

FracImmune$plot="A: Vaccination and Infection History"
FracImmune$Immunity=factor(FracImmune$Immunity,levels=c("Frac. Hybrid","Frac. Inf.","Frac. Vacc.","Frac. Boost" ))
Figure3_df$date=as.Date(Figure3_df$date);FracImmune$date=as.Date(FracImmune$date)
Figure3_df$VE[Figure3_df$VE=="Whole pop."]="Total"
Figure3_df$VE=factor(Figure3_df$VE,levels=c("Hybrid","Prev. Inf.","Vacc.","Total" ))

rcts2=read.csv("case_var_rect.csv");rcts2$xmin=as.Date(rcts2$xmin,format="%m/%d/%Y");rcts2$xmax=as.Date(rcts2$xmax,format="%m/%d/%Y");rcts2$date=as.Date(rcts2$date,format="%m/%d/%Y")
rcts2=rcts2[rcts2$location=="United States",];rcts2$ymin=1.;rcts2$ymax=1.1
rcts2$xmin[rcts2$Variant=="Wild Type & D614G"]=min(Figure3_df$date)
rcts2$xmax[rcts2$Variant=="Delta"]=max(Figure3_df$date)

Fig3A_plot=ggplot()+
  geom_line(data=FracImmune
            ,aes(x=date,y=mean,color=Immunity))+
  scale_color_tron() +
  guides(color=guide_legend("Vacc. and Inf. History"))+
  scale_x_date(expand=c(0,0),date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(expand = c(0,0),limits=c(0,0.6),labels = scales::number_format(accuracy = 0.01))+
  theme_few()+
  labs(y="Fraction of the population",x="")+
  geom_vline(xintercept=as.numeric(Figure3_df$date[Figure3_df$date=="2021-10-01"
                                                   |Figure3_df$date=="2021-11-01"
                                                   |Figure3_df$date=="2021-12-01" ])
             , linetype="dotted", colour="black")+
  geom_vline(xintercept=as.numeric(Figure3_df$date[Figure3_df$date=="2021-09-01"])
             , linetype="dotted", colour="black",linewidth=0.7)+
  geom_rect(data=rcts2[rcts2$Variant=="Wild Type & D614G",],mapping=aes(xmin=xmin,xmax=xmax,ymin=0.54,ymax=0.6),fill="#CCCCCC",col="black")+
  geom_rect(data=rcts2[rcts2$Variant=="Alpha",],mapping=aes(xmin=xmin,xmax=xmax,ymin=0.54,ymax=0.6),fill="#999999",col="black")+
  geom_rect(data=rcts2[rcts2$Variant=="Delta",],mapping=aes(xmin=xmin,xmax=xmax,ymin=0.54,ymax=0.6),fill="#666666",col="black")+
  geom_text(data=rcts2[rcts2$Variant=="Wild Type & D614G",],mapping=aes(x=as.Date("2021-02-14"),y=.573,label="WT & D614G"),size=4.5,color="white")+
  geom_text(data=rcts2[rcts2$Variant=="Alpha",],mapping=aes(x=as.Date("2021-05-15"),y=.573,label=Variant),size=4.5,color="white")+
  geom_text(data=rcts2[rcts2$Variant=="Delta",],mapping=aes(x=as.Date("2021-9-18"),y=.573,label=Variant),size=4.5,color="white")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.title=element_text(size=15),legend.position = "right",
        axis.text=element_text(size=15),legend.text=element_text(size=15),
        legend.title=element_text(size=15),title = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=15))+
  # annotate("text", x=Figure3_df$date[Figure3_df$date=="2021-04-21"],size=4, y=0.5, label= "Alpha 50%") + 
  # annotate("text", x=Figure3_df$date[Figure3_df$date=="2021-07-06"],size=4, y=0.5, label= "Delta 50%") + 
  facet_wrap(vars(plot),ncol=1)

Fig3BC_plot=ggplot()+
  geom_line(data=Figure3_df[Figure3_df$VE!="Frac. Inf."&Figure3_df$VE!="Frac. Vacc."&Figure3_df$VE!="Frac. Boost",]
            ,aes(x=date,y=mean,color=VE))+
  
  geom_ribbon(data=Figure3_df[Figure3_df$VE!="Frac. Inf."&Figure3_df$VE!="Frac. Vacc."&Figure3_df$VE!="Frac. Boost",]
              ,aes(x=date,ymin=lower2,ymax=upper2,fill=VE),alpha=0.5)+
  scale_fill_tron() +
  scale_color_tron() +
  scale_x_date(expand=c(0,0),date_labels="%b",date_breaks  ="1 month")+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  theme_few()+
  labs(y="Vaccine effectiveness",x="")+
  # geom_vline(xintercept=as.numeric(Figure3_df$date[Figure3_df$date=="2021-04-02"
  #                                                  |Figure3_df$date=="2021-06-18"])
  #            , linetype="solid", colour="black")+
  
  
  geom_vline(xintercept=as.numeric(Figure3_df$date[Figure3_df$date=="2021-10-01"
                                                   |Figure3_df$date=="2021-11-01"
                                                   |Figure3_df$date=="2021-12-01" ])
             , linetype="dotted", colour="black")+
  geom_vline(xintercept=as.numeric(Figure3_df$date[Figure3_df$date=="2021-09-01"])
             , linetype="dotted", colour="black",linewidth=0.7)+
  geom_rect(data=rcts2[rcts2$Variant=="Wild Type & D614G",],mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#CCCCCC",col="black")+
  geom_rect(data=rcts2[rcts2$Variant=="Alpha",],mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#999999",col="black")+
  geom_rect(data=rcts2[rcts2$Variant=="Delta",],mapping=aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="#666666",col="black")+
  geom_text(data=rcts2[rcts2$Variant=="Wild Type & D614G",],mapping=aes(x=as.Date("2021-02-14"),y=1.055,label="WT & D614G"),size=4.5,color="white")+
  geom_text(data=rcts2[rcts2$Variant=="Alpha",],mapping=aes(x=as.Date("2021-05-15"),y=1.055,label=Variant),size=4.5,color="white")+
  geom_text(data=rcts2[rcts2$Variant=="Delta",],mapping=aes(x=as.Date("2021-9-18"),y=1.055,label=Variant),size=4.5,color="white")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title=element_text(size=15),legend.position = "right",
        axis.text=element_text(size=15),legend.text=element_text(size=15),
        legend.title=element_text(size=15),title = element_text(size=15),
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(size=15))+
  # annotate("text", x=Figure3_df$date[Figure3_df$date=="2021-04-21"],size=4, y=0.99, label= "Alpha 50%") + 
  # annotate("text", x=Figure3_df$date[Figure3_df$date=="2021-07-06"],size=4, y=0.99, label= "Delta 50%") + 
  facet_wrap(vars(plot),ncol=1)

ggdraw(xlim=c(0,0.65),ylim=c(0,1.6)) +
  draw_plot(Fig3A_plot, x = 0, y = 1.09, width = 0.627, height = 0.4) +
  draw_plot(Fig3BC_plot, x = 0, y = .4, width = 0.593, height = 0.7) +
  draw_plot(Figure3D_plot, x = 0, y = 0, width = 0.532, height = 0.4) 

#Supplemental Plot - Modify data for legend/facets
pars1A_dec$Scenario2=pars1A_nov$Scenario2=pars1A_oct$Scenario2=pars1A_sept$Scenario2="US/D"
pars1B_dec$Scenario2=pars1B_nov$Scenario2=pars1B_oct$Scenario2=pars1B_sept$Scenario2="US/ND-100"
pars1C_dec$Scenario2=pars1C_nov$Scenario2=pars1C_oct$Scenario2=pars1C_sept$Scenario2="US/ND"
pars1D_dec$Scenario2=pars1D_nov$Scenario2=pars1D_oct$Scenario2=pars1D_sept$Scenario2="CA/D"
pars1E_dec$Scenario2=pars1E_nov$Scenario2=pars1E_oct$Scenario2=pars1E_sept$Scenario2="NZ/ND"

pars1A_sept$Month=pars1B_sept$Month=pars1C_sept$Month=pars1D_sept$Month=pars1E_sept$Month="A: September 1, 2021"
pars1A_oct$Month=pars1B_oct$Month=pars1C_oct$Month=pars1D_oct$Month=pars1E_oct$Month="B: October 1, 2021"
pars1A_nov$Month=pars1B_nov$Month=pars1C_nov$Month=pars1D_nov$Month=pars1E_nov$Month="C: November 1, 2021"
pars1A_dec$Month=pars1B_dec$Month=pars1C_dec$Month=pars1D_dec$Month=pars1E_dec$Month="D: December 1, 2021"

#Bind data
Rt_dat_tot=data.frame(rbind(pars1A_sept,pars1B_sept,pars1C_sept,pars1D_sept,pars1E_sept,
                            pars1A_oct,pars1B_oct,pars1C_oct,pars1D_oct,pars1E_oct,
                            pars1A_nov,pars1B_nov,pars1C_nov,pars1D_nov,pars1E_nov,
                            pars1A_dec,pars1B_dec,pars1C_dec,pars1D_dec,pars1E_dec))

#relevel
Rt_dat_tot$Scenario2=factor(Rt_dat_tot$Scenario2,levels=c("US/D","CA/D","NZ/ND","US/ND","US/ND-100"))

#Supplemental Figure
ggplot(data=Rt_dat_tot,aes(x=f_b,y=Rtb),)+
  geom_line(aes(col=Scenario2))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Scenario2),alpha=0.15)+
  theme_few()+scale_x_continuous(limits=c(0, 1.1), expand = c(.0, 0.0),breaks=c(0,0.25,0.5,0.75,1))+
  geom_hline(yintercept=1, linetype="longdash")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Fraction of population receiving 3rd dose")+
  ylab(expression("Effective Reproductive Number ("*R[t]*")"))+
  labs(color="Scenario",fill="Scenario")+scale_fill_manual(values=c("red","yellow3","blue" ,"purple","green"))+
  scale_color_manual(values=c("red","yellow3","blue" ,"purple","green"))+
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position="right" )+
  theme(axis.title=element_text(size=25),#legend.position = c(.1, .80),
        axis.text=element_text(size=15),legend.text=element_text(size=20),
        legend.title=element_text(size=20),panel.background = element_rect(fill = NA, color = "black"))+
  facet_wrap(~Month,nrow = 2)

## Start deaths and infections averted calculation
m=Rt_dat_tot
minmax= m %>% group_by(Scenario,Scenario2, Month) %>% summarize(
  minRt=min(Rtb),maxRt=max(Rtb),ratio=minRt/maxRt,maxB=max(f_b)) 
arrange(minmax,Scenario2,Month)

#CIs for change in Rt calculations - September
rt_ratCI_sept=data.frame(month="Sept",lower=NA,upper=NA,
                         rt_rat=pars1A_sept$Rtb[pars1A_sept$f_b==max(pars1A_sept$f_b)]/pars1A_sept$Rtb,f_b=pars1A_sept$f_b,f_bmax=max(pars1A_sept$f_b))

for (i in 1:nrow(rt_ratCI_sept)) {
  rt_ratCI_sept[i,c("lower","upper")]=quantile(probs = c(0.025,0.975),
                                               x=Rt(cbind(parsA_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_sept,f_b=max(pars1A_sept$f_b)))
                                               /Rt(cbind(parsA_sept[c("R0","f_n","f_v","f_pu","f_pv")],draws_sept,f_b=(i-1)/1000)))
}

#Bind into single dataframe
rt_ratCI=rbind(rt_ratCI_sept,rt_ratCI_oct,rt_ratCI_nov,rt_ratCI_dec)

#Calculation of deaths/infections averted by boosting
#Load Rt/deaths/cases/fraction boosted data
m4=read.csv("Rt_deaths_infections2.csv")
m4$date=as.Date(m4$date)
mindate=as.Date("2021-06-01");maxdate=as.Date("2021-12-31")

#rolling mean
m4$daily_deaths_sm=rollmean(m4$daily_deaths,k=7,fill=NA,align="right")
m4$daily_cases_sm=rollmean(m4$daily_cases,k=7,fill=NA,align="right")
m4$daily_cases_smD=m4$daily_cases_sm*m4$FracDelta
m4=m4[complete.cases(m4),] #remove edge rows w/ NAs from rolling means

#Estimate infections from DEATHS by deconvolution #Distributions from Lewnard et al 2020 medRxiv
infec_deaths=deconvolve_incidence(incidence_data = m4$daily_deaths_sm,
                                  delay=list(inf_shed,shed_sympt,sympt_hosp,hosp_death))
#Daily total infections
m4$daily_inf_est_deaths=
  1*infec_deaths$values[c((-infec_deaths$index_offset+1):length(infec_deaths$values),
                          rep(NA,-infec_deaths$index_offset))]
#Calculate fraction of infections due to delta - using FracDelta from sequences shifted by 11 days (infection to reporting date)
m4$daily_inf_est_deathsD=NA
m4$daily_inf_est_deathsD[1:(nrow(m4)-11)]=m4$daily_inf_est_deaths[1:(nrow(m4)-11)]*
  m4$FracDelta[12:nrow(m4)]

#Reconvolution distribution
nd=1000000 #1M draws
inftodeathV = rweibull(nd,shape=5.983,scale=1.455)+ #delay infection to shedding
  rweibull(nd,shape=0.294,scale= 0.14)+ #shedding to symptom onset
  rgamma(nd,shape=5.078,scale=0.765)+ #symptom onset to hospitalization
  rgamma(nd,shape=sh1,scale=sc1) #fitted distribution hospitalization to death
#count reconvolution deaths by day after rounding day of death
inftodeathH=hist(round(inftodeathV),probability=T,breaks=c(0:100-.5,max(inftodeathV)+.5) )
inftodeathDist=inftodeathH$counts/nd #turn into probability vector length = 100 days

#Calculate reconvoluted deaths
m4$deaths_reconvD=0
for(i in 1:(nrow(m4)-101)) { #loop through days of infection
  #Add deaths occurring on each day after day of infections to appropriate day
  m4$deaths_reconvD[(i+1):(i+101)]=m4$daily_inf_est_deathsD[i]*inftodeathDist+
    m4$deaths_reconvD[(i+1):(i+101)]
}

#Calculate reduced infections
m3=m[m$Scenario2=="US/D",]

#Create column with multiplier for # of infections based on daily decrease due to boosters
BoostDate_sept=as.Date("2021-09-01") #change this for 4 different dates sept, oct, nov, dec, or make loop
BoostDate_oct=as.Date("2021-10-01") 
BoostDate_nov=as.Date("2021-11-01") 
BoostDate_dec=as.Date("2021-12-01") 
m4$mult_dec=1; m4$mult_decL=1; m4$mult_decU=1; 

#add column for observed boosting
m4$f_b=round(m4$Frac_boosted,digits=3)
rt_ratCI$f_b=round(rt_ratCI$f_b,digits=3)

#Calculate Rt deduction and CIs on each day based on observed boosting
rt_ratCI_sept=merge(m4[m4$date>=BoostDate_sept,c("date","f_b")],rt_ratCI[rt_ratCI$month=="Sept",c("lower","upper","rt_rat","f_b")],by=c("f_b"),all.x=T)
rt_ratCI_sept=rt_ratCI_sept %>% arrange(ymd(rt_ratCI_sept$date))

rt_ratCI_oct=merge(m4[m4$date>=BoostDate_oct,c("date","f_b")],rt_ratCI[rt_ratCI$month=="Oct",c("lower","upper","rt_rat","f_b")],all.x=T)
rt_ratCI_oct=rt_ratCI_oct %>% arrange(ymd(rt_ratCI_oct$date))

rt_ratCI_nov=merge(m4[m4$date>=BoostDate_nov,c("date","f_b")],rt_ratCI[rt_ratCI$month=="Nov",c("lower","upper","rt_rat","f_b")],all.x=T)
rt_ratCI_nov=rt_ratCI_nov %>% arrange(ymd(rt_ratCI_nov$date))

rt_ratCI_dec=merge(m4[m4$date>=BoostDate_dec,c("date","f_b")],rt_ratCI[rt_ratCI$month=="Dec",c("lower","upper","rt_rat","f_b")],all.x=T)
rt_ratCI_dec=rt_ratCI_dec %>% arrange(ymd(rt_ratCI_dec$date))

#Calculate change in Rt w/ CIs each day 
m4$drRt_sept[m4$date>=BoostDate_sept]=(rt_ratCI_sept$rt_rat)^(1/3.8)
m4$drRt_oct[m4$date>=BoostDate_oct]=rt_ratCI_oct$rt_rat^(1/3.8)
m4$drRt_nov[m4$date>=BoostDate_nov]=rt_ratCI_nov$rt_rat^(1/3.8)
m4$drRt_dec[m4$date>=BoostDate_dec]=rt_ratCI_dec$rt_rat^(1/3.8)

m4$drRt_septL[m4$date>=BoostDate_sept]=(rt_ratCI_sept$lower)^(1/3.8)
m4$drRt_octL[m4$date>=BoostDate_oct]=rt_ratCI_oct$lower^(1/3.8)
m4$drRt_novL[m4$date>=BoostDate_nov]=rt_ratCI_nov$lower^(1/3.8)
m4$drRt_decL[m4$date>=BoostDate_dec]=rt_ratCI_dec$lower^(1/3.8)

m4$drRt_septU[m4$date>=BoostDate_sept]=rt_ratCI_sept$upper^(1/3.8)
m4$drRt_octU[m4$date>=BoostDate_oct]=rt_ratCI_oct$upper^(1/3.8)
m4$drRt_novU[m4$date>=BoostDate_nov]=rt_ratCI_nov$upper^(1/3.8)
m4$drRt_decU[m4$date>=BoostDate_dec]=rt_ratCI_dec$upper^(1/3.8)

#Create column with multiplier for # of infections based on daily decrease due to boosters
m4$mult_sept=1;m4$mult_septL=1;m4$mult_septU=1; m4$mult_oct=1;m4$mult_octL=1;m4$mult_octU=1;
m4$mult_nov=1;m4$mult_novL=1;m4$mult_novU=1;m4$mult_dec=1;m4$mult_decL=1;m4$mult_decU=1;

#Take product of daily Rt reductions for calculation of reduced infections and deaths each day over time
#September
m_sept=m4[m4$date>=BoostDate_sept,]
for (i in 1:nrow(m_sept)) {
  m_sept$mult_sept[i]=prod(m_sept$drRt_sept[1:i])
  m_sept$mult_septL[i]=prod(m_sept$drRt_septL[1:i])
  m_sept$mult_septU[i]=prod(m_sept$drRt_septU[1:i])
}

m4$mult_sept[m4$date>=BoostDate_sept]=m_sept$mult_sept
m4$mult_septL[m4$date>=BoostDate_sept]=m_sept$mult_septL
m4$mult_septU[m4$date>=BoostDate_sept]=m_sept$mult_septU

#October
m_oct=m4[m4$date>=BoostDate_oct,]
for (i in 1:nrow(m_oct)) {
  m_oct$mult_oct[i]=prod(m_oct$drRt_oct[1:i])
  m_oct$mult_octL[i]=prod(m_oct$drRt_octL[1:i])
  m_oct$mult_octU[i]=prod(m_oct$drRt_octU[1:i])
}

m4$mult_oct[m4$date>=BoostDate_oct]=m_oct$mult_oct
m4$mult_octL[m4$date>=BoostDate_oct]=m_oct$mult_octL
m4$mult_octU[m4$date>=BoostDate_oct]=m_oct$mult_octU

#November
m_nov=m4[m4$date>=BoostDate_nov,]
for (i in 1:nrow(m_nov)) {
  m_nov$mult_nov[i]=prod(m_nov$drRt_nov[1:i])
  m_nov$mult_novL[i]=prod(m_nov$drRt_novL[1:i])
  m_nov$mult_novU[i]=prod(m_nov$drRt_novU[1:i])
}

m4$mult_nov[m4$date>=BoostDate_nov]=m_nov$mult_nov
m4$mult_novL[m4$date>=BoostDate_nov]=m_nov$mult_novL
m4$mult_novU[m4$date>=BoostDate_nov]=m_nov$mult_novU

#December
m_dec=m4[m4$date>=BoostDate_dec,]
for (i in 1:nrow(m_dec)) {
  m_dec$mult_dec[i]=prod(m_dec$drRt_dec[1:i])
  m_dec$mult_decL[i]=prod(m_dec$drRt_decL[1:i])
  m_dec$mult_decU[i]=prod(m_dec$drRt_decU[1:i])
}

m4$mult_dec[m4$date>=BoostDate_dec]=m_dec$mult_dec
m4$mult_decL[m4$date>=BoostDate_dec]=m_dec$mult_decL
m4$mult_decU[m4$date>=BoostDate_dec]=m_dec$mult_decU

#Calculate relative infections averted in each boosting scenario
m4$daily_inf_est_deathsD_reduced_sept=m4$daily_inf_est_deathsD*m4$mult_sept
m4$daily_inf_est_deathsD_reduced_oct=m4$daily_inf_est_deathsD*m4$mult_oct
m4$daily_inf_est_deathsD_reduced_nov=m4$daily_inf_est_deathsD*m4$mult_nov
m4$daily_inf_est_deathsD_reduced_dec=m4$daily_inf_est_deathsD*m4$mult_dec

m4$daily_inf_est_deathsD_reduced_septL=m4$daily_inf_est_deathsD*m4$mult_septL
m4$daily_inf_est_deathsD_reduced_octL=m4$daily_inf_est_deathsD*m4$mult_octL
m4$daily_inf_est_deathsD_reduced_novL=m4$daily_inf_est_deathsD*m4$mult_novL
m4$daily_inf_est_deathsD_reduced_decL=m4$daily_inf_est_deathsD*m4$mult_decL

m4$daily_inf_est_deathsD_reduced_septU=m4$daily_inf_est_deathsD*m4$mult_septU
m4$daily_inf_est_deathsD_reduced_octU=m4$daily_inf_est_deathsD*m4$mult_octU
m4$daily_inf_est_deathsD_reduced_novU=m4$daily_inf_est_deathsD*m4$mult_novU
m4$daily_inf_est_deathsD_reduced_decU=m4$daily_inf_est_deathsD*m4$mult_decU

#calculate reconvoluted deaths with reduced infections
m4$deaths_reconvD_reduced_sept=0;m4$deaths_reconvD_reduced_oct=0;
m4$deaths_reconvD_reduced_nov=0;m4$deaths_reconvD_reduced_dec=0;

m4$deaths_reconvD_reduced_septL=0;m4$deaths_reconvD_reduced_octL=0;
m4$deaths_reconvD_reduced_novL=0;m4$deaths_reconvD_reduced_decL=0;

m4$deaths_reconvD_reduced_septU=0;m4$deaths_reconvD_reduced_octU=0;
m4$deaths_reconvD_reduced_novU=0;m4$deaths_reconvD_reduced_decU=0;
for(i in 1:(nrow(m4)-101)) { #loop through days of infection
  #Add deaths occurring on each day after day of infections to appropriate day
  m4$deaths_reconvD_reduced_sept[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_sept[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_sept[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_oct[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_oct[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_oct[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_nov[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_nov[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_nov[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_dec[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_dec[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_dec[(i+1):(i+101)]
  
  m4$deaths_reconvD_reduced_septL[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_septL[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_septL[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_octL[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_octL[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_octL[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_novL[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_novL[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_novL[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_decL[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_decL[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_decL[(i+1):(i+101)]
  
  m4$deaths_reconvD_reduced_septU[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_septU[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_septU[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_octU[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_octU[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_octU[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_novU[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_novU[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_novU[(i+1):(i+101)]
  m4$deaths_reconvD_reduced_decU[(i+1):(i+101)]=m4$daily_inf_est_deathsD_reduced_decU[i]*inftodeathDist+
    m4$deaths_reconvD_reduced_decU[(i+1):(i+101)]
}


#Bind subsetted dataframes and prep for plotting
m4_septD=data.frame(subset(m4,select=c(date,deaths_reconvD,
                                       deaths_reconvD_reduced_sept
)),month="September") %>% 
  rename(deaths_reconvD_reduced=deaths_reconvD_reduced_sept)

m4_octD=data.frame(subset(m4,select=c(date,deaths_reconvD,
                                      deaths_reconvD_reduced_oct)),month="October") %>% 
  rename(deaths_reconvD_reduced=deaths_reconvD_reduced_oct)

m4_novD=data.frame(subset(m4,select=c(date,deaths_reconvD,
                                      deaths_reconvD_reduced_nov)),month="November") %>% 
  rename(deaths_reconvD_reduced=deaths_reconvD_reduced_nov)

m4_decD=data.frame(subset(m4,select=c(date,deaths_reconvD,
                                      deaths_reconvD_reduced_dec)),month="December") %>% 
  rename(deaths_reconvD_reduced=deaths_reconvD_reduced_dec)

death_FP=rbind(m4_septD,m4_octD,m4_novD,m4_decD)##death four panel figure dataframe
death_FP$plotcol="Deaths"

death_FPL=pivot_longer(death_FP
                       ,cols=c(deaths_reconvD,deaths_reconvD_reduced),
                       names_to="metric",values_to="value")

death_FP$month=factor(death_FP$month,
                      levels = c("September","October","November","December"))

death_FPL$month=factor(death_FPL$month,
                       levels = c("September","October","November","December"))
death_FPL$metric=factor(death_FPL$metric,
                        levels = c("deaths_reconvD","deaths_reconvD_reduced"))

m4_septI=data.frame(subset(m4,select=c(date,
                                       daily_inf_est_deaths,
                                       daily_inf_est_deathsD,
                                       daily_inf_est_deathsD_reduced_sept)),month="September") %>% 
  rename(daily_inf_est_deathsD_reduced=daily_inf_est_deathsD_reduced_sept)

m4_octI=data.frame(subset(m4,select=c(date,
                                      daily_inf_est_deaths,
                                      daily_inf_est_deathsD,
                                      daily_inf_est_deathsD_reduced_oct)),month="October") %>% 
  rename(daily_inf_est_deathsD_reduced=daily_inf_est_deathsD_reduced_oct)

m4_novI=data.frame(subset(m4,select=c(date,
                                      daily_inf_est_deaths,
                                      daily_inf_est_deathsD,
                                      daily_inf_est_deathsD_reduced_nov)),month="November") %>% 
  rename(daily_inf_est_deathsD_reduced=daily_inf_est_deathsD_reduced_nov)

m4_decI=data.frame(subset(m4,select=c(date,
                                      daily_inf_est_deaths,
                                      daily_inf_est_deathsD,
                                      daily_inf_est_deathsD_reduced_dec)),month="December") %>% 
  rename(daily_inf_est_deathsD_reduced=daily_inf_est_deathsD_reduced_dec)

inf_FP=rbind(m4_septI,m4_octI,m4_novI,m4_decI)##death four panel figure dataframe
inf_FP$plotcol="Infections"

inf_FPL=pivot_longer(inf_FP
                     ,cols=c(daily_inf_est_deaths,
                             daily_inf_est_deathsD,
                             daily_inf_est_deathsD_reduced),
                     names_to="metric",values_to="value")

inf_FP$month=factor(inf_FP$month,
                    levels = c("September","October","November","December"))

inf_FPL$month=factor(inf_FPL$month,
                     levels = c("September","October","November","December"))
inf_FPL$metric=factor(inf_FPL$metric,
                      levels = c("daily_inf_est_deaths","daily_inf_est_deathsD",
                                 "daily_inf_est_deathsD_reduced"))

death_inf_FPL=rbind(death_FPL,inf_FPL)

death_inf_FPL$month=factor(death_inf_FPL$month,
                           levels = c("September","October","November","December"))
death_inf_FPL$metric=factor(death_inf_FPL$metric,
                            levels = c("daily_inf_est_deaths","daily_inf_est_deathsD",
                                       "daily_inf_est_deathsD_reduced","deaths_reconvD",
                                       "deaths_reconvD_reduced"))



aesname_IDFP=c("All Infections","Delta Infections","Delta Infections reduced",
               "Delta Deaths","Delta Deaths reduced")
names(aesname_IDFP)=c("daily_inf_est_deaths","daily_inf_est_deathsD",
                      "daily_inf_est_deathsD_reduced","deaths_reconvD",
                      "deaths_reconvD_reduced")

aesname_IFP=c("All Infections","Delta Infections","Delta Infections reduced")
names(aesname_IFP)=c("daily_inf_est_deaths","daily_inf_est_deathsD",
                     "daily_inf_est_deathsD_reduced")

aesname_DFP=c("Delta Deaths","Delta Deaths reduced")
names(aesname_DFP)=c("deaths_reconvD",
                     "deaths_reconvD_reduced")

mindate_FP=as.Date("2021-08-27");maxdate_FP=as.Date("2022-01-30")

#infections/ infections averted plot
inf_FP_plot=ggplot(death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP,])+theme_few()+
  geom_line(data=death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP
                               &death_inf_FPL$metric=="daily_inf_est_deathsD_reduced",],aes(x=date,y=value/0.003,col=metric))+
  geom_line(data=death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP
                               &death_inf_FPL$metric=="daily_inf_est_deathsD",],aes(x=date,y=value/0.003,col=metric))+
  geom_line(data=death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP
                               &death_inf_FPL$metric=="daily_inf_est_deaths",],aes(x=date,y=value/0.003,col=metric))+
  geom_ribbon(data=inf_FP[inf_FP$date>mindate_FP&inf_FP$date<maxdate_FP,],
              aes(x=date,ymin=daily_inf_est_deathsD_reduced/0.003,ymax=daily_inf_est_deathsD/0.003),fill="#FF9933",alpha=0.5)+
  
  scale_y_continuous(
    labels = scales::label_comma(),
  )+
  scale_x_date(expand=c(0,0),date_labels="%b-%d",date_breaks  ="1 month")+
  labs(y="Daily infections",x="",color="")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=16),text = element_text(size=18),legend.position = "bottom",legend.text = element_text(size=15))+
  scale_color_manual(values=c("black","#8B1A1A","#FF3336"),labels=aesname_IFP,breaks=c("daily_inf_est_deaths","daily_inf_est_deathsD",
                                                                                       "daily_inf_est_deathsD_reduced"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~month, nrow=4,scales="free_y")

#deaths/ deaths averted plot
death_FP_plot=ggplot(death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP,])+theme_few()+
  geom_line(data=death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP
                               &death_inf_FPL$metric=="deaths_reconvD_reduced",],aes(x=date,y=value,col=metric))+
  geom_line(data=death_inf_FPL[death_inf_FPL$date>mindate_FP&death_inf_FPL$date<maxdate_FP
                               &death_inf_FPL$metric=="deaths_reconvD",],aes(x=date,y=value,col=metric))+
  geom_ribbon(data=death_FP[death_FP$date>mindate_FP&death_FP$date<maxdate_FP,],
              aes(x=date,ymin=deaths_reconvD_reduced,ymax=deaths_reconvD),fill="#3A8BA5",alpha=0.5)+
  scale_y_continuous(
    labels = scales::label_comma(),
  )+
  scale_x_date(expand=c(0,0),date_labels="%b-%d",date_breaks  ="1 month")+
  labs(y="Daily deaths",x="",color="")+
  theme(axis.title=element_text(size=20),axis.text=element_text(size=16),text = element_text(size=18),legend.position = "bottom",legend.text = element_text(size=15))+
  scale_color_manual(values=c("#00008B","#98F5FF"),labels=aesname_DFP,breaks=c("deaths_reconvD","deaths_reconvD_reduced"))+
  guides(colour = guide_legend(nrow = 1))+
  facet_wrap(~month, nrow=4,scales="free_y")

#Modify plots for final figure
inf_FP_plot2=tag_facet(inf_FP_plot, open = "(", close = ")", tag_pool = c("A","B","C","D"), x = as.Date("2021-09-01"),y = Inf)
death_FP_plot2=tag_facet(death_FP_plot, open = "(", close = ")", tag_pool = c("E","F","G","H"), x = as.Date("2021-09-01"),y = Inf)

inf_FP_plot3=inf_FP_plot2+theme(strip.text = element_text())
death_FP_plot3=death_FP_plot2+theme(strip.text = element_text())

#Figure 4
ggarrange(inf_FP_plot3,death_FP_plot3,nrow=1)

#Deaths and infections averted calculations (using IFR = 0.3%)
m4L=pivot_longer(subset(m4,select=-c(deaths_reconvD_reduced_oct,deaths_reconvD_reduced_nov,
                                     deaths_reconvD_reduced_dec,daily_inf_est_deathsD_reduced_oct,
                                     daily_inf_est_deathsD_reduced_nov,daily_inf_est_deathsD_reduced_dec))
                 ,cols=c(deaths_reconvD_reduced_sept,deaths_reconvD,
                         daily_inf_est_deathsD_reduced_sept,daily_inf_est_deathsD,
                         daily_inf_est_deaths),
                 names_to="metric",values_to="value")

fw_plot=data.frame(metric=c("daily_inf_est_deaths","daily_inf_est_deathsD",
                            "daily_inf_est_deathsD_reduced_sept","deaths_reconvD",
                            "deaths_reconvD_reduced_sept"),plot=c("Infections","Infections","Infections","Deaths","Deaths"))
m4L=merge(m4L,fw_plot)

deaths_D=m4L[m4L$date>mindate&m4L$date<maxdate&m4L$metric=="deaths_reconvD",]
deaths_D_reduced=m4L[m4L$date>mindate&m4L$date<maxdate&m4L$metric=="deaths_reconvD_reduced_sept",]
deaths_D=deaths_D[order(deaths_D$date),]
deaths_D_reduced=deaths_D_reduced[order(deaths_D_reduced$date),]

sum(deaths_D$value)-sum(deaths_D_reduced$value)
sum(m4L$value[m4L$date>mindate&m4L$date<maxdate&m4L$metric=="daily_inf_est_deathsD"])/0.003-
  sum(m4L$value[m4L$date>mindate&m4L$date<maxdate&m4L$metric=="daily_inf_est_deathsD_reduced"])/0.003

deaths_averted=data.frame(month=c("Sept","Oct","Nov","Dec"),lower=NA,mean=NA,upper=NA)
deaths_averted$lower[deaths_averted$month=="Sept"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_septU[m4$date>mindate])
deaths_averted$mean[deaths_averted$month=="Sept"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_sept[m4$date>mindate])
deaths_averted$upper[deaths_averted$month=="Sept"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_septL[m4$date>mindate])

deaths_averted$lower[deaths_averted$month=="Oct"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_octU[m4$date>mindate])
deaths_averted$mean[deaths_averted$month=="Oct"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_oct[m4$date>mindate])
deaths_averted$upper[deaths_averted$month=="Oct"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_octL[m4$date>mindate])

deaths_averted$lower[deaths_averted$month=="Nov"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_novU[m4$date>mindate])
deaths_averted$mean[deaths_averted$month=="Nov"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_nov[m4$date>mindate])
deaths_averted$upper[deaths_averted$month=="Nov"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_novL[m4$date>mindate])

deaths_averted$lower[deaths_averted$month=="Dec"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_decU[m4$date>mindate])
deaths_averted$mean[deaths_averted$month=="Dec"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_dec[m4$date>mindate])
deaths_averted$upper[deaths_averted$month=="Dec"]=sum(m4$deaths_reconvD[m4$date>mindate])-sum(m4$deaths_reconvD_reduced_decL[m4$date>mindate])


deaths_averted[,c("inf_lower","inf_mean","inf_upper")]=deaths_averted[,c("lower","mean","upper")]/0.003

#Deaths averted per million vaccinated
deaths_averted[deaths_averted$month=="Sept",c("deaths_per_mil_vax_lower","deaths_per_mil_vax_mean","deaths_per_mil_vax_upper")]=deaths_averted[deaths_averted$month=="Sept",c("lower","mean","upper")]*1000000/(0.5322761*332858873)
deaths_averted[deaths_averted$month=="Oct",c("deaths_per_mil_vax_lower","deaths_per_mil_vax_mean","deaths_per_mil_vax_upper")]=deaths_averted[deaths_averted$month=="Oct",c("lower","mean","upper")]*1000000/(0.5610537*332858873)
deaths_averted[deaths_averted$month=="Nov",c("deaths_per_mil_vax_lower","deaths_per_mil_vax_mean","deaths_per_mil_vax_upper")]=deaths_averted[deaths_averted$month=="Nov",c("lower","mean","upper")]*1000000/(0.5788164*332858873)
deaths_averted[deaths_averted$month=="Dec",c("deaths_per_mil_vax_lower","deaths_per_mil_vax_mean","deaths_per_mil_vax_upper")]=deaths_averted[deaths_averted$month=="Dec",c("lower","mean","upper")]*1000000/(0.5939408*332858873)

#Fig S1 Plot infections, cases, deaths

#Fig S5 Daily vaccinated plot
daily_vax=case_vax_data #already read in above
daily_vax$date=seq(as.Date('2021-12-13'), as.Date('2020-1-23'), by = "-1 days")
daily_vax$vax_avg=rollmean(daily_vax$Vaccinated,7,fill=NA)

ggplot()+
  geom_line(data=daily_vax[1:343,],aes(x=date,y=vax_avg/1000))+
  ylab("People fully vaccinated (1000s) against COVID-19 (7 day rolling average)")+  
  scale_x_date(expand=c(0,0),date_labels = "%b %y",date_breaks = "1 month")+
  scale_y_continuous(expand=c(0,0),limits =c(0,NA))+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title=element_text(size=15),axis.text=element_text(size=20),
        legend.text=element_text(size=20),legend.title=element_text(size=20))

#Supplemental Figure 1
#Subset data for rt/cases
Rt_death=m4[449:651,]

#Add legend variable
rt=data.frame(Rt_death[c("date","rt","rt_lower","rt_upper","daily_cases_sm","daily_cases_smD")],data="Rt")
cases_delta=data.frame(Rt_death[c("date","rt","rt_lower","rt_upper","daily_cases_sm","daily_cases_smD")],data="Delta Cases")
cases_all=data.frame(Rt_death[c("date","rt","rt_lower","rt_upper","daily_cases_sm","daily_cases_smD")],data="All Cases")

#Bind data and relevel
Rt_death2=rbind(rt,cases_delta,cases_all)
Rt_death2$data=factor(Rt_death2$data,levels=c("All Cases","Delta Cases","Rt"))

##Figure S1 - Warning for cases going over ylim after Omicron invades
ggplot(Rt_death2)+
  geom_hline(yintercept=1, linetype="dashed", color="black")+
  geom_line(data=Rt_death2,aes(x=date,y=rt,color="Rt"))+
  geom_ribbon(data=Rt_death2,aes(x=date,ymin=rt_lower,ymax=rt_upper),fill="royalblue",alpha=0.4)+
  geom_line(data=Rt_death2,aes(x=date,y=daily_cases_sm/250000,color="All Cases"),alpha=0.75,size=1.)+
  geom_line(data=Rt_death2,aes(x=date,y=daily_cases_smD/250000,color="Delta Cases"),alpha=0.75,size=1.)+
  scale_x_date(name= "",expand = c(0, 0),date_labels = "%b-%d",date_breaks = "1 month")+
  scale_y_continuous(
    expand = c(0, 0), limits=c(0,2.2),
    # Features of the first axis
    name = expression("Reproductive Number ("*R[t]*")"),
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*250, name="Cases (1000s; 7-day rolling average)",breaks = c(0,100,200,300,400,500,600,700,800))
  )+
  scale_color_manual(values=c(  "forestgreen","green2","royalblue"))+
  theme_bw()+
  geom_vline(xintercept=as.numeric(Rt_death2$date[Rt_death2$date=="2021-10-01"
                                                  |Rt_death2$date=="2021-11-01"
                                                  |Rt_death2$date=="2021-12-01" ])
             , linetype="dotted", colour="black")+
  geom_vline(xintercept=as.numeric(Rt_death2$date[Rt_death2$date=="2021-09-01"])
             , linetype="dotted", colour="black",linewidth=0.7)+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position="bottom")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title=element_text(size=18),axis.text=element_text(size=25),
        legend.text=element_text(size=25),legend.title=element_blank())

