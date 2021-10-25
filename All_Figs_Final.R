lapply(c("ggplot2","ggthemes","tidyverse","dplyr","lme4"),require,character.only=T) #load packages
setwd("c:/marm/research/covid-19/booster/")
mAll=read.csv("glm_VE_total2.csv") #load VE, nAbs data
m1=read.csv("Nabs ratios Delta WT.csv") #load VE, nAbs data

#Analysis of nAbs ratios relative to Delta
h1=lmer(log(Ratio)~Vaccine+(1|Study),data=m1); h2=lmer(log(Ratio)~1+(1|Study),data=m1);AIC(h1,h2);anova(h2,h1)
Delta_WT_nAbs_Ratio=exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]);Delta_WT_nAbs_Ratio
#Adjust nAbs ratios for Delta using the mixed effects mean from studies for AZ, Pfizer, conv. sera
mAll$Nabs_Ratio[mAll$Delta==1]= mAll$Nabs_Ratio[mAll$Delta==1]/Delta_WT_nAbs_Ratio

#Only show error bars for vaccine averages, so make new columns w/ NAs for individual studies
mAll$Lower_avg=NA;mAll$Upper_avg=NA;mAll$Lower_avg[mAll$Vacc_avg==1]=mAll$Lower[mAll$Vacc_avg==1];mAll$Upper_avg[mAll$Vacc_avg==1]=mAll$Upper[mAll$Vacc_avg==1]

#Create separate data frame for transmission data; #remove transmission data from main data frame
mAllt=mAll[mAll$Groups=="Delta - transmission",]; mAll=mAll[mAll$Groups!="Delta - transmission",]

#waning ratios from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
wane_ratio=1/(.71*497/83+.29*538/41)
boost_ratio=.71*1754/497+.29*1318/538
boost_ratio/wane_ratio

#Table S3: Fitted model for protection vs nAbs
f2=glm(VE~log2(Nabs_Ratio)*Groups,data=mAll[mAll$Vacc_avg==1,],family=binomial,weights=N_eff);summary(f2)

#Calculating predicted VE values all infections Pfizer w/ waning, boosting & adding to mAll
p_w=predict(f2, newdata=data.frame(Nabs_Ratio=wane_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Pfizer",Study="Waned",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=wane_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],VE=plogis(unlist(p_w[1])),Lower_avg=plogis(unlist(p_w[1])-1.96*unlist(p_w[2])),
                     Upper_avg=plogis(unlist(p_w[1])+1.96*unlist(p_w[2])) )

p_b=predict(f2, newdata=data.frame(Nabs_Ratio=boost_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Pfizer",Study="Boosted",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=boost_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],VE=plogis(unlist(p_b[1])),Lower_avg=plogis(unlist(p_b[1])-1.96*unlist(p_b[2])),
                     Upper_avg=plogis(unlist(p_b[1])+1.96*unlist(p_b[2])) )

#labels for plot for Pfizer orig, boosting, waning points
tlab=data.frame(Delta=1,text1=c("Pfizer orig","            Pfizer 8 mo waning","Pfizer boosted"),
                Nabs_Ratio=c(1,wane_ratio,boost_ratio)*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],
                VE=mAll$VE[mAll$Vaccine=="Pfizer"&mAll$Study%in%c("Waned","Pouwels (delta)","Boosted")]-.04)

#Create fitted lines
xval=seq(0.1,5,by=0.01)
mL=data.frame(Nabs_Ratio=rep(xval,4),Groups=rep(c("Delta - infection","Delta - symptomatic","non-Delta - symptomatic","non-Delta - infection"),each=length(xval)),
              Variant_Grouping=rep(c("Delta - All infections","Delta - Symptomatic","non-Delta - Symptomatic","non-Delta - infection"),each=length(xval)),
              Delta=c(rep(1,2*length(xval)),rep(0,2*length(xval))  ) )
mL$VE_pred=predict(f2,newdata=mL, type="response")

#Create CIs for fig 1
mL$VE_logit=unlist(predict(f2,newdata=mL, type="link", se.fit = TRUE)[1])
mL$VE_se=unlist(predict(f2,newdata=mL, type="link", se.fit = TRUE)[2])
mL$upper_bound = plogis(mL$VE_logit + (1.96 * mL$VE_se))
mL$lower_bound = plogis(mL$VE_logit- (1.96 * mL$VE_se))

#Fig 1
#Reorder factors for plotting
mAll$Groups=factor(mAll$Groups,levels = c("non-Delta - symptomatic","non-Delta - infection","Delta - symptomatic","Delta - infection"))
mAll$Variant_Grouping=factor(mAll$Variant_Grouping,levels=c("non-Delta - Study estim.","non-Delta - Vaccine estim.",
  "Delta - Study estim.","Delta - Vaccine estim.","Delta - Prediction"))
tlab2=data.frame(Delta=c(0,1),text1=c("A: non-Delta","B: Delta  "),Nabs_Ratio=c(0.1,0.1),VE=c(0.98,0.98))
mAll$Nabs_Ratio_jit=mAll$Nabs_Ratio+rnorm(nrow(mAll),0,0.02) #jitter xvals
mAll$Nabs_Ratio_jit[mAll$Variant_Grouping=="Delta - Prediction"]=mAll$Nabs_Ratio[mAll$Variant_Grouping=="Delta - Prediction"] #unjitter predicted xvals
ggplot()+
  geom_point(data=mAll,aes(x=Nabs_Ratio_jit,y=VE,col=Groups,shape=Variant_Grouping),size=4)+
  geom_errorbar(data=mAll,aes(x=Nabs_Ratio_jit,ymin=Lower_avg, ymax=Upper_avg,col=Groups), width=0)+
  scale_x_continuous(trans='log2',limits=c(0.1,5),n.breaks=7)+theme_few()+
  scale_shape_manual(values = c(21,19,2,17,0))+#  ylim(0.35,1)+ 
  scale_color_manual(values=c("blue","green","orange","red"))+
  xlab("Ratio of antibody neutralization relative to convalescent")+
  ylab("Protection against infection or disease")+
  theme(axis.title=element_text(size=25),#legend.position = c(.1, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  geom_line(data=mL[mL$Groups=="Delta - symptomatic",],aes(x=Nabs_Ratio,y=VE_pred),col="orange")+
  geom_line(data=mL[mL$Groups=="Delta - infection",],aes(x=Nabs_Ratio,y=VE_pred),col="red")+
  geom_line(data=mL[mL$Groups=="non-Delta - symptomatic",],aes(x=Nabs_Ratio,y=VE_pred),col="blue")+
  geom_line(data=mL[mL$Groups=="non-Delta - infection",],aes(x=Nabs_Ratio,y=VE_pred),col="green")+
  geom_ribbon(data=mL[mL$Groups=="Delta - symptomatic",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="orange",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="Delta - infection",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="red",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="non-Delta - symptomatic",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="blue",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="non-Delta - infection",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="green",alpha=0.2)+
  geom_text(data=tlab,aes(x=Nabs_Ratio,y=VE,label=text1))+
  geom_text(data=tlab2,aes(x=Nabs_Ratio,y=VE,label=text1),size=6,hjust=0)+
  labs(col="Variant-Endpoint", shape="Variant-Grouping")+
  facet_wrap(.~Delta,nrow=2, scales="free_y")+theme(strip.text.x = element_blank())

#Transmission model
f3=glm(VE~log2(Nabs_Ratio),data=mAllt,family=binomial,weights=N_eff);summary(f3)

# Generate predicted min protection vs transmission given infection for boosted nAbs
p_tb=predict(f3,newdata=data.frame(Nabs_Ratio=boost_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Pfizer",Study="Boost_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=boost_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],
                       VE=plogis(unlist(p_tb[1])))
mAllt[mAllt$Study=="Boost_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tb[1])+c(-1,1)*1.96*unlist(p_tb[2]))
p_tw=predict(f3,newdata=data.frame(Nabs_Ratio=wane_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Pfizer",Study="Wane_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=wane_ratio*mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],
                       VE=plogis(unlist(p_tw[1])))
mAllt[mAllt$Study=="Wane_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tw[1])+c(-1,1)*1.96*unlist(p_tw[2]))

#Create fitted lines & CIs for VE for transmission
mLt=data.frame(Nabs_Ratio=xval,Groups="Delta - transmission")
mLt$VE_logit=unlist(predict(f3,newdata=mLt, type="link", se.fit = TRUE)[1])
mLt$VE_se=unlist(predict(f3,newdata=mLt, type="link", se.fit = TRUE)[2])
mLt$fit_line = plogis(mLt$VE_logit)
mLt$upper_bound = plogis(mLt$VE_logit + (1.96 * mLt$VE_se))
mLt$lower_bound = plogis(mLt$VE_logit- (1.96 * mLt$VE_se))

#Fig S4 Plot VE for transmission
ggplot()+
  geom_point(data=mAllt,aes(x=Nabs_Ratio,y=VE,shape=Variant_Grouping),size=4,col="red")+
  scale_shape_manual(values = c(0,17))+
  geom_errorbar(data=mAllt,aes(x=Nabs_Ratio,ymin=Lower_avg, ymax=Upper_avg,col=Groups), width=0)+
  geom_ribbon(data = mLt,aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="red",alpha=0.2)+
  geom_line(data = mLt,aes(x=xval,y=fit_line),color="red")+
  scale_x_continuous(trans='log2',limits=c(min(xval),max(xval)),n.breaks=7)+
  theme_few()+ylim(0,1)+
  xlab("Ratio of antibody neutralization relative to convalescent")+
  ylab("Minimum protection against transmission given infection")+
  theme(axis.title=element_text(size=15),#legend.position = c(.1, .80),
        axis.text=element_text(size=15),
        legend.title=element_blank(),
        legend.position = "none")

### Calculate VE for acquired immunity from prev infection using daily US deaths 
#Fit model for waning rate of nAbs over time
nAbs_inf=data.frame(time=c(0,122,244,365),nabs=c(3.7,1.6,1.1,1.2)) #Havervall data
f1=nls(nabs~c0*exp(c1*time)+c2,start=list(c0=2.1,c1=-.01,c2=1),data=nAbs_inf);summary(f1)
deaths=read.csv("Death.csv")
deaths$pnabs=(predict(f1,newdata=deaths))
deaths$prot_ratio=deaths$pnabs/max(deaths$pnabs)
deaths$VE=predict(f2,newdata=data.frame(Nabs_Ratio=deaths$prot_ratio/exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1])
                                        ,Groups="Delta - infection"),type="response")

#Function time v nABs Pfizer
#Fig S1 Plot pfizer waning over time
nAbs_pfizer=data.frame(time=c(0,0,19.2,24,47.9,74.9,103.9,131.8,158,236),nabs=c(1,1,0.726,0.689,0.371,0.206,0.201,0.167,0.165,0.139),data_from="Falsey") #Falsey Nabs data
f4=nls(nabs~c0*exp(c1*time)+c2,start=list(c0=0.90965429,c1=-0.01383033,c2=0.10804145),data=nAbs_pfizer);summary(f4)

nAbs_pfizer_plot=data.frame(time=0:330,nabs=(predict(f4,newdata=data.frame(time=0:330))),data_from="Predict")
nAbs_pfizer_plot=rbind(nAbs_pfizer_plot,nAbs_pfizer)
nAbs_pfizer_plot$nabs=nAbs_pfizer_plot$nabs/max(nAbs_pfizer$nabs)

ggplot()+
  geom_point(data=nAbs_pfizer_plot[nAbs_pfizer_plot$data_from=="Falsey",],aes(x=time,y=nabs),shape=19)+
  geom_line(data=nAbs_pfizer_plot[nAbs_pfizer_plot$data_from=="Predict",],aes(x=time,y=nabs))+
  xlab("Days since 7 days post vaccination")+
  ylab("Pfizer BioNTech neutralizing antibodies relative to peak")+
  ylim(0,1.05)+
  scale_x_continuous(breaks = seq(0,330,by=30))+
  theme_bw()+
  theme(strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85) )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Moderna waning over time
#nAbs_moderna_prelog=data.frame(time=c(0,76,166),nabs=c(1,0.544,0.236))#Doria-Rose
nAbs_moderna=data.frame(time=c(0,76,166),nabs=c(log(1),log(0.544),log(0.236))) #Doria-Rose
f5=nls(nabs~c0*exp(c1*time)+c2,start=list(c0=4,c1=0.005,c2=-4),
       data=nAbs_moderna,trace = F,control = list(warnOnly=T));summary(f5)
#Model is perfect fit so doesn't meet convergence criteria

#Hybrid immunity/Super nAbs
daily_super=read.csv("Frac_prev_inf_breakdown.csv")
daily_super$pNabs=(predict(f1,newdata=data.frame(time=0:(278))))
daily_super$prot_ratio=daily_super$pNabs/max(daily_super$pNabs)
daily_super$nABs_pfizer=mAll$Nabs_Ratio[mAll$Study=="Pouwels (delta)"&mAll$Vaccine=="Pfizer"]*
  predict(f4,newdata=data.frame(time=daily_super$time_vax)) 
daily_super$nABs_moderna=mAll$Nabs_Ratio[mAll$Study=="Nasreen(delta)"&mAll$Vaccine=="Moderna"]*
  exp(predict(f5,newdata=data.frame(time=daily_super$time_vax)) )
super_pfizer_ratio=6300/800
daily_super$nABs_super=daily_super$nABs_pfizer*super_pfizer_ratio

#CIs for vaccine
#Pfizer and Moderna Waning
nds=10000
qdraws=runif(nds)
y1=plogis(matrix(qnorm(rep(qdraws,each=nrow(daily_super)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_pfizer,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )

qdraws0=runif(nds)
y2=plogis(matrix(qnorm(rep(qdraws0,each=nrow(daily_super)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_moderna,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_moderna,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )
x1=daily_super$Prop_vax_pfizer
x2=daily_super$Prop_vax_moderna
x3=x1%*%y1
x4=x2%*%y2
CI_inf_vax1y=(1.503/2.503)*x3+(1/2.503)*x4
CI_inf_vaxy=as.data.frame(t((1.503/2.503)*x3+(1/2.503)*x4))
#quantile(CI_inf_vax1y, probs = c(0.025,.50,0.975))

qdraws1=runif(nds)
y5=plogis(matrix(qnorm(rep(qdraws1,each=nrow(daily_super)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_pfizer,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )

qdraws2=runif(nds)
y6=plogis(matrix(qnorm(rep(qdraws2,each=nrow(daily_super)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_moderna,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_moderna,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )

x5=x1%*%y5
x6=x2%*%y6
CI_trans_vax1y=y7=(1.503/2.503)*x5+(1/2.503)*x6
CI_trans_vaxy=as.data.frame(t((1.503/2.503)*x5+(1/2.503)*x6))
#quantile(CI_trans_vax1y, probs = c(0.025,.50,0.975))

#Super Waning
qdraws3=runif(nds)
y8=plogis(matrix(qnorm(rep(qdraws3,each=nrow(daily_super)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_super,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_super$nABs_super,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )

z3=daily_super$Prop_super 
CI_inf_super1=z3%*% y8
CI_inf_super=as.data.frame(t(z3%*% y8))
#quantile(CI_inf_super1, probs = c(0.025,.50,0.975))

qdraws4=runif(nds)
y9=plogis(matrix(qnorm(rep(qdraws4,each=nrow(daily_super)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_super,
                                                                                       Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_super$nABs_super,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_super)) )

CI_trans_super1=z3%*% y9
CI_trans_super=as.data.frame(t(z3%*% y9))
#quantile(CI_trans_super1, probs = c(0.025,.50,0.975))

#Waning natural infection
nds=10000
qdraws5=runif(nds)
z0=plogis(matrix(qnorm(rep(qdraws5,each=nrow(deaths)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=deaths$prot_ratio/exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]),
                                                                                                Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=deaths$prot_ratio/exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]),
                                                               
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2]) ),nrow=nrow(deaths)) )

z1=deaths$Percent_inf 
CI_inf1y=z1%*% z0
CI_infy=as.data.frame(t(z1%*% z0))
#quantile(CI_inf1y, probs = c(0.025,.50,0.975))

qdraws6=runif(nds)
z2=plogis(matrix(qnorm(rep(qdraws,each=nrow(deaths)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=deaths$prot_ratio/exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]),
                                                                                               Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=deaths$prot_ratio/exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]),
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2]) ),nrow=nrow(deaths)) )

CI_trans1y=z1%*% z2
CI_transy=as.data.frame(t(z1%*% z2))
#quantile(CI_trans1y, probs = c(0.025,0.50,0.975))

#Fig S3 Plot natural infection nAbs over time
nAbs_inf_plot=data.frame(time=c(0,122,244,365),nabs=c(3.7,1.6,1.1,1.2),data_from="Havervall") 
nAbs_plot=data.frame(time=0:(20*24),nabs=(predict(f1,newdata=data.frame(time=0:(20*24)))),data_from="Predict")
nAbs_plot=rbind(nAbs_inf_plot,nAbs_plot)
nAbs_plot$nabs=nAbs_plot$nabs/max(nAbs_inf_plot$nabs)

ggplot()+
  geom_point(data=nAbs_plot[nAbs_plot$data_from=="Havervall",],aes(x=time,y=nabs),shape=19)+
  geom_line(data=nAbs_plot[nAbs_plot$data_from=="Predict",],aes(x=time,y=nabs))+
  xlab("Days post infection")+
  ylab("Neutralizing antibodies relative to peak")+
  ylim(0,1.05)+
  scale_x_continuous(breaks = seq(0,480,by=30))+
  theme_bw()+
  theme(strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85) )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Fig S2 Plot Moderna waning over time
nAbs_moderna_prelog=data.frame(time=c(0,76,166),nabs=c(1,0.544,0.236),data_from="Doria-Rose")
nAbs_moderna=data.frame(time=c(0:330),nabs=mAll$Nabs_Ratio[mAll$Study=="Nasreen(delta)"&mAll$Vaccine=="Moderna"]*
                          exp(predict(f5,newdata=data.frame(time=c(0:330))))
                          ,data_from="Predict") #Doria-Rose
nAbs_moderna$nabs=nAbs_moderna$nabs/1.4343102
nAbs_moderna_plot=rbind(nAbs_moderna,nAbs_moderna_prelog)

ggplot()+
  geom_point(data=nAbs_moderna_plot[nAbs_moderna_plot$data_from=="Doria-Rose",],aes(x=time,y=nabs),shape=19)+
  geom_line(data=nAbs_moderna_plot[nAbs_moderna_plot$data_from=="Predict",],aes(x=time,y=nabs))+
  xlab("Days since 14 days post 2nd dose")+
  ylab("Moderna neutralizing antibodies relative to peak")+
  ylim(0,1.05)+
  scale_x_continuous(breaks = seq(0,330,by=30))+
  theme_bw()+
  theme(strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85) )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Rt function 
Rt=function(parsR){
  with(parsR,
       R0*((f_n*(1-f_pu))+(f_n*f_pu)*(1-VE_ii)*(1-VE_ti)+(1-f_pv)*(f_v-f_b)*(1-VE_iv)*(1-VE_tv)+(f_pv)*(f_v-f_b)*(1-VE_is)*(1-VE_ts)+(f_b)*(1-VE_ib)*(1-VE_tb)))}

#Scenario A - 56% vacc., R0 = 3.7
parsA=data.frame(
  R0=3.7,
  f_n=0.44, #fraction not vaccinated
  f_v=0.56, #fraction vaccinated
  f_pu=0.686, #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=0.470, #Used cdc data on RR unvaccinated v. vaccinated
  VE_ii=(quantile(CI_inf1y,probs=c(0.5))), #prev infected waned Inf
  VE_ti=(quantile(CI_trans1y,probs=c(0.5))), #prev infected waned trans
  VE_iv=(quantile(CI_inf_vax1y,probs=c(0.5))), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_tv=(quantile(CI_trans_vax1y,probs=c(0.5))), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_is=quantile(CI_inf_super1, probs = c(.5)), #Super VEi
  VE_ts=quantile(CI_trans_super1, probs = c(.5)), #Super VEt
  VE_ib=mAll$VE[mAll$Study=="Boosted"], #boost inf
  VE_tb=mAllt$VE[mAllt$Study=="Boost_Trans"] #boost trans
)
pars1A=data.frame(f_b=seq(0,0.56,by=0.01),Rtb=Rt(cbind(parsA, f_b=seq(0.00,0.56,by=0.01))),Scenario="R0=3.7;  56% vacc., 56.4% prev. inf.")

#Scenario B - 100% vacc., R0 = 7
parsB=parsA; parsB[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,0.564)
pars1B=data.frame(f_b=seq(0,1,by=0.01),Rtb=Rt(cbind(parsB, f_b=seq(0,1,by=0.01))),Scenario="R0=7;  100% vacc., 56.4% prev. inf.")

#Scenario C - 56% vacc., R0 = 7
parsC=parsA;parsC["R0"]=7
pars1C=data.frame(f_b=seq(0,.56,by=0.01),Rtb=Rt(cbind(parsC, f_b=seq(0,.56,by=0.01))),Scenario="R0=7;  56% vacc., 56.4% prev. inf.")

#Scenario D - 75% vacc., R0 3.7, Prev inf = 28%
parsE=parsA;parsE[,c("R0","f_n","f_v","f_pu","f_pv")]=c(3.7,0.25,0.75,0.421,0.238)
pars1E=data.frame(f_b=seq(0,.75,by=0.01),Rtb=Rt(cbind(parsE, f_b=seq(0,.75,by=0.01))),Scenario="R0=3.7;  75% vacc., 28.2% prev. inf.")

#Scenario E - 75% vacc., R0 3.7, Prev inf = 28%
parsF=parsA; parsF[,c("R0","f_n","f_v","f_pu","f_pv")]=c(3.7,0.4,0.6,0.01,0.002)
pars1F=data.frame(f_b=seq(0,0.6,by=0.01),Rtb=Rt(cbind(parsF, f_b=seq(0,.6,by=0.01))),Scenario="R0=3.7;  60% vacc., <1% prev. inf.")

#Calculating CIs for Rt for Fig 2
nd1=10000
draws=data.frame() #generate 10,000 draws from each estimate on logit scale, then transform
draws=cbind("VE_ii"=unlist(CI_infy),
            "VE_ti"=unlist(CI_transy),
            "VE_iv"=unlist(CI_inf_vaxy),
            "VE_tv"=unlist(CI_trans_vaxy),
            "VE_is"=unlist(CI_inf_super),
            "VE_ts"=unlist(CI_trans_super),
            "VE_ib"=(plogis(rnorm(nd1,unlist(p_b[1]),unlist(p_b[2])))),
            "VE_tb"=(plogis(rnorm(nd1,unlist(p_tb[1]),unlist(p_tb[2])))))


for (i in 1:101) { #calculate CIs for each scenario
  if (i<58) { #2 scenarios that have max boosters = 56%
    pars1A[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsA[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1A$f_b[i])))
    pars1C[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsC[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1C$f_b[i])))
  }
  if(i<62){ #Scenario similar to NZ
    pars1F[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsF[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1F$f_b[i])))
    
  }
  if(i<77){ #scenario that has max booster = 75%
    pars1E[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsE[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1E$f_b[i])))
  }
  pars1B[i,c("lower","upper")]=
    quantile(probs = c(0.025,0.975),x=Rt(cbind(parsB[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1B$f_b[i])))
}

#Add point for shots going to unvaxxed
parsD=parsC
parsD[,c("R0","f_n","f_v","f_pv")]=c(7,0.16,0.84,0.662);parsD$f_b=0
pars1D=data.frame(f_b=0,Rtb=Rt(parsD),Scenario="R0=7; f_vacc=84%,")
pars1D[c("lower","upper")]=
  quantile(probs = c(0.025,0.975),x=Rt(cbind(parsD[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1D$f_b)))

Rt_dat=data.frame(rbind(pars1A,pars1B,pars1C,pars1E,pars1F)) #group data into single dataframe

#Fig 2
ggplot(data=Rt_dat,aes(x=f_b,y=Rtb),)+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Scenario),alpha=0.15)+
  geom_point(data=pars1D,aes(x=f_b,y=Rtb),size=5)+
  theme_few()+
  geom_hline(yintercept=1, linetype="longdash")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Fraction of population receiving 3rd dose")+
  ylab("Effective Reproductive Number (Rt)")+
  labs(fill="Scenario")+scale_fill_manual(values=c("red", "grey","blue", "yellow3","green"))+
  scale_color_manual(values=c("red", "grey","blue","yellow3" ,"green"))+
  theme(strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85) )+
  geom_text(data=pars1D,aes(x=f_b,y=Rtb,label="R0 = 7; 84% vacc., no 3rd doses"),size=7,hjust=-0.05)+
  theme(axis.title=element_text(size=25),#legend.position = c(.1, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=20),
        legend.title=element_text(size=20))

min(Rt_dat$f_b[Rt_dat$Scenario=="R0=3.7;  75% vacc., 28.2% prev. inf."&Rt_dat$Rtb<1]) # % boosting needed for Rt<1
min(Rt_dat$f_b[Rt_dat$Scenario=="R0=7;  100% vacc., 56.4% prev. inf."&Rt_dat$Rtb<1])# % boosting needed for Rt<1

#Rt values at beginning & ends of lines
max(pars1A$Rtb);min(pars1A$Rtb);min(pars1A$Rtb)/max(pars1A$Rtb) #R0=3.7;  56% vacc., 56.4% prev. inf.
max(pars1B$Rtb);min(pars1B$Rtb);min(pars1B$Rtb)/max(pars1B$Rtb) #R0=7;  100% vacc., 56.4% prev. inf.
max(pars1C$Rtb);min(pars1C$Rtb);min(pars1C$Rtb)/max(pars1C$Rtb) #R0=7;  56% vacc., 56.4% prev. inf.
max(pars1E$Rtb);min(pars1E$Rtb);min(pars1E$Rtb)/max(pars1E$Rtb) #R0=3.7;  75% vacc., 28.2% prev. inf.
max(pars1F$Rtb);min(pars1F$Rtb);min(pars1F$Rtb)/max(pars1F$Rtb) #R0=3.7;  60% vacc., <1% prev. inf.
pars1D

#Fig S6 Plot infections, cases, deaths
#Crude IFR estimated by dividing total 4.2*total cases through 24 days prior by total deaths
#Infections estimated for day n-24 by dividing deaths on day n by IFR (0.003993)
daily_cases=read.csv("daily_cases.csv")
daily_cases$date=seq(as.Date('2020-1-23'), as.Date('2021-10-15'), by = "1 days")

ggplot()+
  geom_line(data=daily_cases,aes(x=date,y=New.Cases/1000,color="Cases"))+
  geom_line(data=daily_cases,aes(x=date,y=Infections/1000,color="Infections"))+
  geom_line(data=daily_cases,aes(x=date,y=Deaths,color="Deaths"))+
  #xlab("Date")+
  ylab("COVID-19 cases and infections in thousands and COVID-19 deaths")+
  scale_x_date(date_labels = "%b %Y",date_breaks = "2 month")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank() )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(color = "")


#Fig S5 Daily vaccinated plot
daily_fully_vax=read.csv("daily_vaccines.csv")
daily_fully_vax$date=seq(as.Date('2020-12-14'), as.Date('2021-10-22'), by = "1 days")

ggplot()+
  geom_line(data=daily_fully_vax,aes(x=date,y=Fully_Vaccinated))+
  #xlab("Date")+
  ylab("Number of people fully vaccinated against COVID-19")+  
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 month")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

#Fig S7 Rate ratio plot
daily_rate_ratio=read.csv("Rate_ratio.csv")
daily_rate_ratio$date=seq(as.Date('2021-04-04'), as.Date('2021-8-29'), by = "7 days")

ggplot(data=daily_rate_ratio,aes(x=date,y=age_adj_irr))+
  geom_line()+geom_point()+
  ylab("Ratio of COVID-19 cases in unvaccinated v. vaccinated")+
  scale_y_continuous(breaks =seq(0,14,2),limits=c(0,14))+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 month")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank() )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
