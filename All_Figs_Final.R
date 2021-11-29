lapply(c("ggplot2","ggthemes","tidyverse","dplyr","lme4","brms","rstan","mgcv"),require,character.only=T) #load packages
setwd("c:/marm/research/covid-19/booster/")
mAll=read.csv("VE_nAbs_data.csv") #load VE, nAbs data
m1=read.csv("Nabs ratios Delta WT.csv") #load VE, nAbs data

#Table S2 Analysis of nAbs ratios relative to Delta
h1=lmer(log(Ratio)~Vaccine+(1|Study),data=m1); h2=lmer(log(Ratio)~1+(1|Study),data=m1);AIC(h1,h2);anova(h2,h1)
Delta_WT_nAbs_Ratio=exp(summary(lmer(log(Ratio)~1+(1|Study),data=m1))$coefficients[1]);Delta_WT_nAbs_Ratio
#Adjust nAbs ratios for Delta using the mixed effects mean from studies for AZ, Pfizer, conv. sera
mAll$Nabs_Ratio[mAll$Delta==1]= mAll$Nabs_Ratio[mAll$Delta==1]/Delta_WT_nAbs_Ratio

#Only show error bars for vaccine averages, so make new columns w/ NAs for individual studies
mAll$Lower_avg=NA;mAll$Upper_avg=NA;mAll$Lower_avg[mAll$Vacc_avg==1]=mAll$Lower[mAll$Vacc_avg==1];mAll$Upper_avg[mAll$Vacc_avg==1]=mAll$Upper[mAll$Vacc_avg==1]

#Create separate data frame for transmission data; #remove transmission data from main data frame
mAllt=mAll[mAll$Groups=="Delta - transmission",]; mAll=mAll[mAll$Groups!="Delta - transmission",]

#Pfizer waning ratios from Falsey et al 2021 Fig 1A WT, age weighted (71%, 29%)
wane_ratio=1/(.71*497/83+.29*538/41)
boost_ratio=.71*1754/497+.29*1318/538
boost_ratio/wane_ratio
hybrid_ratio=6300/800 #Goel et al
#Moderna boosting @28d relative to post 2nd dose from Chu et al medRxiv 0,181, 209: 1268, 150.2, 1951.7
Mod_boost_ratio=1951.7/1268
Mod_wane_ratio=150.2/1268

Nabs_Pfiz_delta_inf=mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"]

#Table S3: Fitted model for protection vs nAbs
f2=glm(VE~log2(Nabs_Ratio)*Groups,data=mAll[mAll$Vacc_avg==1,],family=binomial,weights=N_eff);summary(f2)

#Calculating predicted VE values all infections Pfizer, Moderna w/ waning, boosting & adding to mAll
p_w=predict(f2, newdata=data.frame(Nabs_Ratio=wane_ratio*Nabs_Pfiz_delta_inf,Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Pfizer",Study="Waned",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=wane_ratio*Nabs_Pfiz_delta_inf,VE=plogis(unlist(p_w[1])),Lower_avg=plogis(unlist(p_w[1])-1.96*unlist(p_w[2])),
                     Upper_avg=plogis(unlist(p_w[1])+1.96*unlist(p_w[2])) )

p_b=predict(f2, newdata=data.frame(Nabs_Ratio=boost_ratio*Nabs_Pfiz_delta_inf,Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Pfizer",Study="Boosted",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=boost_ratio*Nabs_Pfiz_delta_inf,VE=plogis(unlist(p_b[1])),Lower_avg=plogis(unlist(p_b[1])-1.96*unlist(p_b[2])),
                     Upper_avg=plogis(unlist(p_b[1])+1.96*unlist(p_b[2])) )
Nabs_Mod_delta=mean(mAll$Nabs_Ratio[mAll$Vaccine=="Moderna"&mAll$Delta==1])
p_mw=predict(f2, newdata=data.frame(Nabs_Ratio=Mod_wane_ratio*Nabs_Mod_delta,Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Moderna",Study="Waned",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=Mod_wane_ratio*Nabs_Mod_delta,VE=plogis(unlist(p_mw[1])),Lower_avg=plogis(unlist(p_mw[1])-1.96*unlist(p_mw[2])),
                     Upper_avg=plogis(unlist(p_mw[1])+1.96*unlist(p_mw[2])) )

p_mb=predict(f2, newdata=data.frame(Nabs_Ratio=Mod_boost_ratio*Nabs_Mod_delta,Groups="Delta - infection"), type="link",se.fit=T)
mAll=mAll%>% add_row(Vaccine="Moderna",Study="Boosted",Variant_Grouping="Delta - Prediction",Endpoint="All infections",Groups="Delta - infection",
                     Delta=1,Nabs_Ratio=Mod_boost_ratio*Nabs_Mod_delta,VE=plogis(unlist(p_mb[1])),Lower_avg=plogis(unlist(p_mb[1])-1.96*unlist(p_mb[2])),
                     Upper_avg=plogis(unlist(p_mb[1])+1.96*unlist(p_mb[2])) )

#Moderna original protection (not waned or boosted); not currently plotted
p_mo=predict(f2, newdata=data.frame(Nabs_Ratio=Nabs_Mod_delta,Groups="Delta - infection"), type="link",se.fit=T)
plogis(unlist(p_mo[1])+c(-1,0,1)*1.96*unlist(p_mo[2]))

mAll[mAll$Study%in%c("Boosted","Waned","Nasreen(delta)","Pouwels (delta"),]

#labels for plot for Pfizer, Moderna orig, boosting, waning points
tlab=data.frame(Delta=1,text1=c("Pfizer orig","              Pfizer waned","   Pfizer boosted","Moderna waned","Moderna boosted"),
                Nabs_Ratio=c(c(1,wane_ratio,boost_ratio)*Nabs_Pfiz_delta_inf,Mod_wane_ratio*Nabs_Mod_delta,Mod_boost_ratio*Nabs_Mod_delta),
                VE=mAll$VE[mAll$Vaccine%in%c("Pfizer","Moderna")&mAll$Study%in%c("Waned","Pouwels (delta)","Boosted")]
                -c(1,1,-1,1,1)*.04)

#Create fitted lines for Fig 1
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
#library(scales);hue_pal()(4) #color codes for ggplot defaults: red, green, blue, purple
ggplot()+
  geom_point(data=mAll,aes(x=Nabs_Ratio_jit,y=VE,col=Groups,shape=Variant_Grouping),size=4)+
  geom_errorbar(data=mAll,aes(x=Nabs_Ratio_jit,ymin=Lower_avg, ymax=Upper_avg,col=Groups), width=0)+
  scale_x_continuous(trans='log2',limits=c(0.1,5),n.breaks=7)+theme_few()+
  scale_shape_manual(values = c(21,19,2,17,0))+ 
  scale_color_manual(values=c(  "#00BFC4","#7CAE00", "#C77CFF","#F8766D"))+ #colors: b,g,p,r
  xlab("Ratio of antibody neutralization relative to convalescent")+
  ylab("Protection against infection or disease")+
  theme(axis.title=element_text(size=15),#legend.position = "top",#legend.position = c(.1, .80),
        axis.text=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  geom_line(data=mL[mL$Groups=="Delta - symptomatic",],aes(x=Nabs_Ratio,y=VE_pred),col="#C77CFF")+
  geom_line(data=mL[mL$Groups=="Delta - infection",],aes(x=Nabs_Ratio,y=VE_pred),col="#F8766D")+
  geom_line(data=mL[mL$Groups=="non-Delta - symptomatic",],aes(x=Nabs_Ratio,y=VE_pred),col="#00BFC4")+
  geom_line(data=mL[mL$Groups=="non-Delta - infection",],aes(x=Nabs_Ratio,y=VE_pred),col="#7CAE00")+
  geom_ribbon(data=mL[mL$Groups=="Delta - symptomatic",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="#C77CFF",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="Delta - infection",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="#F8766D",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="non-Delta - symptomatic",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="#00BFC4",alpha=0.2)+
  geom_ribbon(data=mL[mL$Groups=="non-Delta - infection",],aes(x=xval,ymin=lower_bound,ymax=upper_bound),fill="#7CAE00",alpha=0.2)+
  geom_text(data=tlab,aes(x=Nabs_Ratio,y=VE,label=text1))+
  geom_text(data=tlab2,aes(x=Nabs_Ratio,y=VE,label=text1),size=5,hjust=0)+
  labs(col="Variant-Endpoint", shape="Variant-Grouping")+
  facet_wrap(.~Delta,nrow=2, scales="free_y")+theme(strip.text.x = element_blank())

#Transmission model - Fig S1
f3=glm(VE~log2(Nabs_Ratio),data=mAllt,family=binomial,weights=N_eff);summary(f3)

# Generate predicted protection vs transmission given infection for waned, boosted nAbs
p_tb=predict(f3,newdata=data.frame(Nabs_Ratio=boost_ratio*Nabs_Pfiz_delta_inf,Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Pfizer",Study="Boost_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=boost_ratio*Nabs_Pfiz_delta_inf,
                       VE=plogis(unlist(p_tb[1])))
mAllt[mAllt$Study=="Boost_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tb[1])+c(-1,1)*1.96*unlist(p_tb[2]))

p_tw=predict(f3,newdata=data.frame(Nabs_Ratio=wane_ratio*Nabs_Pfiz_delta_inf,Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Pfizer",Study="Wane_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=wane_ratio*Nabs_Pfiz_delta_inf,
                       VE=plogis(unlist(p_tw[1])))
mAllt[mAllt$Study=="Wane_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tw[1])+c(-1,1)*1.96*unlist(p_tw[2]))

p_tmb=predict(f3,newdata=data.frame(Nabs_Ratio=Mod_boost_ratio*Nabs_Mod_delta,Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Moderna",Study="Mod_Boost_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=Mod_boost_ratio*Nabs_Mod_delta,
                       VE=plogis(unlist(p_tmb[1])))
mAllt[mAllt$Study=="Mod_Boost_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tmb[1])+c(-1,1)*1.96*unlist(p_tmb[2]))

p_tmw=predict(f3,newdata=data.frame(Nabs_Ratio=Mod_wane_ratio*Nabs_Mod_delta,Groups="Delta - transmission"),type="link",se.fit=T)
mAllt=mAllt%>% add_row(Vaccine="Moderna",Study="Mod_Wane_Trans",Variant_Grouping="Delta - Prediction",Endpoint="Transmission",Groups="Delta - transmission",
                       Nabs_Ratio=Mod_wane_ratio*Nabs_Mod_delta,
                       VE=plogis(unlist(p_tmw[1])))
mAllt[mAllt$Study=="Mod_Wane_Trans",c("Lower_avg","Upper_avg")]=plogis(unlist(p_tmw[1])+c(-1,1)*1.96*unlist(p_tmw[2]))


tlabT=data.frame(text1=c("Astrazeneca orig","Pfizer orig","Pfizer boosted","Pfizer waned","Moderna boosted","Moderna waned"),
                Nabs_Ratio=c(mAllt$Nabs_Ratio[mAllt$Vaccine=="Astra"],c(1,boost_ratio,wane_ratio)*Nabs_Pfiz_delta_inf,Mod_boost_ratio*Nabs_Mod_delta,Mod_wane_ratio*Nabs_Mod_delta),
                VE=mAllt$VE-c(0.02,rep(.03,5)))

#Create fitted lines & CIs for VE for transmission Fig S1
mLt=data.frame(Nabs_Ratio=xval,Groups="Delta - transmission")
mLt$VE_logit=unlist(predict(f3,newdata=mLt, type="link", se.fit = TRUE)[1])
mLt$VE_se=unlist(predict(f3,newdata=mLt, type="link", se.fit = TRUE)[2])
mLt$fit_line = plogis(mLt$VE_logit)
mLt$upper_bound = plogis(mLt$VE_logit + (1.96 * mLt$VE_se))
mLt$lower_bound = plogis(mLt$VE_logit- (1.96 * mLt$VE_se))

#Fig S1 Plot VE for transmission
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
        axis.text=element_text(size=15),legend.title=element_blank(),
        legend.position = "none")+
  geom_text(data=tlabT,aes(x=Nabs_Ratio,y=VE,label=text1),size=4,hjust=0)
  

m2=read.csv("nAbs_waning.csv") #Relative (0-1) antibody titers over time

#Fit w/ nls to set priors; model g5 is best
g1=nls(nAbs~c0*exp(c1*Day)+c2,start=list(c0=.8,c1=-0.02,c2=.2),
       data=m2,trace = F,control = list(warnOnly=T));summary(g1)
g2=nls(nAbs~c0*exp(c1*Day+c3*Day*Moderna+c4*Infection*Day)+c2,
       start=list(c0=.8,c1=-0.02,c2=.2,c3=0.01,c4=0.01),
       data=m2,trace = F,control = list(warnOnly=T));summary(g2)
g3=nls(nAbs~(c0+c7*Infection+c8*Moderna)*exp(c1*Day+c3*Day*Moderna+c4*Infection*Day)+c2+c5*Moderna+c6*Infection,
       start=list(c0=.8,c1=-0.02,c2=.2,c3=0.01,c4=0.01,c5=0.01,c6=0.01,c7=0.01,c8=0.01),
       data=m2,trace = F,control = list(warnOnly=T));summary(g3)
g4=nls(nAbs~(c0+c7*Infection)*exp(c1*Day+c3*Day*Moderna)+c2+c6*Infection,
       start=list(c0=.8,c1=-0.02,c2=.2,c3=0.01,c6=0.01,c7=0.01),
       data=m2,trace = F,control = list(warnOnly=T));summary(g4)
g5=nls(nAbs~(c0+c7*Infection)*exp(c1*Day+c3*Day*Moderna+c4*Infection*Day)+c2+c6*Infection,#best model
       start=list(c0=.8,c1=-0.02,c2=.2,c3=0.01,c4=0.01,c6=0.01,c7=0.01),
       data=m2,trace = F,control = list(warnOnly=T));summary(g5) #Table S4
AIC(g1,g2,g3,g4,g5)

#Fit bayesian model to estimate CIs for waning antibody models
options(mc.cores = parallel::detectCores()) #speeds up stan
#priors set to be much larger than 95% CIs from nls fit (g5)
pr1 = prior(normal(.9, 1), nlpar = "c0")+ prior(normal(-.02, .02), nlpar = "c1") +  
  prior(normal(0.1, .3), nlpar = "c2")+prior(normal(0.001, 0.03), nlpar = "c3")+
  prior(normal(0.01, 0.03), nlpar = "c4")+prior(normal(0.1, .3), nlpar = "c6")+
  prior(normal(-0.1, 0.3), nlpar = "c7")

#model fit often has a small number of divergent transitions & sometimes transitions which exceed max treedepth
#but model runs without either produce nearly identical parameter estimates and CIs
f6 = brm(bf(nAbs~(c0+c7*Infection)*exp(c1*Day+c3*Day*Moderna+c4*Infection*Day)+c2+c6*Infection, 
            c0+c1+c2+c3+c4+c6+c7 ~ 1, nl = TRUE),data = m2, prior = pr1,
         chains = 4, iter = 5000,control = list(adapt_delta = 0.8))
summary(f6);fixef(f6)

#Create fitted values for Fig 2
p1=data.frame(Day=rep(0:365,3))
prval1=predict(f6,newdata=data.frame(Day=p1$Day,Moderna=rep(c(0,1,0),each=366),
                                     Infection=rep(c(0,0,1),each=366)))
p1$nAbs=prval1[,"Estimate"];p1$lower=prval1[,"Q2.5"];p1$upper=prval1[,"Q97.5"]
p1$Vaccine=rep(c("Pfizer","Moderna", "Infection"),each=366)

#Making absolute neutralizing antibody predicted lines - rescaling of relative values
q1=p1;q1$response="absolute_nAbs";p1$response="relative_nAbs"
q1[q1$Vaccine=="Pfizer",c("nAbs","lower","upper")]=q1[q1$Vaccine=="Pfizer",c("nAbs","lower","upper")]*Nabs_Pfiz_delta_inf
q1[q1$Vaccine=="Infection",c("nAbs","lower","upper")]=q1[q1$Vaccine=="Infection",c("nAbs","lower","upper")]*mAll$Nabs_Ratio[mAll$Vaccine=="Convalescent"&mAll$Study=="Pouwels (delta)"]
q1[q1$Vaccine=="Moderna",c("nAbs","lower","upper")]=q1[q1$Vaccine=="Moderna",c("nAbs","lower","upper")]*mAll$Nabs_Ratio[mAll$Study=="Nasreen(delta)"&mAll$Vaccine=="Moderna"]

#Making absolute neutralizing antibody points - rescaling of relative values
n1=m2;
n1[n1$Vaccine=="Pfizer",c("nAbs")]=n1[n1$Vaccine=="Pfizer",c("nAbs")]*Nabs_Pfiz_delta_inf
n1[n1$Vaccine=="Infection",c("nAbs")]=n1[n1$Vaccine=="Infection",c("nAbs")]*mAll$Nabs_Ratio[mAll$Vaccine=="Convalescent"&mAll$Study=="Pouwels (delta)"]
n1[n1$Vaccine=="Moderna",c("nAbs")]=n1[n1$Vaccine=="Moderna",c("nAbs")]*mAll$Nabs_Ratio[mAll$Study=="Nasreen(delta)"&mAll$Vaccine=="Moderna"]
m2$response="relative_nAbs";n1$response="absolute_nAbs"
n1=rbind(m2,n1)

#Making protection vs time dataframe
r1=q1;r1$response="VE";
names(r1)[2]="nAbs2" #renaming nAbs in r1 to nAbs2 b/c nAbs is column used for Y in plot
r1=rbind(r1,data.frame(Day=0:365,nAbs2=NA,lower=NA,upper=NA,Vaccine="Hybrid",response="VE"))
r1$nAbs2[r1$Vaccine=="Hybrid"]=hybrid_ratio*r1$nAbs2[r1$Vaccine=="Pfizer"]

#Column header for protection is "nAbs" to facilitate plotting
r1$nAbs=predict(f2,newdata=data.frame(Nabs_Ratio=r1$nAbs2,Groups="Delta - infection"),type="response")
r1$VE_logit=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=r1$nAbs2,Groups="Delta - infection"),type="link",se.fit=T)[1])
r1$VE_se=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=r1$nAbs2,Groups="Delta - infection"),type="link",se.fit=T)[2])
r1$upper=plogis(r1$VE_logit + (1.96 * r1$VE_se))
r1$lower=plogis(r1$VE_logit - (1.96 * r1$VE_se))
r1=r1[,c("Day","nAbs","lower","upper","Vaccine","response")] #eliminating working columns to bind to other dataframes
p1=rbind(p1,q1,r1) #combining absolute and relative antibodies and protection into 1 dataframe for plotting

#Panel labels
tlab2=data.frame(response=c("absolute_nAbs","relative_nAbs","VE"),
                 text1=c("A: Absolute neutralizing antibody titers","B: Relative neutralizing antibody titers",
                         "C: Protection against infection"),Day=c(300),nAbs=c(1.5,1.1,.93))
#make point for boosted
br1=mAll$Study=="Boosted"
bp=data.frame(Day=365,nAbs=mAll$VE[br1],lower=mAll$Lower_avg[br1],upper=mAll$Upper_avg[br1],response="VE")
#Figure 2
ggplot()+
  geom_point(data=n1,aes(x=Day,y=nAbs,col=Vaccine),size=3)+
  geom_ribbon(data=p1,aes(x=Day,y=nAbs,ymin=lower, ymax=upper,fill=Vaccine),colour=NA,alpha=0.1)+
  geom_line(data=p1,aes(x=Day,y=nAbs,col=Vaccine))+theme_few()+
  xlab("Days since infection/vaccination")+ylab("Neutralizing antibody titres or Protection against infection")+
  labs(col="Immunity", fill="Immunity")+
  theme(axis.title=element_text(size=20),#legend.position = c(.1, .80),
        axis.text=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=20))+
  facet_wrap(.~response,nrow=3, scales="free_y")+
  theme(strip.text.x = element_blank())+
  geom_text(data=tlab2,aes(x=Day,y=nAbs,label=text1),size=5)+
  geom_point(data=bp,aes(x=Day,y=nAbs),size=3,col=c("#C77CFF","#00BFC4"))+
  geom_text(data=bp,aes(x=Day,y=nAbs+0.01*c(1,-1),label=c("Pfizer Boosted","Moderna Boosted")),size=4,hjust=1.1)+
  geom_errorbar(data=bp,aes(x=Day,ymin=lower, ymax=upper), width=0,col="#C77CFF")+
  theme(legend.position="top")

#Read in data set with daily cases, new fully vaccinated, and deaths w/ row 1 being most recent data
case_vax_data=read.csv("Cases_Vax_Death.csv") #load case and vax data

#Start date (first day with data) and end date (last day with data), first day vaccine given, and population
start_date = '2020-01-23'; end_date = '2021-11-21';first_vax="2020-12-14";population=332858873

#Fraction of mRNA vaccines given that are pfizer
frac_pfizer=(1.503/2.503);frac_moderna=1-frac_pfizer

#Add columns with date, number of infections, and proportion of 2nd doses given on a specific day
case_vax_data$date=seq(as.Date(end_date), as.Date(start_date), by = "-1 days")

#Daily and cumulative infections
case_vax_data$inf=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$inf[i]=case_vax_data$Cases[i]*4
}

case_vax_data$cum_inf=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_inf[i]=sum(case_vax_data$inf[i:nrow(case_vax_data)])
}

#Proportion vaccinated each day
case_vax_data$prop_vax=case_vax_data$Vaccinated/sum(case_vax_data$Vaccinated)

#Create cumulative infected and cumulative fully vaccinated columns
case_vax_data$cum_vax=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_vax[i]=sum(case_vax_data$Vaccinated[i:nrow(case_vax_data)])
}

#Create cumulative uninfected and cumulative unvaccinated columns
case_vax_data$cum_unvax=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_unvax[i]=(population-case_vax_data$cum_vax[i])
}

case_vax_data$cum_uninf=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$cum_uninf[i]=(population-case_vax_data$cum_inf[i])
}

#Daily case ratios of unvaccinated/vaccinated using CDC data
#Load file
vacc_irr=read.csv("irr_unvax_vax.csv")

#GAM
g0=gam(age_adj_irr ~s(time),
       data = vacc_irr,method = "REML");summary(g0)

#Plot data with model
ggplot(data = vacc_irr,
       aes(x=time, y=age_adj_irr)) + geom_point()+
  geom_smooth(method="gam",formula = y ~ s(x))

#Build data frame using predict
fitted_irr=data.frame(time=seq(-110,232,by=1),IRR=predict(g0,newdata=data.frame(time=seq(-110,232,by=1))))

#Reverse time order
#attach(fitted_irr)
fitted_irr = fitted_irr[order(-fitted_irr$time),]
#detach(fitted_irr)

#Remove time
case_ratios=data.frame(IRR=fitted_irr["IRR"])

#Start and end dates of case ratio data (leave as is if using CDC data)
case_ratio_start="2020-12-14";case_ratio_end="2021-11-21"

case_ratios$date=as.Date(NA)
#Create dataframe with case ratios of unvaccinated/vaccinated 
case_ratios$date=seq(as.Date(case_ratio_end), as.Date(case_ratio_start), by = "-1 days")
#case_ratios_recent = data.frame(IRR=IRR_now,date=seq(as.Date(end_date), as.Date(case_ratio_end)+1, by = "-1 days"))
case_ratios_prevax = data.frame(IRR=0, date= seq(as.Date(case_ratio_start)-1, as.Date(start_date), by = "-1 days"))
case_ratios_all=rbind(case_ratios,case_ratios_prevax)

#Add case ratios to main dataframe
case_vax_data=cbind(case_vax_data,IRR=case_ratios_all$IRR)

#Create columns splitting infections between unvaccinated and vaccinated
case_vax_data$daily_inf_in_vax=NA
case_vax_data$daily_inf_in_unvax=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$daily_inf_in_vax=(case_vax_data$inf/(case_vax_data$IRR*(case_vax_data$cum_unvax/case_vax_data$cum_vax)+1))
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
case_vax_data$daily_super=NA

for (i in 1:nrow(case_vax_data)) {
  case_vax_data$daily_super[nrow(case_vax_data)+1-i]=(case_vax_data$cum_unvax_pos[nrow(case_vax_data)+2-i]/(case_vax_data$cum_unvax_neg[nrow(case_vax_data)+2-i]+case_vax_data$cum_unvax_pos[nrow(case_vax_data)+2-i])*case_vax_data$Vaccinated[nrow(case_vax_data)+1-i]+case_vax_data$daily_inf_in_vax[nrow(case_vax_data)+1-i])
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

#Create columns with proportion of all individuals with hybrid immunity, vaccinated/seronegative, and
#unvacccinated/seropositive (calculated two ways: 1. using case ratios and 2. using death data)
case_vax_data$prop_super=(case_vax_data$daily_super/sum(case_vax_data$daily_super))
case_vax_data$prop_unvax_pos=(case_vax_data$daily_unvax_pos/sum(case_vax_data$daily_unvax_pos))
case_vax_data$prop_vax_neg_adj=(case_vax_data$daily_vax_neg_adj/sum(case_vax_data$daily_vax_neg_adj))

#Create dataframe to calculate VEs
#Pull relevant columns from case_vax_data
VE_data=data.frame(time=c(rep(0,15),1:nrow(case_vax_data)),prop_inf_unvax=NA,prop_super=NA,prop_pfizer=NA,prop_moderna=NA)
VE_data$prop_inf_unvax[8:(nrow(case_vax_data)+8-24)]=case_vax_data$prop_unvax_pos[24:nrow(case_vax_data)]
VE_data$prop_super[15:(nrow(case_vax_data)+15-1)]=case_vax_data$prop_super
VE_data$prop_pfizer[8:(nrow(case_vax_data)+8-1)]=case_vax_data$prop_vax_neg_adj
VE_data$prop_moderna[1:(nrow(case_vax_data)+1-1)]=case_vax_data$prop_vax_neg_adj
VE_data[is.na(VE_data)] = 0

#Deduction in nAbs from Delta
Nabs_Pfiz_delta_inf=mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"]

#Add columns with predicted nAbs for each class
VE_data$nAbs_inf_unvax=predict(g5,newdata=data.frame(Day=VE_data$time,Infection=1,Moderna=0))/Delta_WT_nAbs_Ratio  
VE_data$nAbs_pfizer=Nabs_Pfiz_delta_inf*
  predict(g5,newdata=data.frame(Day=VE_data$time,Infection=0,Moderna=0))
VE_data$nAbs_moderna=mAll$Nabs_Ratio[mAll$Study=="Nasreen(delta)"&mAll$Vaccine=="Moderna"]*
  predict(g5,newdata=data.frame(Day=VE_data$time,Infection=0,Moderna=1))
VE_data$nAbs_super=VE_data$nAbs_pfizer*hybrid_ratio

#CIs for vaccine immunity w/ waning 
#CIs for vaccine
#Pfizer and Moderna Waning
nds=10000
qdraws=runif(nds)
y1=plogis(matrix(qnorm(rep(qdraws,each=nrow(VE_data)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_pfizer,
                                                                                                Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )

qdraws0=runif(nds)
y2=plogis(matrix(qnorm(rep(qdraws0,each=nrow(VE_data)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_moderna,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_moderna,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )
x1=VE_data$prop_pfizer
x2=VE_data$prop_moderna
x3=x1%*%y1
x4=x2%*%y2
CI_inf_vax1y=(1.503/2.503)*x3+(1/2.503)*x4
CI_inf_vaxy=as.data.frame(t((1.503/2.503)*x3+(1/2.503)*x4))
#quantile(CI_inf_vax1y, probs = c(0.025,.50,0.975))

qdraws1=runif(nds)
y5=plogis(matrix(qnorm(rep(qdraws1,each=nrow(VE_data)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_pfizer,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )

qdraws2=runif(nds)
y6=plogis(matrix(qnorm(rep(qdraws2,each=nrow(VE_data)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_moderna,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_moderna,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )

x5=x1%*%y5
x6=x2%*%y6
CI_trans_vax1y=y7=(1.503/2.503)*x5+(1/2.503)*x6
CI_trans_vaxy=as.data.frame(t((1.503/2.503)*x5+(1/2.503)*x6))
#quantile(CI_trans_vax1y, probs = c(0.025,.50,0.975))

#Super Waning
qdraws3=runif(nds)
y8=plogis(matrix(qnorm(rep(qdraws3,each=nrow(VE_data)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_super,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_super,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )

z3=VE_data$prop_super 
CI_inf_super1=z3%*% y8
CI_inf_super=as.data.frame(t(z3%*% y8))
quantile(CI_inf_super1, probs = c(0.025,.50,0.975))

qdraws4=runif(nds)
y9=plogis(matrix(qnorm(rep(qdraws4,each=nrow(VE_data)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_super,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_super,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(VE_data)) )

CI_trans_super1=z3%*% y9
CI_trans_super=as.data.frame(t(z3%*% y9))
quantile(CI_trans_super1, probs = c(0.025,.50,0.975))

#Waning natural infection
nds=10000
qdraws5=runif(nds)
z0=plogis(matrix(qnorm(rep(qdraws5,each=nrow(VE_data)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_inf_unvax,
                                                                                                 Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_inf_unvax,
                                                               
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2]) ),nrow=nrow(VE_data)) )

z1=VE_data$prop_inf_unvax  
CI_inf1y=z1%*% z0
CI_infy=as.data.frame(t(z1%*% z0))
#quantile(CI_inf1y, probs = c(0.025,.50,0.975))

qdraws6=runif(nds)
z2=plogis(matrix(qnorm(rep(qdraws,each=nrow(VE_data)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_inf_unvax,
                                                                                                Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=VE_data$nAbs_inf_unvax,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2]) ),nrow=nrow(VE_data)) )

CI_trans1y=z1%*% z2
CI_transy=as.data.frame(t(z1%*% z2))
#quantile(CI_trans1y, probs = c(0.025,0.50,0.975))


#Read in NZ data
daily_NZ_vax=read.csv("daily_NZ.csv")
daily_NZ_vax$nAbs_pfizer=Nabs_Pfiz_delta_inf*
  predict(g5,newdata=data.frame(Day=daily_NZ_vax$time,Infection=0,Moderna=0)) 

#New Zealand waning in vaccinated
nds=10000
qdraws=runif(nds)
n1=plogis(matrix(qnorm(rep(qdraws,each=nrow(daily_NZ_vax)),mean=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_NZ_vax$nAbs_pfizer,
                                                                                                     Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f2,newdata=data.frame(Nabs_Ratio=daily_NZ_vax$nAbs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_NZ_vax)) )

m3=daily_NZ_vax$prop_vax
CI_inf_vax1_NZ=m3%*%n1
CI_inf_vaxy_NZ=as.data.frame(t(m3%*%n1))
quantile(CI_inf_vax1_NZ, probs = c(0.025,.50,0.975))

qdraws1=runif(nds)
n2=plogis(matrix(qnorm(rep(qdraws1,each=nrow(daily_NZ_vax)),mean=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_NZ_vax$nAbs_pfizer,
                                                                                                      Groups="Delta - infection"),type="link",se.fit=T)[1]),
                       sd=unlist(predict(f3,newdata=data.frame(Nabs_Ratio=daily_NZ_vax$nAbs_pfizer,
                                                               Groups="Delta - infection"),type="link",se.fit=T)[2])),nrow=nrow(daily_NZ_vax)) )

CI_trans_vax1y_NZ=m3%*%n2
CI_trans_vaxy_NZ=as.data.frame(t(m3%*%n2))
quantile(CI_trans_vax1y_NZ, probs = c(0.025,.50,0.975))



parsA=data.frame(
  R0=3.7,
  f_n=0.41, #fraction not vaccinated
  f_v=0.59, #fraction vaccinated
  f_pu=0.715, #fraction prev infected in US by vacccination status - from 44M cases, 4.2 inf/case
  f_pv=0.464, #Used cdc data on RR unvaccinated v. vaccinated
  VE_ii=(quantile(CI_inf1y,probs=c(0.5))), #prev infected waned Inf
  VE_ti=(quantile(CI_trans1y,probs=c(0.5))), #prev infected waned trans
  VE_iv=(quantile(CI_inf_vax1y,probs=c(0.5))), #mAll$VE[mAll$Study=="Waned"], #vacc waned infect
  VE_tv=(quantile(CI_trans_vax1y,probs=c(0.5))), #mAllt$VE[mAllt$Study=="Waned_Trans"], #vacc waned tran
  VE_is=quantile(CI_inf_super1, probs = c(.5)), #Super VEi
  VE_ts=quantile(CI_trans_super1, probs = c(.5)), #Super VEt
  VE_ib=mAll$VE[mAll$Study=="Boosted"&mAll$Vaccine=="Pfizer"], #boost inf
  VE_tb=mAllt$VE[mAllt$Study=="Boost_Trans"] #boost trans
)

#Rt function 
Rt=function(parsR){
  with(parsR,
       R0*((f_n*(1-f_pu))+(f_n*f_pu)*(1-VE_ii)*(1-VE_ti)+(1-f_pv)*(f_v-f_b)*(1-VE_iv)*(1-VE_tv)+(f_pv)*(f_v-f_b)*(1-VE_is)*(1-VE_ts)+(f_b)*(1-VE_ib)*(1-VE_tb)))}


pars1A=data.frame(f_b=seq(0,0.59,by=0.01),Rtb=Rt(cbind(parsA, f_b=seq(0.00,0.59,by=0.01))),Scenario="R0=3.7;  59% vacc., 56.9% prev. inf.")

#Scenario B - 100% vacc., R0 = 7
parsB=parsA; parsB[,c("R0","f_n","f_v","f_pu","f_pv")]=c(7,0,1,0,0.569)
pars1B=data.frame(f_b=seq(0,1,by=0.01),Rtb=Rt(cbind(parsB, f_b=seq(0,1,by=0.01))),Scenario="R0=7.0;  100% vacc., 56.9% prev. inf.")

#Scenario C - 56% vacc., R0 = 7
parsC=parsA;parsC["R0"]=7
pars1C=data.frame(f_b=seq(0,.59,by=0.01),Rtb=Rt(cbind(parsC, f_b=seq(0,.59,by=0.01))),Scenario="R0=7.0;  59% vacc., 56.9% prev. inf.")

#Scenario D - 75% vacc., R0 3.7, Prev inf = 28%
parsE=parsA;parsE[,c("R0","f_n","f_v","f_pu","f_pv")]=c(3.7,0.25,0.75,0.421,0.238)
pars1E=data.frame(f_b=seq(0,.75,by=0.01),Rtb=Rt(cbind(parsE, f_b=seq(0,.75,by=0.01))),Scenario="R0=3.7;  75% vacc., 28.5% prev. inf.")

#Scenario E - 69% vacc., R0 3.7, Prev inf = 28%
parsF=parsA; parsF[,c("R0","f_n","f_v","f_pu","f_pv","VE_iv","VE_tv")]=c(3.7,0.31,0.69,0.021,0.005,(quantile(CI_inf_vax1_NZ,probs=c(0.5))),(quantile(CI_trans_vax1y_NZ,probs=c(0.5))))
pars1F=data.frame(f_b=seq(0,0.69,by=0.01),Rtb=Rt(cbind(parsF, f_b=seq(0,.69,by=0.01))),Scenario="R0=3.7;  69% vacc., 1% prev. inf.")



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

draws_NZ=data.frame()
draws_NZ=cbind("VE_ii"=unlist(CI_infy),
               "VE_ti"=unlist(CI_transy),
               "VE_iv"=unlist(CI_inf_vaxy_NZ),
               "VE_tv"=unlist(CI_trans_vaxy_NZ),
               "VE_is"=unlist(CI_inf_super),
               "VE_ts"=unlist(CI_trans_super),
               "VE_ib"=(plogis(rnorm(nd1,unlist(p_b[1]),unlist(p_b[2])))),
               "VE_tb"=(plogis(rnorm(nd1,unlist(p_tb[1]),unlist(p_tb[2])))))


for (i in 1:101) { #calculate CIs for each scenario
  if (i<61) { #2 scenarios that have max boosters = 56%
    pars1A[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsA[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1A$f_b[i])))
    pars1C[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsC[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1C$f_b[i])))
  }
  if(i<77){ #scenario that has max booster = 75%
    pars1E[i,c("lower","upper")]=
      quantile(probs = c(0.025,0.975),x=Rt(cbind(parsE[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1E$f_b[i])))
  }
  pars1B[i,c("lower","upper")]=
    quantile(probs = c(0.025,0.975),x=Rt(cbind(parsB[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1B$f_b[i])))
}

for (i in 1:70) {
  pars1F[i,c("lower","upper")]=
    quantile(probs = c(0.025,0.975),x=Rt(cbind(parsF[c("R0","f_n","f_v","f_pu","f_pv")],draws_NZ,f_b=pars1F$f_b[i])))
  
}


#Add point for shots going to unvaxxed
parsD=parsC
parsD[,c("R0","f_n","f_v","f_pv")]=c(7,0.115,0.885,0.668);parsD$f_b=0
pars1D=data.frame(f_b=0,Rtb=Rt(parsD),Scenario="R0=7.0; f_vacc=88.5%,")
pars1D[c("lower","upper")]=
  quantile(probs = c(0.025,0.975),x=Rt(cbind(parsD[c("R0","f_n","f_v","f_pu","f_pv")],draws,f_b=pars1D$f_b)))

Rt_dat=data.frame(rbind(pars1A,pars1B,pars1C,pars1E,pars1F)) #group data into single dataframe

#Fig 3
ggplot(data=Rt_dat,aes(x=f_b,y=Rtb),)+
  geom_line(aes(col=Scenario))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Scenario),alpha=0.15)+
  geom_point(data=pars1D,aes(x=f_b,y=Rtb),size=5)+
  theme_few()+scale_x_continuous(limits=c(0, 1.01), expand = c(.01, 0)) +
  geom_hline(yintercept=1, linetype="longdash")+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Fraction of population receiving 3rd dose")+
  ylab("Effective Reproductive Number (Rt)")+
  labs(fill="Scenario")+scale_fill_manual(values=c("red", "grey","blue", "yellow3","green"))+
  scale_color_manual(values=c("red", "grey","blue","yellow3" ,"green"))+
  theme(strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85) )+
  geom_text(data=pars1D,aes(x=f_b,y=Rtb,label="R0 = 7; 88.5% vacc., no 3rd doses"),size=7,hjust=-0.05)+
  theme(axis.title=element_text(size=25),#legend.position = c(.1, .80),
        axis.text=element_text(size=25),legend.text=element_text(size=20),
        legend.title=element_text(size=20))

max(pars1A$Rtb);min(pars1A$Rtb);min(pars1A$Rtb)/max(pars1A$Rtb) #R0=3.7;  59% vacc., 56.9% prev. inf.
max(pars1B$Rtb);min(pars1B$Rtb);min(pars1B$Rtb)/max(pars1B$Rtb) #R0=7.0;  100% vacc., 56.9% prev. inf.
min(pars1B$f_b[pars1B$Rtb<1]) #frac boosted needed for Rt<1 for CA scenario
max(pars1C$Rtb);min(pars1C$Rtb);min(pars1C$Rtb)/max(pars1C$Rtb) #R0=7.0;  59% vacc., 56.9% prev. inf.
max(pars1E$Rtb);min(pars1E$Rtb);min(pars1E$Rtb)/max(pars1E$Rtb) #R0=3.7;  75% vacc., 28.2% prev. inf.
min(pars1E$f_b[pars1E$Rtb<1]) #frac boosted needed for Rt<1 for CA scenario
max(pars1F$Rtb);min(pars1F$Rtb);min(pars1F$Rtb)/max(pars1F$Rtb) #R0=3.7;  69% vacc., 1% prev. inf.

#Estimate Rt for 100% vaccinated by no waning
parsG=parsB; parsG[,c("VE_iv","VE_tv")]=c(unlist(predict(f2, newdata=data.frame(Nabs_Ratio=mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - infection"), type="response")),
                                          unlist(predict(f3, newdata=data.frame(Nabs_Ratio=mAll$Nabs_Ratio[mAll$Vaccine=="Pfizer"&mAll$Study=="Pouwels (delta)"],Groups="Delta - infection"), type="response")))
parsG$f_b=0
pars1G=data.frame(f_b=0,Rtb=Rt(parsG),Scenario="R0=7.0;  100% vacc., 56.4% prev. inf., no waning")



#Fig S3 Plot infections, cases, deaths
#Crude IFR estimated by dividing total 4.2*total cases through 24 days prior by total deaths
#Infections estimated for day t-24 by dividing deaths on day n by IFR (0.003993)
#daily_cases=read.csv("Cases_Vax_Death.csv") #already read in above
daily_cases$Infections=NA
for (i in 1:(nrow(daily_cases)-23)) {
  daily_cases$Infections[i+23]=as.integer(daily_cases$Deaths[i]/0.003993)
}
#daily_cases[is.na(daily_cases)] = 0
daily_cases$date=seq(as.Date('2021-11-21'), as.Date('2020-1-23'), by = "-1 days")

ggplot()+
  geom_line(data=daily_cases,aes(x=date,y=Cases/1000,color="Cases"))+
  geom_line(data=daily_cases,aes(x=date,y=Infections/1000,color="Infections"))+
  geom_line(data=daily_cases,aes(x=date,y=Deaths,color="Deaths"))+
  ylab("COVID-19 deaths or cases and infections in 1000s")+
  scale_x_date(date_labels = "%b %y",date_breaks = "2 month")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank() )+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  labs(color = "")+
  theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=20))


#Fig S2 Daily vaccinated plot
ggplot()+
  geom_line(data=daily_cases[1:343,],aes(x=date,y=Vaccinated/1000))+
  ylab("People fully vaccinated against COVID-19, 1000s")+  
  scale_x_date(date_labels = "%b %y",date_breaks = "1 month")+
  theme_bw()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank())+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=20))


#Fig S4 Rate ratio plot
vacc_irr$date=seq(as.Date('2021-09-27'), as.Date('2021-04-04'), by = "-7 days")
ggplot(data = vacc_irr,aes(x=date, y=age_adj_irr)) +
  geom_point(size=4)+geom_smooth(method="gam",formula = y ~ s(x),col="black")+
  ylab("Ratio of COVID-19 cases in unvaccinated v. vaccinated")+
  scale_y_continuous(breaks =seq(0,16,2),limits=c(0,16))+
  scale_x_date(date_labels = "%b %Y",date_breaks = "1 month")+
  theme_few()+
  theme(text = element_text(size=16),strip.text.x = element_blank(),strip.background = element_rect(colour="white", fill="white"),
        legend.position=c(.8,.85),axis.title.x=element_blank() )+
  theme(axis.title=element_text(size=25),axis.text=element_text(size=25),
        legend.text=element_text(size=20),legend.title=element_text(size=20))

#Table S5 values
#R0=3.7;  59% vacc., 56.9% prev. inf.
parsA$f_pv*parsA$f_v #Fraction vaccinated and infected
(1-parsA$f_pv)*parsA$f_v #Fraction vaccinated and not infected
parsA$f_pu*parsA$f_n #Fraction unvaccinated and infected
(1-parsA$f_pu)*parsA$f_n #Fraction unvaccinated and not infected
#R0=7;  100% vacc., 56.9% prev. inf.
parsB$f_pv*parsB$f_v #Fraction vaccinated and infected
(1-parsB$f_pv)*parsB$f_v #Fraction vaccinated and not infected
parsB$f_pu*parsB$f_n #Fraction unvaccinated and infected
(1-parsB$f_pu)*parsB$f_n #Fraction unvaccinated and not infected
#R0=7;  59% vacc., 56.9% prev. inf.
parsC$f_pv*parsC$f_v #Fraction vaccinated and infected
(1-parsC$f_pv)*parsC$f_v #Fraction vaccinated and not infected
parsC$f_pu*parsC$f_n #Fraction unvaccinated and infected
(1-parsC$f_pu)*parsC$f_n #Fraction unvaccinated and not infected
#R0=3.7;  75% vacc., 28.5% prev. inf.
parsE$f_pv*parsE$f_v #Fraction vaccinated and infected
(1-parsE$f_pv)*parsE$f_v #Fraction vaccinated and not infected
parsE$f_pu*parsE$f_n #Fraction unvaccinated and infected
(1-parsE$f_pu)*parsE$f_n #Fraction unvaccinated and not infected
#R0=3.7;  69% vacc., 1% prev. inf.
parsF$f_pv*parsF$f_v #Fraction vaccinated and infected
(1-parsF$f_pv)*parsF$f_v #Fraction vaccinated and not infected
parsF$f_pu*parsF$f_n #Fraction unvaccinated and infected
(1-parsF$f_pu)*parsF$f_n #Fraction unvaccinated and not infected

#Table S6 values: mean & CIs
quantile(CI_inf1y,probs=c(0.025,0.5,0.975)) #prev infected waned Inf
quantile(CI_trans1y,probs=c(0.025,0.5,0.975)) #prev infected waned trans
quantile(CI_inf_vax1y,probs=c(0.025,0.5,0.975)) #vacc waned infect US Scenarios
quantile(CI_trans_vax1y,probs=c(0.025,0.5,0.975)) #vacc waned tran US Scenarios
quantile(CI_inf_vax1_NZ, probs = c(0.025,.50,0.975)) #vacc waned infect NZ Scenarios
quantile(CI_trans_vax1y_NZ, probs = c(0.025,.50,0.975))#vacc waned tran NZ Scenarios
quantile(CI_inf_super1, probs = c(0.025,0.5,0.975)) #Hybrid waned VEi 
quantile(CI_trans_super1, probs = c(0.025,0.5,0.975)) #Hybrid waned VEt
mAll$Lower_avg[mAll$Study=="Boosted"&mAll$Vaccine=="Pfizer"];mAll$VE[mAll$Study=="Boosted"&mAll$Vaccine=="Pfizer"];mAll$Upper_avg[mAll$Study=="Boosted"&mAll$Vaccine=="Pfizer"] #Boosted VEi
mAllt$Lower_avg[mAllt$Study=="Boost_Trans"];mAllt$VE[mAllt$Study=="Boost_Trans"];mAllt$Upper_avg[mAllt$Study=="Boost_Trans"] #Boosted VEt
