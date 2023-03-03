#library(nimble)
#library(parallel)
#library(R0)
library(dplyr)
library(bayesplot)
library(rstan)
library(ggplot2)
library(lubridate)
library(scales)
library(loo)
set.seed(10120)




gen_interval<-c(0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1)
death_interval<-c(rep(0.04,5),rep(0.09,6),rep(0.05,4),rep(0.01,6))
len_gen_int<-length(gen_interval)
len_death_int<-length(death_interval)

ifr<-0.02
# Define time varying effective reproduction number
Ret <- function(i,cp,R) { 
  return(R[sum(i>cp)+1])
}


routbreak <- function(n=100, Ret,cp=c(60,100,150,200),R=c(1.5,0.95,1.35,0.8,1.8),gen_interval,death_interval ,ifr, initial_cases = 100,N=10000000) {
  # Set up time series of incident cases
  y <- rep(0, n)
  deaths <- rep(0,n)
  len_gen_int<-length(gen_interval)
  len_death_int<-length(death_interval)
  y[1:len_gen_int]<-initial_cases
  deaths[1:len_death_int]<-0
  
  # Loop over all time points
  
  
  for (i in (len_gen_int+1):n) {
    y[i] <- rpois(1,Ret(i,cp,R)*((N-sum(y[1:(i-1)]))/N)*sum(y[(i-len_gen_int):(i-1)]*gen_interval))
  }
  print(ifr)
  for(i in (len_death_int+1):n){
    deaths[i]<-rpois(1,ifr*sum(y[(i-len_death_int):(i-1)]*death_interval))
  }
  
  
  # Data frame with the result. Assume we start on 15th of Feb
  res <- data.frame( y=y,deaths=deaths)
  
  #Done
  return(res)
}

# Generate an outbreak (no stochasticity, just the difference equation)
n_cases<-250
pop<-100000000
out <- routbreak(n=n_cases, Ret=Ret, gen_interval=gen_interval,death_interval=death_interval,initial_cases = 1000, ifr=ifr,N=pop)
#out <- out %>% mutate(ratio = y/lag(y))

cases<-as.integer(out$y)
death_cases<-as.integer(out$deaths)
sum(death_cases)/sum(cases)
I<-rep(1,n_cases)

for(i in 1:len_gen_int){
  I[i]<-sum(cases[1:i])}
for(i in (len_gen_int+1):n_cases){
  I[i]<-sum(cases[(i-len_gen_int):(i-1)]*gen_interval)}
plot(cases,xlab = 'Days')
#lines(Cases_pred,col='red')
plot(death_cases,xlab = 'Days')

plot(I,ylab='Active Infectious',xlab = 'Days')

plot(cases,xlab = 'Days')
lines(death_cases)

gen_interval2<-c(rep(0,length(death_cases)-7),0.1, 0.1, 0.2, 0.2, 0.2, 0.1, 0.1)
death_interval2<-c(rep(0,length(death_cases)-21),rep(0.04,5),rep(0.09,6),rep(0.05,4),rep(0.01,6))
ifr2<-rep(0.02,3)
ifr_cp<-c(100,150)
data<-list(N0=6,reported_cases=cases,Nobs=length(death_cases),pop=pop,deaths=death_cases,s_int_rev1=gen_interval2,s_int_rev2=death_interval2,ifr=ifr2,ifr_cp=ifr_cp,epi_start=21)

inits_4cp=list(list(e=rep(20,3),T1=50,R0=3,effect=rep(1,4),prediction0=rep(10,6)),list(e=rep(20,3),T1=50,R0=3,effect=rep(1,4),prediction0=rep(10,6)),list(e=rep(20,3),T1=50,R0=3,effect=rep(1,4),prediction0=rep(10,6)))

inits_5cp=list(list(e=rep(20,4),T1=50,R0=3,effect=rep(1,5),prediction0=rep(10,6)),list(e=rep(20,4),T1=50,R0=3,effect=rep(1,5),prediction0=rep(10,6)),list(e=rep(20,4),T1=50,R0=3,effect=rep(1,5),prediction0=rep(10,6)))

inits_3cp=list(list(e=rep(20,2),T1=50,R0=3,effect=rep(1,3),prediction0=rep(10,6)),list(e=rep(20,2),T1=50,R0=3,effect=rep(1,3),prediction0=rep(10,6)),list(e=rep(20,2),T1=50,R0=3,effect=rep(1,3),prediction0=rep(10,6)))
inits_6cp=list(list(e=rep(20,5),T1=50,R0=3,effect=rep(1,6),prediction0=rep(10,6)),list(e=rep(20,5),T1=50,R0=3,effect=rep(1,6),prediction0=rep(10,6)),list(e=rep(20,5),T1=50,R0=3,effect=rep(1,6),prediction0=rep(10,6)))

fit_4cp<-stan(
  file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\deaths SIR model 4 changepoints.stan",  # Stan program
  data = data,    # named list of data
  chains = 3,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 30000,            # total number of iterations per chain
  thin=1,           # total number of iterations per chain
  cores  = 3,              # number of cores (could use one per chain)
  refresh = 500,             # per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  init=inits_4cp )
save(fit_4cp,file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\simulation\\SIM_4cp.RData")

fit_5cp<-stan(
  file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\deaths SIR model 5 changepoints.stan",  # Stan program
  data = data,    # named list of data
  chains = 3,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 30000,            # total number of iterations per chain
  thin=1,           # total number of iterations per chain
  cores  = 3,              # number of cores (could use one per chain)
  refresh = 500,             # per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  init=inits_5cp )
save(fit_5cp,file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\simulation\\SIM_5cp.RData")




fit_3cp<-stan(
  file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\deaths SIR model 3 changepoints.stan",  # Stan program
  data = data,    # named list of data
  chains = 3,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 30000,            # total number of iterations per chain
  thin=1,           # total number of iterations per chain
  cores  = 3,              # number of cores (could use one per chain)
  refresh = 500,             # per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  init=inits_3cp )  
save(fit_3cp,file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\simulation\\SIM_3cp.RData")


fit_6cp<-stan(
  file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\deaths SIR model 6 changepoints.stan",  # Stan program
  data = data,    # named list of data
  chains = 3,             # number of Markov chains
  warmup = 15000,          # number of warmup iterations per chain
  iter = 30000,            # total number of iterations per chain
  thin=1,           # total number of iterations per chain
  cores  = 3,              # number of cores (could use one per chain)
  refresh = 500,             # per how many iterations progress is shown
  control=list(max_treedepth =15,adapt_delta=0.8),
  init=inits_6cp )  

save(fit_6cp,file="C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\simulation\\SIM_6cp.RData")

summ_fit_3cp<-summary(fit_3cp,pars=c('T1','T2','T3'))
summ_fit_5cp<-summary(fit_5cp,pars=c('T1','T2','T3','T4','T5'))
summ_fit_4cp<-summary(fit_4cp,pars=c('T1','T2','T3','T4'))
traceplot(fit_3cp,pars=c('T1','T2','T3'))
traceplot(fit_4cp,pars=c('T1','T2','T3','T4'))
traceplot(fit_5cp,pars=c('T1','T2','T3','T4','T5'))
traceplot(fit_6cp,pars=c('T1','T2','T3','T4','T5','T6'))

Rt_fit<-list()
Rt_fit[[1]]<-summary(fit_3cp,pars=c('Rt'),use_cache=F)$summary
Rt_fit[[2]]<-summary(fit_4cp,pars=c('Rt'),use_cache=F)$summary
Rt_fit[[3]]<-summary(fit_5cp,pars=c('Rt'),use_cache=F)$summary
Rt_fit[[4]]<-summary(fit_6cp,pars=c('Rt'),use_cache=F)$summary

est_cases<-list()
est_cases[[1]]<-summary(fit_3cp,pars=c('prediction'),use_cache=F)$summary
est_cases[[2]]<-summary(fit_4cp,pars=c('prediction'),use_cache=F)$summary
est_cases[[3]]<-summary(fit_5cp,pars=c('prediction'),use_cache=F)$summary
est_cases[[4]]<-summary(fit_6cp,pars=c('prediction'),use_cache=F)$summary

mean_deaths<-list()
mean_deaths[[1]]<-summary(fit_3cp,pars=c('E_deaths'))$summary
mean_deaths[[2]]<-summary(fit_4cp,pars=c('E_deaths'))$summary
mean_deaths[[3]]<-summary(fit_5cp,pars=c('E_deaths'))$summary
mean_deaths[[4]]<-summary(fit_6cp,pars=c('E_deaths'))$summary


log_lik3 <- extract_log_lik(fit_3cp, parameter_name = "log_likelihood")
log_lik4 <- extract_log_lik(fit_4cp, parameter_name = "log_likelihood")
log_lik5 <- extract_log_lik(fit_5cp, parameter_name = "log_likelihood")
log_lik6 <- extract_log_lik(fit_6cp, parameter_name = "log_likelihood")

waic(log_lik3)
waic(log_lik4)
waic(log_lik5)
waic(log_lik6)

loo(log_lik3)
loo(log_lik4)
loo(log_lik5)
loo(log_lik6)

#plot for estimated cases
plot_filepath<-"C:\\Users\\pbarb\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\rstan final runs countries and states\\simulation\\plots\\"
plot_index<-c('3cp','4cp','5cp','6cp')
dates<-1:250

traceplot(fit_3cp,pars=c('T1','T2','T3'))
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_3cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_3cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

traceplot(fit_4cp,pars=c('T1','T2','T3','T4'))
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_4cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_4cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

traceplot(fit_5cp,pars=c('T1','T2','T3','T4','T5'))
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_5cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_5cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

traceplot(fit_6cp,pars=c('T1','T2','T3','T4','T5','T6'))
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_6cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_traceplot_6cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)


stan_ac(fit_3cp,pars=c('T1','T2','T3'))
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_3cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_3cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

stan_ac(fit_4cp,pars=c('T1','T2','T3','T4'))
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_4cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_4cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

stan_ac(fit_5cp,pars=c('T1','T2','T3','T4','T5'))
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_5cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_5cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)

stan_ac(fit_6cp,pars=c('T1','T2','T3','T4','T5','T6'))
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_6cp.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_autocorrelation_6cp.jpeg'),device='jpeg',dpi=300,width=10,height=9)


for(i in 1:4){
  
  dat_cases<-data.frame(real_cases=cases,est_cases=est_cases[[i]][,6],est_cases_low=est_cases[[i]][,4],est_cases_high=est_cases[[i]][,8])
  cols <- c("Median with 95% Cr.I."="black")
  shapes <- c(" "=17)
  
  print(ggplot() +
          # points
          geom_point(data=dat_cases, aes(x=dates, y=real_cases,shape=' '), colour="black",size=4) + 
          # red line
          geom_line(data=dat_cases, aes(x=dates, y=est_cases,colour="Median with 50% Cr.I."),
                    size=0.75)+
          geom_ribbon(aes(ymin = est_cases_low, ymax = est_cases_high,x=dates),data=dat_cases,size=0.75, alpha = 0.5, 
                      linetype="dashed",color="black", fill = "black")+
          scale_shape_manual(name="Reported Cases",values=shapes, 
                             guide = guide_legend(override.aes=aes(fill=NA)))+
          scale_color_manual(name="Estimated Cases",values=cols, 
                             guide = guide_legend(override.aes=aes(fill=NA)))+
          ylab("Daily Infections") + xlab("Days")+theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 30)))
  ggsave(filename=paste0(plot_filepath,'SIM_cases_',plot_index[i],'_blacked.tiff'),device='tiff',dpi=300,width=10,height=9)
  ggsave(filename=paste0(plot_filepath,'SIM_cases_',plot_index[i],'_blacked.jpeg'),device='jpeg',dpi=300,width=10,height=9)
  
  }

#plot for deaths
Deaths_real<-death_cases
for(i in 1:4){
  dat_deaths<-data.frame(real_deaths=Deaths_real,est_deaths=mean_deaths[[i]][,6],est_deaths_low=mean_deaths[[i]][,5],est_deaths_high=mean_deaths[[i]][,7])
  cols <- c("Median with 50% Cr.I."="red")
  shapes <- c(" "=17)
  
  print(ggplot() +
          # points
          geom_point(data=dat_deaths, aes(x=dates, y=real_deaths,shape=' '), colour="black") + 
          # red line
          geom_line(data=dat_deaths, aes(x=dates, y=est_deaths,colour="Median with 50% Cr.I."),
                    size=0.5)+
          geom_ribbon(aes(ymin = est_deaths_low, ymax = est_deaths_high,x=dates),data=dat_deaths, alpha = 0.25, 
                      linetype="dashed",color="red", fill = "red")+
          scale_shape_manual(name="Reported Deaths",values=shapes, 
                             guide = guide_legend(override.aes=aes(fill=NA)))+
          scale_color_manual(name="Estimated Deaths",values=cols, 
                             guide = guide_legend(override.aes=aes(fill=NA)))+
          ylab("Daily Deaths") + xlab("Days")+theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 30)))
  ggsave(filename=paste0(plot_filepath,'SIM_deaths_',plot_index[i],'.tiff'),device='tiff',dpi=300,width=10,height=9)
  ggsave(filename=paste0(plot_filepath,'SIM_deaths_',plot_index[i],'.jpeg'),device='jpeg',dpi=300,width=10,height=9)
  
  }
#plot for Rt
Rt_real<-c(rep(1.5,60),rep(0.95,40),rep(1.35,50),rep(0.8,50),rep(1.8,50))#cp=c(60,100,150,200),R=c(1.5,0.95,1.35,0.8,1.8)
for(i in 1:4){
  dat_Rt<-data.frame(real_Rt=Rt_real,est_Rt=Rt_fit[[i]][,6],est_Rt_low=Rt_fit[[i]][,4],est_Rt_high=Rt_fit[[i]][,8])
  cols <- c("Median estimate with 95% Cr.I."="black",'Real'='black')
  print(ggplot() +
          # points
          geom_line(data=dat_Rt, aes(x=dates, y=real_Rt, colour="Real"),size=1) + 
          # red line
          geom_line(data=dat_Rt, aes(x=dates, y=est_Rt,colour="Median estimate with 50% Cr.I."),
                    size=1,alpha=1,linetype="dashed",colour='black')+
          geom_ribbon(aes(ymin = est_Rt_low, ymax = est_Rt_high,x=dates),data=dat_Rt,size=1, alpha = 0.25, 
                      linetype="dashed",
                      colour="black",fill='black')+
          scale_color_manual(name=expression(R[c](t)),values=cols, 
                             guide = guide_legend(override.aes=aes(fill=NA)))+
          ylab("Reproduction Number") + xlab("Days")+theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 30)))
  ggsave(filename=paste0(plot_filepath,'SIM_control_Rt_',plot_index[i],'_blacked.tiff'),device='tiff',dpi=300,width=10,height=9)
  ggsave(filename=paste0(plot_filepath,'SIM_control_Rt_',plot_index[i],'_blacked.jpeg'),device='jpeg',dpi=300,width=10,height=9)
  
}

