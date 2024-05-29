library(nimble)
library(parallel)
library(R0)
library(dplyr)
library(bayesplot)
#library(rstan)
library(ggplot2)
library(splines2)
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


data<-list(deaths=death_cases)
constants<- list(N = length(death_cases),death_interval=death_interval,len_death_int=length(death_interval),ifr=ifr)
inits <- list(cases=rpois(length(data$deaths),1000)  )  


start_time_global<-Sys.time()

DPcode_poisson <- nimbleCode(
  { 
    
    for(i in (len_death_int+1):N){
      mu_deaths[i]<-ifr*sum(cases[(i-len_death_int):(i-1)]*death_interval[1:len_death_int])
      deaths[i] ~ dpois(mu_deaths[i])
    }
    for(i in 1:30) {
      cases[i] ~ dpois(lambda)

    }
    lambda~dgamma(1,0.000001)  
    
    for(i in (31):N) {
      cases[i] ~ T(dnorm(cases[i-1],sd=sqrt(sigma2)),0,)
    }
    sigma2~dinvgamma((1/2)*nu,1/2) 
    nu~dgamma(0.2,0.1)
  }
)


DPmodel_poisson <- nimbleModel(DPcode_poisson, data = data,inits = inits, constants = constants) 
  
CDPmodel_MCMCConf_poisson <- configureMCMC(DPmodel_poisson, print = TRUE,enableWAIC = TRUE)
  
  
CDPmodel_MCMCConf_poisson$addMonitors(c('cases','mu_deaths','sigma2','nu'))
  
DPmodel_MCMC_poisson <- buildMCMC(CDPmodel_MCMCConf_poisson)
  
CDPmodel_poisson<- compileNimble(DPmodel_poisson,showCompilerOutput = TRUE)
  
CDPmodel_MCMC_poisson <- compileNimble(DPmodel_MCMC_poisson, project = DPmodel_poisson)
  
CDPmodel_samples_poisson<- runMCMC(CDPmodel_MCMC_poisson, nburnin = 50000, niter = 100000,thin=50, setSeed = 1,summary=T,WAIC=T)
  



end_time<-Sys.time()
print(end_time-start_time_global)



save(CDPmodel_samples_poisson,file='C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\Nimble\\Final Renewal equations - NB likelihood easy sim\\Estimation of cases\\Final simple model alone, similar timeframe with the other models\\sim_simple_model_student.RData')

est_cases_median<-CDPmodel_samples_poisson$summary[1:250,2]
est_cases_low<-CDPmodel_samples_poisson$summary[1:250,4]
est_cases_high<-CDPmodel_samples_poisson$summary[1:250,5]

dat_cases<-data.frame(real_cases=cases,est_cases=est_cases_median,est_cases_low=est_cases_low,est_cases_high=est_cases_high)
cols <- c("Median with 95% Cr.I."="red")
shapes <- c(" "=17)

ggplot() +
  # points
  geom_point(data=dat_cases, aes(x=1:nrow(dat_cases), y=real_cases,shape=' '), colour="black") + 
  # red line
  geom_line(data=dat_cases, aes(x=1:nrow(dat_cases), y=est_cases,colour="Median with 95% Cr.I."),
            size=1)+
  geom_ribbon(aes(ymin = est_cases_low, ymax = est_cases_high,x=1:nrow(dat_cases)),data=dat_cases, 
              linetype="dashed",
              color="red",alpha=0.25,fill='red')+
  scale_shape_manual(name="Real Cases",values=shapes, 
                     guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_color_manual(name="Estimated Cases",values=cols, 
                     guide = guide_legend(override.aes=aes(fill=NA)))+
  ylab("Daily infections") + xlab("Days")   +theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 25)) 

plot_filepath<-"C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\Nimble\\Final Renewal equations - NB likelihood easy sim\\Estimation of cases\\Final simple model alone, similar timeframe with the other models\\plots\\"
plot_index<-c('df45','df35','df25','df15')

ggsave(filename=paste0(plot_filepath,'SIM_cases_95CrI_nosm_st_t.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_cases_95CrI_nosm_st_t.jpeg'),device='jpeg',dpi=300,width=10,height=9)


#Deaths fits
est_deaths_median<-CDPmodel_samples_poisson$summary[(1*250+2):(2*250+1),2]
est_deaths_low<-CDPmodel_samples_poisson$summary[(1*250+2):(2*250+1),4]
est_deaths_high<-CDPmodel_samples_poisson$summary[(1*250+2):(2*250+1),5]



dat_deaths<-data.frame(real_deaths=death_cases,est_deaths=est_deaths_median,est_deaths_low=est_deaths_low,est_deaths_high=est_deaths_high)
cols <- c("Median with 95% Cr.I."="red")
shapes <- c(" "=17)

ggplot() +
  # points
  geom_point(data=dat_deaths, aes(x=1:nrow(dat_deaths), y=real_deaths,shape=' '), colour="black") + 
  # red line
  geom_line(data=dat_deaths, aes(x=1:nrow(dat_deaths), y=est_deaths,colour="Median with 95% Cr.I."),
            size=1)+
  geom_ribbon(aes(ymin = est_deaths_low, ymax = est_deaths_high,x=1:nrow(dat_deaths)),data=dat_deaths, alpha = 0.1, 
              linetype="dashed",
              color="red",alpha=0.5,fill='red')+
  scale_shape_manual(name="Real Deaths",values=shapes, 
                     guide = guide_legend(override.aes=aes(fill=NA)))+
  scale_color_manual(name="Estimated Deaths",values=cols, 
                     guide = guide_legend(override.aes=aes(fill=NA)))+
  ylab("Daily deaths") + xlab("Days") +theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 25)) 

ggsave(filename=paste0(plot_filepath,'SIM_deaths_95CrI_student_t.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'SIM_deaths_95CrI_student_t.jpeg'),device='jpeg',dpi=300,width=10,height=9)

#spline smooth
level_smoothing<-c(45,35,25,17)
for(i in 1:length(level_smoothing)){
  sm_cases_median<-smooth.spline(est_cases_median,df=level_smoothing[i])
  sm_cases_low<-smooth.spline(est_cases_low,df=level_smoothing[i])
  sm_cases_high<-smooth.spline(est_cases_high,df=level_smoothing[i])
  

  dat_cases<-data.frame(real_cases=cases,est_cases=sm_cases_median$y,est_cases_low=sm_cases_low$y,est_cases_high=sm_cases_high$y)
  cols <- c("Median with 95% Cr.I."="red")
  shapes <- c(" "=17)
  
  print(ggplot() +
    # points
    geom_point(data=dat_cases, aes(x=1:nrow(dat_cases), y=real_cases,shape=' '), colour="black") + 
    # red line
    geom_line(data=dat_cases, aes(x=1:nrow(dat_cases), y=est_cases,colour="Median with 95% Cr.I."),
              size=1)+
    geom_ribbon(aes(ymin = est_cases_low, ymax = est_cases_high,x=1:nrow(dat_cases)),data=dat_cases, 
                linetype="dashed",
                color="red",alpha=0.25, fill = "red")+
    scale_shape_manual(name="Real Cases",values=shapes, 
                       guide = guide_legend(override.aes=aes(fill=NA)))+
    scale_color_manual(name="Estimated Cases",values=cols, 
                       guide = guide_legend(override.aes=aes(fill=NA)))+
    ylab("Daily infections") + xlab("Days" )   +theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 25)) )

  ggsave(filename=paste0(plot_filepath,'SIM_cases_95CrI_st_t',plot_index[i],'.tiff'),device='tiff',dpi=300,width=10,height=9)
  ggsave(filename=paste0(plot_filepath,'SIM_cases_95CrI_st_t',plot_index[i],'.jpeg'),device='jpeg',dpi=300,width=10,height=9)
  
}

