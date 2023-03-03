library(nimble)
library(parallel)
library(R0)
library(dplyr)
library(bayesplot)
library(rstan)
library(ggplot2)
library(ggpubr)

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

constants <- list(N = length(cases),Trunc=18,gen_interval=gen_interval,len_gen_int=length(gen_interval),pop=pop)


data=list(cases=cases)
run_MCMC_NB <- function(seed,data,constants,nburnin,niter,thin) {
  library(nimble)
  inits<-  list(R = runif(constants$Trunc,0,3),
                v = rbeta(constants$Trunc-1, 1, 1),
                alpha = runif(1,0,10),
                Rt=runif(length(data$cases),0,3),
                z=sample(1:constants$Trunc,length(data$cases),replace = T)
  )  
  
    DPcode<- nimbleCode(
    {
      
      for(i in (len_gen_int+1):N) {
        m[i]<-Rt[i]*((pop-sum(cases[1:(i-1)]))/pop)*sum(cases[(i-len_gen_int):(i-1)]*gen_interval[1:len_gen_int])
        p[i] <- r/(m[i] + r)
        cases[i] ~ dnegbin(p[i],r)
      }
      r ~ dcat(pi[1:1000])
      for (i in 1:1000) {pi[i] <- 1/1000}
      
      
      for(i in 1:(N)) {
        z[i] ~ dcat(w[1:Trunc])
        Rt[i]<-R[z[i]]
      }
      for(i in 1:Trunc) {
        R[i] ~ T(dnorm(0, sd = 3),0,)
      }
      for(i in 1:(Trunc-1)) { # stick-breaking variables
        v[i] ~ dbeta(1, alpha)
      }
      w[1:Trunc] <- stick_breaking(v[1:(Trunc-1)]) # stick-breaking weights
      alpha~dgamma(1,1)
      
    }
  )
  
  
  DPmodel <- nimbleModel(DPcode, data = data,inits = inits, constants = constants) 
  
  CDPmodel_MCMCConf <- configureMCMC(DPmodel, print = TRUE,enableWAIC = TRUE)
  
  
  CDPmodel_MCMCConf$addMonitors(c('Rt','z','w'))
  
  DPmodel_MCMC<- buildMCMC(CDPmodel_MCMCConf)
  
  CDPmodel<- compileNimble(DPmodel,showCompilerOutput = TRUE)
  
  CDPmodel_MCMC <- compileNimble(DPmodel_MCMC, project = DPmodel)
  
  CDPmodel_samples <- runMCMC(CDPmodel_MCMC, nburnin = nburnin, niter = niter,thin=thin, setSeed = seed,summary=T,WAIC=T)
  
  return(CDPmodel_samples)
}

this_cluster <- makeCluster(8)

start_time<-Sys.time()
  
run<-parLapply(cl = this_cluster, X = seq(1,80,10), 
               fun = run_MCMC_NB,
               data = data,constants=constants,nburnin = 50000, niter =100000,thin=10)
end_time<-Sys.time()
stopCluster(this_cluster)  


save(run,file='C:\\Users\\barbounakis\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\Nimble\\Final Renewal equations - NB likelihood easy sim\\DP\\stick\\DP_NB_sim.RData')


#plot for Rt
Rt_real<-c(rep(1.5,60),rep(0.95,40),rep(1.35,50),rep(0.8,50),rep(1.8,50))#cp=c(60,100,150,200),R=c(1.5,0.95,1.35,0.8,1.8)

get_results<-function(chain_output,chains){
  
  nrows<-dim(chain_output[[1]]$samples)[1]
  ncols<-dim(chain_output[[1]]$samples)[2]
  parameters<-matrix(NA,nrow = length(chains)*nrows, ncol=ncols)
  
  for(j in 1:length(chains)){
    parameters[(1+(j-1)*(nrows)):(j*nrows),]<-chain_output[[chains[j]]]$samples
  }
  colnames(parameters)<-colnames(chain_output[[1]]$samples)
  
  
  mean_vector<-apply(parameters,2,mean)
  sd_vector<-apply(parameters,2,sd)
  median_vector<-apply(parameters,2,quantile,probs=0.5,na.rm=T)
  
  quant_low<-apply(parameters,2,quantile,probs=0.025,na.rm=T)
  quant_high<-apply(parameters,2,quantile,probs=0.975,na.rm=T)
  ess_bulk_vector<-apply(parameters,2,ess_bulk)
  
  ess_tail_vector<-apply(parameters,2,ess_tail)
  Rhat_vector<-apply(parameters,2,Rhat)
  
  results<-round(cbind(mean_vector,sd_vector,median_vector,quant_low,quant_high,ess_bulk_vector,ess_tail_vector,Rhat_vector),2)
  return(results)
}

Rt_results_all_chains<-get_results(run,chains=1:8)




#dates<-1:250
Rt_real<-c(rep(1.5,60),rep(0.95,40),rep(1.35,50),rep(0.8,50),rep(1.8,50))#cp=c(60,100,150,200),R=c(1.5,0.95,1.35,0.8,1.8)
dates<-1:250
cases_daily<-cases
#plot for Rt-NB Likelihood
length(cases_daily)
starting_day<-1
dat_Rt<-list()
grph<-list()
for(j in 1:8){
  dat_Rt[[j]]<-data.frame(real_Rt=Rt_real,est_Rt=run[[j]]$summary[(18+starting_day):(18+length(cases_daily)),2],est_Rt_low=run[[j]]$summary[(18+starting_day):(18+length(cases_daily)),4],est_Rt_high=run[[j]]$summary[(18+starting_day):(18+length(cases_daily)),5])
  cols <- c("Median estimate with 95% Cr.I."="red",'Real'='black')
  grph[[j]]<-ggplot() +
    # points
    geom_line(data=dat_Rt[[j]], aes(x=dates[(starting_day):length(cases_daily)], y=real_Rt, colour="Real"),size=1) + 
    # red line
    geom_line(data=dat_Rt[[j]], aes(x=dates[(starting_day):length(cases_daily)], y=est_Rt,colour="Median estimate with 95% Cr.I."),
              size=1)+
    geom_ribbon(aes(ymin = est_Rt_low, ymax = est_Rt_high,x=dates[(starting_day):length(cases_daily)]),data=dat_Rt[[j]], alpha = 0.5, 
                linetype="dashed",
                colour="red",fill='red')+ geom_hline(yintercept=1, linetype="dashed", color = "red",size=1)+
    scale_color_manual(name=expression(R[c](t)),values=cols, 
                       guide = guide_legend(override.aes=aes(fill=NA)))+
    ylab("") + xlab("Days")+ylim(c(0,3)) +theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 30))
  
  
  
}


print(ggarrange(grph[[1]], grph[[2]], grph[[3]], grph[[4]],   
                labels = paste0('chain ',1:4),label.x = 0.1,
                label.y = 1,legend='none',common.legend = TRUE,
                ncol = 2, nrow = 2))
plot_filepath<-"C:\\Users\\pbarb\\OneDrive - aueb.gr\\PhD\\Coding\\Covid-19\\Branching process\\Nimble\\Final Renewal equations - NB likelihood easy sim\\DP\\Stick\\Deaths & cases 2 frameworks FINAL - 2\\plots Rt estimations\\simulation\\"
ggsave(filename=paste0(plot_filepath,'Daily_infections_DP_multiple_chains.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'Daily_infections_DP_multiple_chains.jpeg'),device='jpeg',dpi=300,width=10,height=9)



#plot for  control Rt-NB Likelihood
length(cases_daily)
starting_day<-1
dat_Rt<-data.frame(real_Rt=Rt_real,est_Rt=Rt_results_all_chains[(18+starting_day):(18+length(cases_daily)),3],est_Rt_low=Rt_results_all_chains[(18+starting_day):(18+length(cases_daily)),4],est_Rt_high=Rt_results_all_chains[(18+starting_day):(18+length(cases_daily)),5])
cols <- c("Median estimate with 95% Cr.I."="black",'Real'='black')
print(ggplot() +
        # points
        geom_line(data=dat_Rt, aes(x=dates[(starting_day):length(cases_daily)], y=real_Rt, colour="Real" ),size=1) + 
        #red line
        geom_line(data=dat_Rt, aes(x=dates[(starting_day):length(cases_daily)], y=est_Rt,colour="Median estimate with 95% Cr.I."),
                  size=0.5,linetype="dashed")+
        geom_ribbon(aes(ymin = est_Rt_low, ymax = est_Rt_high,x=dates[(starting_day):length(cases_daily)]),data=dat_Rt, alpha = 0.5, 
                    linetype="dashed",
                    colour="black",fill='black')+ geom_hline(yintercept=1, linetype="dashed", color = "black",size=0.5)+
        scale_color_manual(name=expression(R[c](t)),values=cols, 
                           guide = guide_legend(override.aes=aes(fill=NA)))+
        ylab("Reproduction Number") + xlab("Days")+ylim(c(0,3)) +theme(axis.text.x=element_text(angle=60, hjust=1))+ theme(legend.position = "none")+ theme(text = element_text(size = 30))
)
ggsave(filename=paste0(plot_filepath,'Daily_infections_DP_Rt_chains_combined_blacked.tiff'),device='tiff',dpi=300,width=10,height=9)
ggsave(filename=paste0(plot_filepath,'Daily_infections_DP_Rt_chains_combined_blacked.jpeg'),device='jpeg',dpi=300,width=10,height=9)


