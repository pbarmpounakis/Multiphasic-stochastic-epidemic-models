data{  
  int<lower=0> N0; // number of days to impute
  
  int<lower=1> Nobs; // number of days observed must be <= N
    
  int<lower=0> deaths[Nobs]; //daily deaths 

  int<lower=0> pop; //population
   
  int<lower=0> reported_cases[Nobs]; //daily reported cases 

  vector<lower=0>[Nobs] s_int_rev1 ; // serial interval reverse order vector for generation distribution
  
  vector<lower=0>[Nobs] s_int_rev2 ; // serial interval reverse order vector for infection to death distribution
      
  real<lower=0,upper=1> ifr[4]; // infection-fatality ratio
  
  int<lower=0> ifr_cp[3];// changepoint for ifr

  int<lower=0> epi_start;// day to start for deaths
}



parameters{
real<lower=0> R0;
vector[6] effect;
vector<lower=0>[N0] prediction0;
real<lower=0> phi;
real<lower=0> tau;
real<lower=3,upper=Nobs-3> T1;
real<lower=0,upper=100> e[5]; //extra time for the changepoint number 2, 3, 4, 5, 6

}

transformed parameters{

real<lower=T1,upper=Nobs-3> T2;
real<lower=T2,upper=Nobs-3> T3;
real<lower=T3,upper=Nobs-3> T4;
real<lower=T4,upper=Nobs-3> T5;
real<lower=T5,upper=Nobs-3> T6;

vector<lower=0>[Nobs] Rt ;
vector<lower=0>[Nobs] Rt_adj ;
vector<lower=0>[Nobs] cumul_sum ;
vector<lower=0>[Nobs] prediction;
vector<lower=0>[Nobs] E_deaths;


T2=T1+e[1];
T3=T2+e[2];
T4=T3+e[3];
T5=T4+e[4];
T6=T5+e[5];



prediction[1:N0]=prediction0;
cumul_sum[1:N0]=cumulative_sum(prediction[1:N0]);

for( t in 1:Nobs){

Rt[t]=R0*exp(-inv_logit(5* (t-T1))*inv_logit(5* (T2-t))*effect[1]-inv_logit(5* (t-T2))*inv_logit(5* (T3-t))*effect[2]-inv_logit(5* (t-T3))*inv_logit(5* (T4-t))*effect[3]-inv_logit(5*(t-T4))*inv_logit(5*(T5-t))*effect[4]-inv_logit(5*(t-T5))*inv_logit(5*(T6-t))*effect[5]-inv_logit(5*(t-T6))*effect[6]);
}
 
Rt_adj[1:N0+1]=Rt[1:N0+1];


for(i in (N0+1):Nobs-1){
real convolution=dot_product(prediction[1:i-1] ,s_int_rev1[Nobs-i+2:Nobs]);
prediction[i]=Rt_adj[i]*convolution;
cumul_sum[i]=cumul_sum[i-1]+prediction[i];
Rt_adj[i+1]=((pop-cumul_sum[i])/pop)*Rt[i+1];

}
prediction[Nobs]=Rt_adj[Nobs]*dot_product(prediction[1:Nobs-1] ,s_int_rev1[2:Nobs]);
cumul_sum[Nobs]=cumul_sum[Nobs-1]+prediction[Nobs];

E_deaths[1]=ifr[1]*prediction[1]*s_int_rev2[Nobs];
for(i in 2:Nobs){
if (i<=ifr_cp[1])
E_deaths[i]=ifr[1]*dot_product(prediction[1:i-1] ,s_int_rev2[Nobs-i+2:Nobs]);
else if (i>ifr_cp[1] && i<=ifr_cp[2])
E_deaths[i]=ifr[2]*dot_product(prediction[1:i-1] ,s_int_rev2[Nobs-i+2:Nobs]);
else if (i>ifr_cp[2] && i<=ifr_cp[3])
E_deaths[i]=ifr[3]*dot_product(prediction[1:i-1] ,s_int_rev2[Nobs-i+2:Nobs]);
else if (i>ifr_cp[3])
E_deaths[i]=ifr[4]*dot_product(prediction[1:i-1] ,s_int_rev2[Nobs-i+2:Nobs]);
}


}

model{
R0 ~ normal(3,1);
effect ~ normal(0,3);
tau ~ exponential(0.03);
prediction0 ~ exponential(1/tau);
phi ~ cauchy(0,5);
deaths[epi_start:Nobs] ~ neg_binomial_2(E_deaths[epi_start:Nobs], phi);

}
generated quantities{
real Deviance;
vector[Nobs-epi_start+1]  log_likelihood;
vector[Nobs] rep_true_ratio;

for(i in 1:Nobs){
rep_true_ratio[i]= reported_cases[i] / prediction[i];
}


for( i in 1:(Nobs-epi_start+1)){
log_likelihood[i]= neg_binomial_2_lpmf(deaths[epi_start+i-1]|E_deaths[epi_start+i-1],phi);
}

Deviance=-2*sum(log_likelihood);
}
