rm(list=ls())
set.seed(12345)

N=10000 # population size per scenario

##Scenario 1
###Low SES inequality+week SES health relation
s1=rlnorm(N, meanlog=0, sdlog=0.30)
error_s1=rnorm(N,mean=0,sd=1)
health_s1=5+0.20*s1+error_s1
health_s1=health_s1-min(health_s1)+0.01

##Scenario 2
###Moderate SES ineqaulity+moderate SES health relation
s2=rlnorm(N, meanlog=0,sdlog=0.80)
error_s2=rnorm(N,mean=0,sd=1)
health_s2=5+0.50*s2+error_s2
health_s2=health_s2-min(health_s2)+0.01

## Scenario 3
###High SES inequality/heavy upper tail
u_s3=runif(N)
xm=1
alpha=2
s3=xm/(u_s3^(1/alpha))
error_s3=rnorm(N,mean=0,sd=1)
health_s3=5+0.80*s3+error_s3
health_s3=health_s3-min(health_s3)+0.01

## Scenario s4
###bimodal society
group_s4=rbinom(N,size=1,prob=0.50)
ses_s4=rep(NA,N)
ses_s4[group_s4==0]=rlnorm(sum(group_s4==0),
				meanlog=0,sdlog=0.30)
ses_s4[group_s4==1]=rlnorm(sum(group_s4==1),
				meanlog=2,sdlog=0.50)
error_s4=rnorm(N,mean=0,sd=1)
health_s4=5+0.50*ses_s4+error_s4
health_s4=health_s4-min(health_s4)+0.01

##Scenario 5
### Nonlinear logramithmic relation
s5=rlnorm(N,meanlog=0,sdlog=0.80)
error_s5=rnorm(N,mean=0,sd=1)
health_s5=5+1.00*log(s5)+error_s5
health_s5=health_s5-min(health_s5)+0.01

## Scenario 6
###Threshold relation
s6=rlnorm(N,meanlog=0,sdlog=0.80)
error_s6=rnorm(N,mean=0,sd=1)
cutoff_s6=median(s6)
health_s6=5+1.50*ifelse(s6>cutoff_s6,1,0)+error_s6
health_s6=health_s6-min(health_s6)+0.01

##Scenario 7
###Bounded proportion health outcome
s7=rlnorm(N,meanlog=0,sdlog=0.80)
eta_s7=-1+0.80*as.numeric(scale(log(s7)))
mu_s7=exp(eta_s7)/(1+exp(eta_s7))
phi=20
health_s7=rbeta(N,
			shape1=mu_s7*phi,
			shape2=(1-mu_s7)*phi)

## Scneario 8
###Binary health outcome
s8=rlnorm(N,meanlog=0,sdlog=0.80)
eta_s8=-1+0.80*as.numeric(scale(log(s8)))
p_s8=exp(eta_s8)/(1+exp(eta_s8))
health_s8=rbinom(N,size=1,prob=p_s8)

##Relative SES variable
rel_ses_s1 = s1 / mean(s1)
rel_ses_s2 = s2 / mean(s2)
rel_ses_s3 = s3 / mean(s3)
rel_ses_s4 = ses_s4 / mean(ses_s4)
rel_ses_s5 = s5 / mean(s5)
rel_ses_s6 = s6 / mean(s6)
rel_ses_s7 = s7 / mean(s7)
rel_ses_s8 = s8 / mean(s8)

## Final population data
library(data.table)
pop = data.table(
  id = 1:N,
  ses_s1 = s1,
  rel_ses_s1 = rel_ses_s1,
  health_s1 = health_s1,
  ses_s2 = s2,
  rel_ses_s2 = rel_ses_s2,
  health_s2 = health_s2,
  ses_s3 = s3,
  rel_ses_s3 = rel_ses_s3,
  health_s3 = health_s3,
  ses_s4 = ses_s4,
  rel_ses_s4 = rel_ses_s4,
  health_s4 = health_s4,
  ses_s5 = s5,
  rel_ses_s5 = rel_ses_s5,
  health_s5 = health_s5,
  ses_s6 = s6,
  rel_ses_s6 = rel_ses_s6,
  health_s6 = health_s6,
  ses_s7 = s7,
  rel_ses_s7 = rel_ses_s7,
  health_s7 = health_s7,
  ses_s8 = s8,
  rel_ses_s8 = rel_ses_s8,
  health_s8 = health_s8
)

## Check data
head(pop)
summary(pop)
colnames(pop)

## Save data
save_path <- "D:/Lab/SES-Related Health Inequality/final SES paper/R files/SES_simulation_population.csv"
fwrite(pop, save_path)

cat("Data saved successfully.\n")
cat("File location:\n")
cat(save_path, "\n")
