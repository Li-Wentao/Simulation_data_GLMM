# This R code is written to generate simulated data for distributed trivairate GLMM
setwd("/Users/jiayito/Dropbox/000_UPenn_Research/000_project/000_xiaoqian/new/simulated_data")

library(MASS)

####################################
#    Values of parameters
####################################

expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(p){log(p/(1-p))}

# coefficients of all of the variables (log odds ratio); length = 10; 1 intercept + 9 variables
beta0 = c(-1.5,0.1,-0.5,-0.3,0.4,-0.2,-0.25,0.35,-0.1,0.5) # from -0.5 to 0.5
# true value of odds ratio: exp(beta0)

####################################
#     Data generation function
####################################

data.generator <- function(num_sites, site_size, sd_mu, seed, beta){
  
  set.seed(seed)
  
  # true sensitivity
  sens.true = 0.6 
  
  # true specificity
  spec.true = 0.9 
  
  # total number of observations
  N = num_sites * site_size
  
  # variables: X (indepedent from each other)
  x1 = rep(1, N) # fixed intercept
  x2 =  rbinom(N,1,prob=0.1)
  x3 =  rbinom(N,1,prob=0.3)
  x4 =  rbinom(N,1,prob=0.5)
  x5 =  rnorm(N, 0, 0.5)
  x6 =  rnorm(N, 0, 1)
  x7 =  rnorm(N, 0, 1.5)
  x8 =  runif(N, min = -0.5, max = 0.5)
  x9 =  runif(N, min = -0.7, max = 0.7)
  x10 =  runif(N, min = -1, max = 1)
  x = cbind(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10)
  
  # random error
  rand_err = rnorm(N, 0, 1)
  
  # create the ID
  ID = rep(1:num_sites, each = site_size)
  
  # generate random intercept using trivariate normal distribution (covariance matrix of mu_i1, mu_i2, mu_i3)
  Sigma = matrix(c(sd_mu, 0, 0,  
                   0, sd_mu, 0, 
                   0, 0, sd_mu), 
                 nrow=3, byrow=TRUE)
  mymean = c(0,0,0) # mu_i0, mu_i1, mu_i2, which follows trivariate normal distribution with mean (0,0,0)
  myRE = mvrnorm(n = num_sites, mymean, Sigma)  
  
  # g(\pi)= logit(\pi)
  logitPi =  x %*% beta + rand_err + rep(myRE[, 1],each=site_size)
  
  # g(se)
  beta1 = qlogis(sens.true)
  logitSe = beta1 + myRE[,2] 
  
  # g(sp)= logit(\pi)
  beta2 = qlogis(spec.true)
  logitSp = beta2 + myRE[,3]
  
  
  # D: true status
  D = rbinom(N, 1, prob = plogis(logitPi))
  
  # number of patients with outcome = 1 and 0
  PiY = rowSums(matrix(D, ncol = site_size, nrow = num_sites, byrow = TRUE)) 
  
  # number of patients whose test result is 1 conditioning on positive test result
  SeY = rbinom(num_sites, PiY, prob=plogis(logitSe)) 
  
  # number of patients whose test result is 0 conditioning on negative test result
  SpY = rbinom(num_sites, (site_size-PiY), prob=plogis(logitSp)) 
  
  # number of disease patient in each hospital, a vector with length num_sites
  Nd = PiY 
  
  # number of not disease patient in each hospital, a vector with length num_sites
  Nn = site_size-PiY 
  
  # my dataset
  mydat = data.frame(cbind(ID, # ID
                           rep(site_size, each = num_sites*site_size), # site size
                           logitPi, # logit Pi
                           rep(Nd, each = site_size), # true disease
                           rep(SeY, each = site_size), # test disease
                           rep(Nn, each = site_size), # true not disease
                           rep(SpY, each = site_size), # test not disease
                           x)) # variables
  names(mydat) = c("Site_ID", 
                   "Site_sample_size",
                   "LogitPi", 
                   "Number_of_true_disease",
                   "Number_of_test_positive_among_true_disease",
                   "Number_of_true_not_disease",
                   "Number_of_test_negative_among_true_not_disease",
                   "X1", "X2", "X3", "X4", "X5", "X6", "X7", "X8", "X9","X10")
  
  return(mydat)
}



####################################
#     Data generation: under each setting, we generated 20 datasets.
####################################

# generate 20 seeds, which are used to generate the datasets
set.seed(1234)
seed_list = sample(1:10000, 20)
seed_list

####################################
#        Setting One: 
#    number of sites: 2
#    number of variables: 10
#    variance of intercept: small
####################################

for(i in 1:20){
  dataset1 = data.generator(num_sites = 2, site_size = 500, sd_mu = 1, seed = seed_list[i], beta = beta0)
  write.csv(dataset1, paste("./Setting_1/2Sites_500PatientsEachSite_SmallVar1_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Two: 
#    number of sites: 2
#    number of variables: 10
#    variance of intercept: large
####################################


for(i in 1:20){
  dataset2 = data.generator(num_sites = 2, site_size = 500, sd_mu = 2, seed = seed_list[i], beta = beta0)
  write.csv(dataset2, paste("./Setting_2/2Sites_500PatientsEachSite_LargeVar4_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Three: 
#    number of sites: 10
#    number of variables: 10
#    variance of intercept: small
####################################

for(i in 1:20){
  dataset3 = data.generator(num_sites = 10, site_size = 500, sd_mu = 1, seed = seed_list[i], beta = beta0)
  write.csv(dataset3, paste("./Setting_3/10Sites_500PatientsEachSite_SmallVar1_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Four: 
#    number of sites: 10
#    number of variables: 10
#    variance of intercept: large
####################################

for(i in 1:20){
  dataset4 = data.generator(num_sites = 10, site_size = 500, sd_mu = 2, seed = seed_list[i], beta = beta0)
  write.csv(dataset4, paste("./Setting_4/10Sites_500PatientsEachSite_LargeVar4_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Five: 
#    number of sites: 2
#    number of variables: 10
#    variance of intercept: small
####################################

for(i in 1:20){
  dataset5 = data.generator(num_sites = 2, site_size = 30, sd_mu = 1, seed = seed_list[i], beta = beta0)
  write.csv(dataset5, paste("./Setting_5/2Sites_30PatientsEachSite_SmallVar1_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Six: 
#    number of sites: 2
#    number of variables: 10
#    variance of intercept: large
####################################


for(i in 1:20){
  dataset6 = data.generator(num_sites = 2, site_size = 30, sd_mu = 2, seed = seed_list[i], beta = beta0)
  write.csv(dataset6, paste("./Setting_6/2Sites_30PatientsEachSite_LargeVar4_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Seven: 
#    number of sites: 10
#    number of variables: 10
#    variance of intercept: small
####################################

for(i in 1:20){
  dataset7 = data.generator(num_sites = 10, site_size = 30, sd_mu = 1, seed = seed_list[i], beta = beta0)
  write.csv(dataset7, paste("./Setting_7/10Sites_30PatientsEachSite_SmallVar1_Dataset", i, ".csv", sep = ""))
}



####################################
#        Setting Eight: 
#    number of sites: 10
#    number of variables: 10
#    variance of intercept: large
####################################

for(i in 1:20){
  dataset8 = data.generator(num_sites = 10, site_size = 30, sd_mu = 2, seed = seed_list[i], beta = beta0)
  write.csv(dataset8, paste("./Setting_8/10Sites_30PatientsEachSite_LargeVar4_Dataset", i, ".csv", sep = ""))
}



