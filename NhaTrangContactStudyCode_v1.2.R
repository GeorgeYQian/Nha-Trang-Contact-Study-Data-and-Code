######################################################################################################
# This file contains the code required to obtain the results shown in the paper
# (the cubic p-spline fit to age and Bayesian logistic regression)
#
# Please address any questions to me (george): george.qian08@alumni.imperial.ac.uk
#
# G Qian, 2021 on behalf of the authors of:
# "Pneumococcal exposure routes for infants, a nested cross-sectional survey in Nha Trang, Vietnam"
######################################################################################################


###################################
### load libraries
###################################
#Import libraries
library(tidyverse)
library(rjags)
library(binom)
library(runjags)
library(mmcc)
library(splines)
library(reshape2)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)


###################################
### load data
###################################
#Load data from the contact study
df_contacts <- read.csv("df_contacts_public.csv")

#Define age grouping and arrange infant contacts into these groups
#Here, age groups are simply 1 year 'bands' (e.g. 1 yr old, 2 yrs etc)
brks <- 0:100
df_contacts$ctagegp <- cut(df_contacts$contact_age, brks, right = F)

#Load relevant data from the carriage study
df_infants <- read.csv("df_infants_public.csv")

#The Carriage Data here uses the data from Carla's Thesis.
#The following rates were calculated using Table 6.1 in Carla Talarico's thesis (page 179),
#available online: https://deepblue.lib.umich.edu/handle/2027.42/64706)
CarriageDataThesis3 <- read.csv("CarriageDataCarlaHalfYear.csv")

# Carriage prevalence observed in the trial
RateCarInf46months <- 0.1478873
RateCarInf68months <- 0.1839378
RateCarInf810months <- 0.2358722
RateCarInf1012months <- 0.27897
RateCarTod1214months <- 0.2769231
RateCarTod1416months <- 0.352381
RateCarTod1618months <- 0.3449477
RateCarTod1820months <- 0.3556851
RateCarTod2022months <- 0.4134276
RateCarTod2224months <- 0.4205298
carriageR <-
  c(
    RateCarInf46months,
    RateCarInf68months,
    RateCarInf810months,
    RateCarInf1012months,
    RateCarTod1214months,
    RateCarTod1416months,
    RateCarTod1618months,
    RateCarTod1820months,
    RateCarTod2022months,
    RateCarTod2224months
  ) #read.csv("CarlaDataCarriageRate.csv")
age <-  c(5 / 12, 7 / 12, 9 / 12, 11 / 12, 13 / 12, 15 / 12, 17 / 12, 19 /
            12, 21 / 12, 23 / 12)
NTrang4thUnder2 <- data.frame(age, carriageR)

# some data wrangling to get individual level data into format needed by model
CarriageDataThesis3 <- CarriageDataThesis3 %>%
  rename(CarriageStatus = CarriageStatusH,
         AllAges = AllAgesH)

CarriageDataThesis3_ <- CarriageDataThesis3 %>%
  count(AllAges, CarriageStatus) %>%
  pivot_wider(names_from = CarriageStatus,
              values_from = n,
              values_fill = 0) %>%
  mutate(N = `0` + `1`) %>%
  mutate(y = `1`,
         p = y / N) %>%
  nest(data = c(y, N)) %>%
  mutate(CI = map(data, ~ prop.test(.x$y, .x$N)$conf.int)) %>%
  mutate(CI = map(CI, ~ data.frame(
    value = .x, key = c("L", "U")
  ))) %>%
  mutate(CI = map(CI, ~ pivot_wider(
    .x, values_from  = value,
    names_from = key
  ))) %>%
  unnest(CI)

###################################
### First Jags Model (cubic spline)
###################################

bspline <-
  function(x,
           K,
           bdeg = 3,
           cyclic = FALSE,
           xl = min(x),
           xr = max(x)) {
    x <- as.matrix(x, ncol = 1)
    
    ndx <- K - bdeg
    
    # p-spline algorithm, as outlined in Eilers and Marx (1996)
    dx <- (xr - xl) / ndx
    t <- xl + dx * (-bdeg:(ndx + bdeg))
    T <- (0 * x + 1) %*% t
    X <- x %*% (0 * t + 1)
    P <- (X - T) / dx
    B <- (T <= X) & (X < (T + dx))
    r <- c(2:length(t), 1)
    
    for (k in 1:bdeg) {
      B <- (P * B + (k + 1 - P) * B[, r]) / k
      
    }
    
    B <- B[, 1:(ndx + bdeg)]
    
    if (cyclic == 1) {
      for (i in 1:bdeg) {
        B[, i] <- B[, i] + B[, K - bdeg + i]
      }
      B <- B[, 1:(K - bdeg)]
    }
    
    return(B)
  }

K = 6
knots <- c(0, 3, 6, 9, 12, 18, 24)

X <- bs(
  x = unlist(CarriageDataThesis3_$AllAges),
  knots = knots,
  intercept = FALSE
)

x.pred <- seq(0, 99, length.out = 100)

X.pred <- bs(x = x.pred,
             knots = knots,
             intercept = FALSE)

makeQ = function(degree, K, epsilon = 1e-3) {
  x <- diag(K)
  E <- diff(x, differences = degree)
  return(t(E) %*% E + x * epsilon)
}

Q <- makeQ(2, nrow(X))

attributes(X) <- attributes(X)["dim"]
attributes(X.pred) <- attributes(X.pred)["dim"]

#Prepare 
dat <- CarriageDataThesis3_ %>%
  unnest(data) %>%
  select(y, N) %>% as.list %>%
  append(.,
         values = list(
           K = ncol(X),
           Q = Q,
           X.pred = X.pred,
           X      = X,
           n = nrow(X),
           m = nrow(X.pred),
           beta.0 = matrix(data = rep(0, ncol(X)),
                           ncol = 1)
         ))

#Implement RJags to find the cubic spline fit
mod <-
  jags.model(file = "spline_model.R",
             data = dat,
             n.adapt = 5e3)

number_iter = 1000
burn <- jags.samples(mod, variable.names = "p.rep", n.iter = number_iter)
pred <- coda.samples(mod, "p.rep", n.iter = number_iter)
post <- coda.samples(mod, c("beta", "beta.00"), n.iter = number_iter)

# model diagnostics
statsBeta <- tidy(post)
#plot(post)

#Now we can plot the spline
Figure2a <- tidy(pred) %>%
  bind_cols(data.frame(x = x.pred)) %>%
  ggplot(data = ., aes(x = x, y = mean)) +
  # geom_point(data= CarriageDataThesis3_,
  #            aes(x = AllAges, y= p),color="red") +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`),
              alpha = 0.5,
              fill = "lightskyblue") +
  geom_line() +
  geom_point(data = CarriageDataThesis3_,
             aes(x = AllAges, y = p),
             color = "red",
             size = 3) + geom_linerange(data = CarriageDataThesis3_, aes(
               x = AllAges,
               ymin = L,
               ymax = U,
               y = p
             )) +
  xlab("Age") +
  theme_bw() +
  ylab("Prevalence") +
  scale_x_continuous(trans = "sqrt", #we will use a square root scale on the x-axis to make the lower ages more prominent
                     breaks = c(1:10) ^ 2) +
  geom_point(data = data.frame(x = knots, mean = 0), color = "grey") + xlab("Age") +
  geom_point(data = NTrang4thUnder2,
             aes(age, carriageR),
             colour = "darkgreen",
             size = 3) + scale_fill_discrete(breaks = c("trt1", "ctrl", "trt2")) + labs(y = "Carriage Prevalence", x = "Age (years)")
Figure2a


#####################
### Define Data Block
#####################

#Now take the output of this spline and use them to calculate PEI.
#These will be inputs to the second JAGS model

#For each set of beta parameters in 'post':
# {Get the formula of the cubic spline
#  Calculate the carriage prevalence at all 101 ages (0-100)
#  Calculate PEI (only the prevalence changes across samples of 'post')
# }


pred_prevalence <- 
  pred %>%
  mcmc_to_dt() %>%
  pivot_wider(names_from = iteration, values_from = value, names_prefix = "post ") %>%
  mutate(ctagegp = cut(0:99,breaks = 0:100, right=F), .before="post 1") %>%
  select(-c(chain,parameter)) 

PEI <-
  df_contacts %>%
  tibble() %>%
  select(-c(X, contact_serial_no, contact_age)) %>%
  left_join(pred_prevalence, by = "ctagegp") %>%
  pivot_longer(cols = starts_with("post"), names_to = "posterior") %>%
  group_by(infant_serial_no, posterior) %>% summarise(PEI = 1 - prod(1-value)) 

Figure2b <- 
  PEI %>%
  left_join(df_infants[,c("infant_serial_no","infant_pneumo_status")], by = "infant_serial_no") %>%
  filter(!is.na(infant_pneumo_status)) %>%
  ggplot(aes(x = as.factor(infant_pneumo_status), y = PEI)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width=.07, fill="orange") +
  ylab("PEI value") + xlab("") +
  scale_y_continuous(labels = scales::percent) +
  scale_x_discrete(labels=c("0" = "non-carriers", "1" = "carriers")) +
  theme_classic() 

Figure2 <-
  Figure2a + 
  annotation_custom(ggplotGrob(Figure2b),
                    xmin = 3.1,
                    xmax = 10,
                    ymin = .38,
                    ymax = .85)
ggsave("Figures/Figure2.pdf", unit = "cm", width=20, height=12)                    



###############################################################################
### Second Jags Model (Logistic Regression: 
### Y -> Pneumo Carriage, Covariates -> PEI, Infant Age, Levels -> Commune)
###############################################################################

#Read data in preparation for Rjags -> collect these in a list called 'jdat'

PEI_av <-
  PEI %>% 
  group_by(infant_serial_no) %>% summarise(m = mean(PEI, na.rm = T),
                                           s = sd(PEI, na.rm = T)) %>%
  left_join(df_infants, by = "infant_serial_no") %>%
  filter(m>0)


#The commune IDs should be converted to factors
PEI_av$CommuneFactor <- as.factor(PEI_av$commune_id)


jdat <-
  list(
    "pneu_carr" = PEI_av$infant_pneumo_status,
    "mu" = PEI_av$m,
    "sigma" = PEI_av$s,
    "age"= PEI_av$infant_age_months,
    "commune" = PEI_av$CommuneFactor
  )

#Define model
jcode <- "model{
	#Likelihood
	#As we have assumed a beta distribution for PEI (exp.prop), we need to specify a[i] and b[i] appropriately
	for (i in 1:length(pneu_carr)){
    pneu_carr[i] ~ dbern(prob_carr[i])
    a[i] <- ((1 - mu[i]) / sigma[i] ^ 2 - 1 / mu[i]) * mu[i] ^ 2
    b[i] <- a[i] * (1 / mu[i] - 1)
    #Now, exp.prop will take a beta distribution
    exp.prop[i] ~ dbeta(a[i], b[i])
    #And finally, specify the logistic regression equation for infant carriage
    logit(prob_carr[i]) =  beta0 + beta1 * exp.prop[i] + beta2[commune[i]] + beta3 * age[i]
	}
	#Define our priors
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
  
  #There are 27 communes that we include in the contact study and so...
  
  for(i in 1:27) {
      beta2[i] ~  dnorm(beta2_mu, beta2_sigma)
  }

  #If we include beta0, then beta2_mu by definition should be 0 or we get identifiability issues
  beta2_mu <- 0
  #We assume the standard deviation follows a gamma distribution 
  beta2_sigma ~ dgamma(1,1)

  #Prior for beta3 (coefficient of infant age)
  beta3 ~ dnorm(0,1)
  
}"


#Fit and draw posterior samples
mcmc.length = 10000
jmod = jags.model(
  textConnection(jcode),
  data = jdat,
  n.chains = 4,
  n.adapt = 1000
)
update(jmod)
jpos = coda.samples(jmod,  c("beta0","beta1","beta2","beta3"), n.iter = mcmc.length)
plot(jpos) # check convergence

#Show the coefficients of the logistic regression
print(beta0_est <-
        jpos[[1]][, "beta0"] %>% quantile(probs = c(.5, .025, .975)))
print(beta1_est <-
        jpos[[1]][, "beta1"] %>% quantile(probs = c(.5, .025, .975)))

print( beta1_est <- jpos[[1]][,"beta2[1]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[2]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[3]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[4]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[5]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[6]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[7]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[8]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[9]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[10]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[11]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[12]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[13]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[14]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[15]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[16]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[17]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[18]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[19]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[20]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[21]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[22]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[23]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[24]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[25]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[26]"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta2[27]"] %>% quantile(probs=c(.5,.025,.975)) )

print( beta1_est <- jpos[[1]][,"beta3"] %>% quantile(probs=c(.5,.025,.975)) )


