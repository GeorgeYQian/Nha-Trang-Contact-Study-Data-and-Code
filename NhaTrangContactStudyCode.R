

#####################################
# This file contains the code required to obtain the results shown in the paper
# (the cubic p-spline fit to age and Bayesian logistic regression)
#
# Please address any questions to me (george): george.qian08@alumni.imperial.ac.uk
#
# G Qian, 2021 on behalf of the authors of: 
# "Pneumococcal exposure routes for infants, a nested cross-sectional survey in Nha Trang, Vietnam"
######################################
set.seed(31415)

###################################
### Load libraries
###################################
#Import libraries
library(tidyverse)
library(rjags)
library(binom)
library(runjags)
library(splines)
library(reshape2)
library(gghalves)
library(knitr)
options(dplyr.summarise.inform = FALSE)

if (!require(mmcc)){
    prompt <- menu(c("Yes", "No"), title="Do you wish to install the mmcc package from github?")
    
    if (prompt == 1L){
        remotes::install_github("njtierney/mmcc")
        library(mmcc)
    }
    
}


###################################
### Load data
###################################
# Load data from the contact study
df_contacts <- 
    read.csv("contacts_public.csv") %>%
    mutate(ctagegp = cut(contact_age, 0:100, right = F))

# Load relevant data from the carriage study
df_infants <- read.csv("infants_public.csv")

# CarriageData here uses the data from Carla's Thesis.
#The following rates were calculated using Table 6.1 in Carla Talarico's thesis (page 179),
#available online: https://deepblue.lib.umich.edu/handle/2027.42/64706)
CarriageDataThesis3 <- read.csv("Talarico.csv")

# Carriage prevalence observed in the trial
NTrang4thUnder2 <- 
    tibble (age = c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23)/12,
            carriageR = c(0.147, 0.183, 0.235, 0.278, 0.276, 0.352, 0.344, 0.355, 0.413, 0.420))

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


#######################################
### Figure 1 (contact age distribution)
#######################################

Figure1 <- df_contacts %>%
    mutate(contact_age = floor(contact_age)) %>%
    group_by(contact_age,contact_time) %>% 
    summarise(N = length(contact_age)) %>%
    ggplot(aes(x = contact_age, 
               y = N, 
               fill = factor(contact_time),
               group = factor(contact_time))) +
    geom_bar(stat="identity")+
    xlab("Contact Age (years)") + 
    ylab("Number of Contacts") + 
    scale_fill_manual(values  = rev(tail(RColorBrewer::brewer.pal(4,"Reds"),-1)), 
                      name    = "Contact Time", 
                      labels  = c("Long (>60 minutes)",
                                  "Medium (5-60 minutes)",
                                  "Short (<5 minutes)")) +  
    theme_classic() +
    theme(legend.position = c(0.8, 0.8)) 

ggsave("Figures/Figure1.pdf", plot = Figure1, unit = "cm", width=15, height=9)


###################################
### First Jags Model (cubic spline)
###################################

K = 6
knots <- c(0, 3, 6, 9, 12, 18, 24)

X <- bs(
    x = unlist(CarriageDataThesis3_$AllAges),
    knots = knots,
    intercept = FALSE
)

x.pred <- seq(0, 99, length.out = 100)

X.pred <- bs(x         = x.pred,
             knots     = knots,
             intercept = FALSE)

makeQ = function(degree, K, epsilon = 1e-3) {
    x <- diag(K)
    E <- diff(x, differences = degree)
    return(t(E) %*% E + x * epsilon)
}

Q <- makeQ(2, nrow(X))

attributes(X)      <- attributes(X)["dim"]
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

# Model diagnostics
statsBeta <- tidy(post)

# Now we can plot the spline
Figure2a <- tidy(pred) %>%
    bind_cols(data.frame(x = x.pred)) %>%
    ggplot(data = ., aes(x = x, y = mean)) +
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
    xlab("Age of infant (years)") +
    theme_classic() +
    ylab("Prevalence") +
    #we will use a square root scale on the x-axis to make the lower ages more prominent
    scale_x_continuous(trans = "sqrt",
                       breaks = c(0:10) ^ 2) +
    geom_point(data = data.frame(x = knots, mean = 0),
               color = "grey") +
    xlab("Age") +
    geom_point(data = NTrang4thUnder2,
               aes(age, carriageR),
               colour = "darkgreen",
               size = 3) + 
    scale_fill_discrete(breaks = c("trt1", "ctrl", "trt2")) + 
    labs(y = "Carriage Prevalence", x = "Age (years)") 


#####################
### Define Data Block
#####################

#Now take the output of this spline and use them to calculate PEI.
#These will be inputs to the second JAGS model

#For each set of beta parameters in 'post':
#  Calculate the carriage prevalence at all 101 ages (i.e. 0-100) from the cubic spline
#  Calculate PEI 
#  Use as input into the PEI Jags


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
    left_join(pivot_longer(pred_prevalence,
                           cols = starts_with("post"),
                           names_to = "posterior"),
              by = "ctagegp") %>%
    mutate(posterior = parse_number(posterior)) %>%
    group_by(infant_serial_no, posterior) %>% summarise(PEI = 1 - prod(1-value)) 

PEI_summary_by_carrier_status <- PEI %>%
    left_join(df_infants[,c("infant_serial_no","infant_pneumo_status")],
              by = "infant_serial_no") %>%
    filter(!is.na(infant_pneumo_status)) %>%
    mutate(infant_pneumo_status = ifelse(infant_pneumo_status == 1L,
                                         "Carriers",
                                         "Non-carriers")) %>%
    group_by_at(.vars = vars(-c(posterior, PEI))) %>%
    nest %>%
    mutate(Q = map(data, ~quantile(.x$PEI, c(0.025,
                                             0.5,
                                             0.975), na.rm = T))) %>%
    unnest_wider(Q) %>%
    ungroup

Figure2b <- PEI_av %>%
    filter(!is.na(infant_pneumo_status)) %>%
    mutate(
        infant_age_months = factor(infant_age_months),
        infant_age_months = fct_recode(infant_age_months, "12+" = "12",
                                       "12+" = "13"),
        infant_pneumo_status = ifelse(infant_pneumo_status == 1L,
                                      "Carriers",
                                      "Non-carriers")) %>%
    
    {ggplot(data = .,
            aes(x = factor(infant_age_months),
                y = m)) +
            geom_half_violin(data = filter(., infant_pneumo_status == "Non-carriers"),
                             aes(fill = "Non-carriers"),
                             color  = NA, side = 'l') +
            geom_half_violin(data = filter(., infant_pneumo_status == "Carriers"),
                             aes(fill = "Carriers"),
                             color  = NA, side = 'r')} +
    theme_classic() +
    scale_y_continuous(limits = c(0,1),
                       labels = scales::percent,
                       name = "Mean PEI") +
    xlab("Age of infant (years)") +
    scale_fill_manual(name = NULL,
                      values = c("Non-carriers" = "#BCBDDC",
                                 "Carriers"     = "#756BB1")) +
    theme(legend.position = 'top')


Figure2 <-
    Figure2a + 
    annotation_custom(ggplotGrob(Figure2b),
                      xmin = 3.1,
                      xmax = 10,
                      ymin = .35,
                      ymax = .95)

ggsave("Figures/Figure2.pdf", plot = Figure2, unit = "cm", width=15, height=9)                   


###########################################################################
### Second Jags Model (Logistic Regression: Y -> Pneumo Carriage, X -> PEI)
###########################################################################

# Read data in preparation for Rjags -> collect these in a list called 'jdat'

PEI_av <-
    PEI %>% 
    group_by(infant_serial_no) %>% summarise(m = mean(PEI, na.rm = T),
                                             s = sd(PEI, na.rm = T)) %>%
    left_join(df_infants, by = "infant_serial_no") %>%
    filter(m>0) %>%
    mutate(commune_id = as.factor(commune_id))

jdat <-
    list(
        "pneu_carr" = PEI_av$infant_pneumo_status,
        "mu" = PEI_av$m,
        "sigma" = PEI_av$s,
        "age"= PEI_av$infant_age_months,
        "commune" = PEI_av$commune_id
    )

#Define model
jcode <- "model{

	#Define Likelihood
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
  beta0 ~ dnorm(0,1e-3)
  beta1 ~ dnorm(0,1e-3)
  
  #There are 27 communes that we include in the contact study and so...
  
  for(i in 1:27) {
      beta2[i] ~  dnorm(beta2_mu, beta2_sigma)
  }

  #If we include beta0, then beta2_mu by definition should be 0 or we get identifiability issues
  beta2_mu <- 0
  
  #We assume the standard deviation follows a gamma distribution 
  beta2_sigma ~ dgamma(2,2)

  #Prior for beta3 (coefficient of infant age)
  beta3 ~ dnorm(0,1e-3)
  
}"


#Fit and draw posterior samples
mcmc.length = 10000
jmod = jags.model(
    textConnection(jcode),
    data = jdat,
    n.chains = 4,
    n.adapt = 1000
)
# update(jmod) # what's this line for?
jpos = coda.samples(jmod,  c("beta0","beta1","beta2","beta3"), n.iter = mcmc.length)

# show regression coefficients
tidy(jpos)


###########################################################################
### Figure 3
###########################################################################


# Calculate exposure across all age groups 
# Exposure from any one age group is defined as:
# sum(number of contacts in age group*P(carriage|age))

exposure <-
    df_contacts %>%
    tibble() %>%
    select(-c(X, contact_serial_no)) %>%
    left_join(pred_prevalence, by = "ctagegp") %>%
    pivot_longer(cols = starts_with("post"), names_to = "posterior") %>%
    mutate(ctagegp = cut(contact_age, breaks = c(seq(0,65, by=5),100), right=F)) %>%
    group_by(ctagegp, posterior) %>% summarise(exp = sum(value)) %>%
    filter(!is.na(ctagegp)) %>%  #Now obtain all posterior estimates to find confidence intervals
    group_by(posterior) %>% mutate(exp_prop = exp / sum(exp)) %>%
    group_by(ctagegp) %>% summarise(med = median(exp_prop),
                                    lo  = quantile(exp_prop, 0.025),
                                    hi  = quantile(exp_prop, 0.975))

s = seq(0,60, by=5)

ages = c(paste0(s,"-",s+4),"Over 65")

# Hone in on the proportion of exposure due to 0-15 year olds alone
# (with each year group in a separate bin)

exposure_fine <-
    df_contacts %>%
    tibble() %>%
    select(-c(X, contact_serial_no)) %>%
    left_join(pred_prevalence, by = "ctagegp") %>%
    pivot_longer(cols = starts_with("post"), names_to = "posterior") %>%
    group_by(contact_age, posterior) %>% summarise(exp = sum(value)) %>%
    filter(!is.na(exp)) %>%
    group_by(posterior) %>% mutate(exp_prop = exp / sum(exp)) %>%
    group_by(contact_age) %>% summarise(med = median(exp_prop),
                                        lo = quantile(exp_prop, 0.025),
                                        hi = quantile(exp_prop, 0.975))


# Plot Figure 3a (proportion of exposure due to every age group)

Figure3a <- 
    exposure %>%
    ggplot(aes(x = ctagegp, y = med, ymin = lo, ymax = hi)) +
    geom_bar(stat = "identity", fill = "#9E9AC8", color = 'black') +
    geom_linerange(color = "black", lwd=1.2) + 
    ylab("Proportion of exposure") + xlab("Age (years)") + 
    scale_x_discrete(labels = ages) +
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust = 1, angle = 45))


# Plot Figure 3b (showing the proportion of exposure due to 0-15 year olds)

Figure3b <- 
    exposure_fine %>%
    filter(contact_age < 15) %>%
    ggplot(aes(x = contact_age, y = med, ymin = lo, ymax = hi)) +
    geom_bar(stat = "identity", fill = "#9E9AC8", color = 'black') +
    geom_linerange(color = "black", lwd=1.2) + 
    scale_x_continuous(breaks = 0:14) +
    theme_classic(base_size = 8) +
    theme(axis.title = element_blank())

Figure3 <-
    Figure3a + 
    annotation_custom(ggplotGrob(Figure3b),
                      xmin = 4,
                      xmax = 14,
                      ymin = .2,
                      ymax = .55)

ggsave("Figures/Figure3.pdf", plot = Figure3, unit = "cm", width=15, height=9)                    

## Mann-Whitney U tests of PEI by infant age

PEI_summary_by_carrier_status %>% 
    left_join(select(PEI_av, infant_serial_no, infant_age_months)) %>%
    select(infant_pneumo_status, `50%`, Age = infant_age_months) %>%
    group_by(infant_pneumo_status, Age) %>%
    nest %>%
    spread(infant_pneumo_status, data) %>%
    mutate(mwu = map2(.x = Carriers,
                      .y = `Non-carriers`,
                      .f = ~wilcox.test(x = .x$`50%`, y = .y$`50%`))) %>%
    transmute(p = map_dbl(.x = mwu, 'p.value')) %>%
    mutate(Sig. = cut(p,
                      c(0,0.001, 0.01, 0.05, 0.1, 1),
                      labels = c("***", "**", "*", ".", ""))) %>%
    knitr::kable(digits = 3, format = 'simple')

