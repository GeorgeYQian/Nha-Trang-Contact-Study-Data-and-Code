
#####################################

# Dear Reader,

# This file contains the code required to obtain the results shown in the paper
# (the cubic p-spline fit to age and Bayesian logistic regression)

# Please do address any questions to me (george): george.qian08@alumni.imperial.ac.uk
# Thank you for your time and interest in this project. We hope you find it helpful.

# G Qian, 2021 on behalf of the authors of: "Pneumococcal exposure routes for infants, a nested cross-sectional survey in Nha Trang, Vietnam"

######################################


######################################
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

######################################


#Load data from the contact study
df_contacts <- read.csv("df_contacts_public.csv")

#Define age grouping and arrange infant contacts into these groups
#Here, age groups are simply 1 year 'bands' (e.g. 1 yr old, 2 yrs etc)
brks <- 0:100
df_contacts$ctagegp <- cut(df_contacts$contact_age,brks, right=F)


#Load relevant data from the carriage study
df_infants <- read.csv("df_infants_public.csv")



###################################
### First Jags Model (cubic spline)
###################################

# CarriageData here uses the data from Carla's Thesis.
CarriageDataThesis3<-read.csv("CarriageDataCarlaHalfYear.csv")
#The following rates were calculated using Table 6.1 in Carla Talarico's thesis (page 179),
#available online: https://deepblue.lib.umich.edu/handle/2027.42/64706)
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

carriageR<- c(RateCarInf46months,RateCarInf68months,RateCarInf810months,RateCarInf1012months,RateCarTod1214months,RateCarTod1416months,RateCarTod1618months,RateCarTod1820months,RateCarTod2022months,RateCarTod2224months) #read.csv("CarlaDataCarriageRate.csv")
age <-  c(5/12,7/12,9/12,11/12,13/12,15/12,17/12,19/12,21/12,23/12)

NTrang4thUnder2 <- data.frame(age,carriageR)

CarriageDataThesis3<-CarriageDataThesis3 %>% 
  rename(
    CarriageStatus = CarriageStatusH,
    AllAges= AllAgesH
  )

CarriageDataThesis3_ <- CarriageDataThesis3 %>%
  count(AllAges, CarriageStatus) %>%
  pivot_wider(names_from = CarriageStatus,
              values_from = n, values_fill = 0) %>%
  mutate(N = `0` + `1`) %>%
  mutate(y = `1`,
         p = y/N) %>%
  nest(data = c(y,N)) %>%
  mutate(CI = map(data, ~prop.test(.x$y, .x$N)$conf.int)) %>%
  mutate(CI = map(CI, ~data.frame(value = .x, key = c("L", "U")))) %>%
  mutate(CI = map(CI, ~pivot_wider(.x, values_from  = value,
                                   names_from = key))) %>%
  unnest(CI)

CarriageDataThesis3_ %>%
  ggplot(data = ., aes(x = AllAges, y = p)) + 
  geom_pointrange(aes(ymin = L, ymax = U)) +
  geom_line()

bspline <- function(x, K, bdeg=3, cyclic=FALSE, xl=min(x), xr=max(x)){
  x <- as.matrix(x,ncol=1)
  
  ndx <- K - bdeg
  
  # p-spline algorithm, as outlined in Eilers and Marx (1996)
  dx <- (xr - xl) / ndx
  t <- xl + dx * (-bdeg:(ndx+bdeg))
  T <- (0 * x + 1) %*% t
  X <- x %*% (0 * t + 1)
  P <- (X - T) / dx
  B <- (T <= X) & (X < (T + dx))
  r <- c(2:length(t), 1)
  
  for (k in 1:bdeg){
    B <- (P * B + (k + 1 - P) * B[ ,r]) / k; 
  }
  
  B <- B[,1:(ndx+bdeg)]
  
  if (cyclic == 1){
    for (i in 1:bdeg){
      B[ ,i] <- B[ ,i] + B[ ,K-bdeg+i]    
    }
    B <- B[ , 1:(K-bdeg)]
  }
  
  return(B)
}

K = 6

knots <- c(0, 3, 6, 9, 12, 18, 24)

X <- bs(x = unlist(CarriageDataThesis3_$AllAges), 
        knots = knots,
        intercept = FALSE)

x.pred <- seq(0, 99, length.out=100)

X.pred <- bs(x = x.pred, 
             knots = knots, intercept = FALSE)

makeQ = function(degree, K, epsilon=1e-3){
  x <- diag(K)
  E <- diff(x, differences=degree)
  return( t(E) %*% E + x*epsilon)
} 

Q <- makeQ(2, nrow(X))


attributes(X) <- attributes(X)["dim"]
attributes(X.pred) <- attributes(X.pred)["dim"]

#Prepare the data for RJags to use
dat <- CarriageDataThesis3_ %>%
  unnest(data) %>%
  select(y, N) %>% as.list %>%
  append(., values = list(K = ncol(X), Q = Q, 
                          X.pred = X.pred, 
                          X      = X,
                          n = nrow(X),
                          m = nrow(X.pred),
                          beta.0 = matrix(data = rep(0, ncol(X)), 
                                          ncol = 1)))

#Implement RJags to find the cubic spline fit
mod <- jags.model(file = "spline_model.R", data = dat, n.adapt = 5e3)

burn <- jags.samples(mod, variable.names = "p.rep", n.iter = 1e4)

pred <- coda.samples(mod, "p.rep", n.iter = 1e4)

post <- coda.samples(mod, c("beta", "beta.00"), n.iter = 1e4)

#Now we can plot the cubic spline
#Define dot size for plotting
sizeDot <- 3
  
PlotSpline <- tidy(pred) %>%
  bind_cols(data.frame(x = x.pred)) %>%
  ggplot(data= ., aes(x=x, y = mean)) +
  # geom_point(data= CarriageDataThesis3_,
  #            aes(x = AllAges, y= p),color="red") +
  geom_ribbon(aes(ymin = `2.5%`,
                  ymax = `97.5%`),
              alpha = 0.5,
              fill = "lightskyblue") +
  geom_line() + 
  geom_point(data= CarriageDataThesis3_,
             aes(x = AllAges, y= p),color="red", size=sizeDot) + geom_linerange(data=CarriageDataThesis3_, aes(x=AllAges,ymin = L, ymax = U, y=p))+
  xlab("Age") +
  theme_bw() +
  ylab("Prevalence") +  
  scale_x_continuous(trans = "sqrt", #we will use a square root scale on the x-axis to make the lower ages more prominent 
                     breaks = c(1:10)^2) +
  geom_point(data = data.frame(x = knots, mean=0), color="grey")+xlab("Age") +geom_point(data=NTrang4thUnder2,aes(age,carriageR),colour="darkgreen", size=sizeDot)+scale_fill_discrete(breaks=c("trt1","ctrl","trt2"))+ labs(y= "Carriage Prevalence", x = "Age (years)")
#N.B. NTrang4thUnder2 contains data from Carla Talarico's Thesis, which are contained in the file 'CarriageDataCarlaHalfYear.csv')

#Show plot
PlotSpline


statsBeta<-tidy(post)
plot(post)



#####################
### Define Data Block   
#####################

#Now take the output of this spline and use them to calculate PEI.
#These will be inputs to the second JAGS model

#For each set of beta parameters in 'post':
# {Get the formula of the cubic spline
#  Calculate the carriage prevalence at all 61 ages (0-60)
#  Calculate PEI (remember, only the prevalence changes across samples of 'post')
#  Generate df_expFINAL dataframe for each run of for loop
#  Use as input into the PEI Jags/6

# }

number_iter = 10000
max_age = 100 #Maximum age of contacts (set to 100 here)
#pred[[1]] is a matrix [number_iter=10000, max_age=100] containing the runs of the MCMC (10000) and 100 corresponding to age
pred[[1]][number_iter, max_age]

predStore<-tidy(pred) %>%
  bind_cols(data.frame(x = x.pred))

AgeVector<-predStore$x

#mcmc_to_dt is a neat function that we can access here: https://rdrr.io/github/njtierney/mmcc/man/mcmc_to_dt.html 
pred_wider<-pred %>% 
  mcmc_to_dt()   %>% 
  pivot_wider(names_from = parameter, values_from = value)

#Calculate PEI for all N:

#First, set up the vectors required
carriageIter1 = pred_wider[1,-(1:2)] 
carriageIter1ExtVec <- unlist(carriageIter1, use.names=FALSE)
AllAges <- 0:max_age - 1
carriageIterVectors <- vector()

#Now, we will use a loop to obtain the values
for (r in 1:number_iter) {
  carriageVector <- unlist(pred_wider[r,-(1:2)],use.names=FALSE)
  carriageIterVectors <- rbind(carriageIterVectors,carriageVector)
}

#Set up the age bins 
carriageIterVectorsAge <- as.data.frame(t(carriageIterVectors))
carriageIterVectorsAge$agegroups <- levels(cut(df_contacts$contact_age,brks, right=F))
storePEIValues <- ((carriageIterVectorsAge))

#change column names now
#define a list of varying "varname"
varname <- c( 'Iteration','ctagegp')
#Find how many times "varname" repeats itself
n <- c(10000, 1) 

#Replace column name
names(storePEIValues)[1:ncol(storePEIValues)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), varname, n))
names(storePEIValues)[10001] <- "ctagegp"
#Check
names(storePEIValues)

#Keep the contact age at the beginning of the dataframe
storePEIValuesRenamed<-storePEIValues %>%
  select(ctagegp, everything())

#Rename the 1st column (contact's age group)
names(storePEIValuesRenamed)[1] <- "ctagegp"
#Check to see if column names make sense
names(storePEIValuesRenamed)

#Set up another dataframe that gives us the number of contacts of each age for each infant
#and append this to the previous dataframe
store_PEI_Age = df_contacts %>% 
  group_by(contact_serial_no ,ctagegp) %>% summarise(N = length(contact_age)) %>% 
  inner_join(y=storePEIValuesRenamed, by="ctagegp" )

#First, test by calculating PEI for a single iteration
ProbNonInfected     <- 1 - store_PEI_Age[,4]
ContactsPerAgeGroup <- store_PEI_Age$N
ProbEscape_df       <- ProbNonInfected ^ ContactsPerAgeGroup
#This is the probability of 'escape' from contacts in each age group
store_PEI_Age['escape'] <- ProbEscape_df
#Now specify the chance of having at least 1 effective contact (i.e. PEI)
ProbExposure_df <- store_PEI_Age %>%
  group_by(contact_serial_no)   %>%
  summarise(exp.prop = 1 - prod(escape))

#Now calculate all PEI values
ProbExposureStore_df <- as.data.frame(ProbExposure_df$contact_serial_no)
names(ProbExposureStore_df [1])<- "ctslno"
gc()

store_PEI_Age %>% 
  pivot_longer(!c("ctslno","ctagegp", "N"), names_to = "Posterior", values_to = "carriagePrev") %>%
  mutate(escape = (1-carriagePrev)^N) %>%
  group_by(contact_serial_no) %>% summarise(exp.prop = 1 - prod(escape)) -> ProbExposureStore_df

#For loop that calculates PEI by taking number_iter samples of carriage rates from the spline fit above 
for (i in 1:number_iter) {
  df_tmp = store_PEI_Age[,c(1,2,3,i+3)]
  ProbNonInfected     <- 1 - df_tmp[,4]
  ContactsPerAgeGroup <- df_tmp$N
  ProbEscape_df       <- ProbNonInfected ^ ContactsPerAgeGroup
  
  df_tmp['escape'] <- ProbEscape_df[1]
  
  tempo<- df_tmp %>%
    group_by(contact_serial_no)   %>%
    summarise(exp.prop = 1 - prod(escape))
  
  #ProbExposureStore_df$PEI[i]<-tempo
  
  ProbExposureStore_df<-cbind(ProbExposureStore_df,tempo[,2])
}

#We notice that ProbExposureStore_df is a dataframe whose dimensions are 1583 by number_iter+1

#Change column names now (makes things easier to handle later)
#Define a list of varying "varname"
varname <- c( 'contact_serial_no','PEI')
#define how many times above "varname" repeat itself
n <- c(1, 10000) 

#Replace column name
names(ProbExposureStore_df)[1:ncol(ProbExposureStore_df)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep="_"), varname, n))
#Drop the first column (lots of names)
ProbExposureStoreOnlyNum <- ProbExposureStore_df[ -c(1) ]
#Calculate the mean PEI value for each individual (i.e. calculate row mean)
MeanPEI<-rowMeans(ProbExposureStoreOnlyNum,na.rm=TRUE)
#Calculate the standard deviation of the PEI values
stdPEI<-apply(ProbExposureStoreOnlyNum,1,sd)
#Median of individuals' PEI
medianPEI<-apply(ProbExposureStoreOnlyNum,1,median)



#Append the mean and standard deviation of each infant's PEI values to their carriage status
df_infants_PEI<-data.frame(df_infants,MeanPEI,stdPEI)


ProbExposureStoreOnlyNumReshape <- melt(ProbExposureStoreOnlyNum)  #the function melt reshapes it from wide to long
ProbExposureStoreOnlyNumReshape$rowid <- 1:1583
#Inspect the first 10 rows
head(ProbExposureStoreOnlyNumReshape, 10)

#Take a look at one individual's PEI values by plotting them with a histogram
PEI1<-ProbExposureStoreOnlyNumReshape %>% slice(1:10001)
ggplot(PEI1, aes(x=value)) + geom_histogram(binwidth=0.1)
# We can also use a boxplot
pPEI1 <- ggplot(PEI1%>%slice(1:1000), aes(x=variable,y=value)) + 
  geom_violin()
pPEI1  + 
  geom_boxplot(width=0.1, fill="orange") +  
  theme(
    strip.background = element_blank(),
  ) + ylab("PEI value") + xlab("Individual A")+theme(  plot.title = element_text(hjust = 0.5), #axis.title.x=element_blank(Raccoon), 
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank(),
                                                       strip.text = element_text(size=15),
                                                       text = element_text(size = 15)) +
  scale_y_continuous(limits = c(0,1))+ theme(
    # Hide panel borders and remove grid lines
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black")
    
  )


###########################################################################
### Second Jags Model (Logistic Regression: Y -> Pneumo Carriage, X -> PEI)
###########################################################################

#Read data in preparation for Rjags -> collect these in a list called 'jdat'

jdat <- list("pneu_carr" = df_infants_PEI$infant_pneumo_status, "mu" = df_infants_PEI$MeanPEI, "sigma"= df_infants_PEI$stdPEI)

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
    logit(prob_carr[i]) =  beta0 + beta1 * exp.prop[i]
	}
	
	#priors
  beta0 ~ dnorm(0,1)
  beta1 ~ dnorm(0,1)
}"


#Fit and draw posterior samples
mcmc.length=10000
jmod = jags.model(textConnection(jcode), data=jdat, n.chains=4, n.adapt=1000) 
update(jmod)
jpos = coda.samples(jmod, c("beta0","beta1"), n.iter=mcmc.length)
plot(jpos) # check convergence

#Show the coefficients of the logistic regression
print( beta0_est <- jpos[[1]][,"beta0"] %>% quantile(probs=c(.5,.025,.975)) )
print( beta1_est <- jpos[[1]][,"beta1"] %>% quantile(probs=c(.5,.025,.975)) )

