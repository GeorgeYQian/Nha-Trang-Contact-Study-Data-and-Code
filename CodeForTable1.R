#Packages required
install.packages('utf8')
#require(tidyverse)
library(tidyverse)
library(rjags)
library(binom)
library(runjags)
library(mmcc)
library(splines)
library(MASS)
library(lme4)
options(dplyr.summarise.inform = FALSE)



# load data
df_contacts <- read.csv("df_contacts_public.csv")
df_infants <- read.csv("df_infants_public.csv")


#Calculate number of contacts each infant has had
df_contactsNumber<-as.data.frame(table(df_contacts$infant_serial_no))
df_infants$number_contacts <- df_contactsNumber$Freq
df_infants$infant_pneumo_status <- df_infants$infant_pneumo_status

#Gender

#Males
infant_male <- subset(df_infants, gender == 1)
dim(infant_male)[1]
mean(infant_male$number_contacts)
sd(infant_male$number_contacts)

#Females
infant_female <- subset(df_infants, gender == 2)
dim(infant_female)[1]
mean(infant_female$number_contacts)
sd(infant_female$number_contacts)

#Negative Binomial Regression to find the Odds
negBinGender<- (glm.nb(number_contacts ~ gender, data=df_infants, link=log))
#Display Ratio
exp(negBinGender$coefficients[2])


#Age-adjusted analysis
df_infants <- df_infants %>% 
  mutate(agegp = if_else( infant_age_months<8 , 0, 1)) 


AdjustedRatioGender = lmer(number_contacts ~ gender +  agegp + (1 | commune_id), data = df_infants)
 summary(AdjustedRatioGender)


#Family information
#No siblings in family
infant_no_siblings<- subset(df_infants, siblings == 0)
dim(infant_no_siblings)[1]
mean(infant_no_siblings$number_contacts)
sd(infant_no_siblings$number_contacts)


#1 or more siblings
infant_siblings<- subset(df_infants, siblings != 0)
dim(infant_siblings)[1]
mean(infant_siblings$number_contacts)
sd(infant_siblings$number_contacts)


df_infants <- df_infants %>% 
  mutate(sibling_status = if_else(siblings ==0 , 0, 1)) 

#Negative Binomial Regression to find the Odds
negBinSiblings<- (glm.nb(number_contacts ~ sibling_status, data=df_infants, link=log))
#Display Ratio
exp(negBinSiblings$coefficients[2])



#Number of people in the household

#2 to 4 people
infant_household2To4 <- subset(df_infants, household_size>=2 & household_size<=4)
dim(infant_household2To4)[1]
mean(infant_household2To4$number_contacts)
sd(infant_household2To4$number_contacts)

#over 4 people
infant_householdOver4 <- subset(df_infants, so103>4)
dim(infant_householdOver4)[1]
mean(infant_householdOver4$number_contacts)
sd(infant_householdOver4$number_contacts)


df_infants <- df_infants %>% 
  mutate(HH_status = if_else(household_size != 4, 0, 1)) #if anybody in the HH does not have a degree, then 1, otherwise 2


#Negative Binomial Regression to find the Odds
negBinHouseholdSize<- (glm.nb(number_contacts ~ HH_status, data=df_infants, link=log))
#Display Ratio
exp(negBinHouseholdSize$coefficients[2])





#Employment status of primary caretaker

#Not in employment
infant_no_employment <- subset(df_infants,caregiver_employment==1)
dim(infant_no_employment)[1]
mean(infant_no_employment$number_contacts)
sd(infant_no_employment$number_contacts)


#In employment
infant_employment <- subset(df_infants,caregiver_employment!=1)
dim(infant_employment)[1]
mean(infant_employment$number_contacts)
sd(infant_employment$number_contacts)

#Negative Binomial Regression to find the Odds
negBinEmployment<- (glm.nb(number_contacts ~ caregiver_employment, data=df_infants, link=log))
#Display Ratio
exp(negBinEmployment$coefficients[2])


# Highest level of education
#Has not reached uni level
infant_education <- subset(df_infants,household_education!=4)
dim(infant_education)[1]
mean(infant_education$number_contacts)
sd(infant_education$number_contacts)

#Uni level
infant_degree <- subset(df_infants,household_education==4)
dim(infant_degree)[1]
mean(infant_degree$number_contacts)
sd(infant_degree$number_contacts)


df_infants <- df_infants %>% 
  mutate(degree_status = if_else(household_education != 4, 1, 2)) #if anybody in the HH does not have a degree, then 1, otherwise 2

#Negative Binomial Regression to find the Odds
negBinEducation<- (glm.nb(number_contacts ~ degree_status, data=df_infants, link=log))
#Display Ratio
exp(negBinEducation$coefficients[2])



#Infant's mobility status

#Can Sit (At least)
infant_sit <- subset(df_infants,infant_mobility!=1)
dim(infant_sit)[1]
mean(infant_sit$number_contacts)
sd(infant_sit$number_contacts)

#Can't yet sit
infant_no_sit<- subset(df_infants,infant_mobility==1)
dim(infant_no_sit)[1]
mean(infant_no_sit$number_contacts)
sd(infant_no_sit$number_contacts)


#
df_infants <- df_infants %>% 
  mutate(sit_status = if_else(infant_mobility != 1, 1, 2)) 
#Negative Binomial Regression to find the Odds
negBinSit<- (glm.nb(number_contacts ~ sit_status, data=df_infants, link=log))
#Display Ratio
exp(negBinSit$coefficients[2])



#Can Crawl (at least)

infant_crawl <- subset(df_infants,infant_mobility!=1&infant_mobility!=2)
dim(infant_crawl)[1]
mean(infant_crawl$number_contacts)
sd(infant_crawl$number_contacts)

#Can't yet crawl
infant_no_crawl<- subset(df_infants,infant_mobility==1|infant_mobility==2)
dim(infant_no_crawl)[1]
mean(infant_no_crawl$number_contacts)
sd(infant_no_crawl$number_contacts)


df_infants <- df_infants %>% 
  mutate(crawl_status = if_else(infant_mobility == 1|infant_mobility==2, 1, 2)) 
#Negative Binomial Regression to find the Odds
negBinCrawl<- (glm.nb(number_contacts ~ crawl_status, data=df_infants, link=log))
#Display Ratio
exp(negBinCrawl$coefficients[2])


#Can Walk

infant_walk <- subset(df_infants,infant_mobility==4)
dim(infant_walk)[1]
mean(infant_walk$number_contacts)
sd(infant_walk$number_contacts)

#Can't yet walk
infant_no_walk<- subset(df_infants,infant_mobility!=4)
dim(infant_no_walk)[1]
mean(infant_no_walk$number_contacts)
sd(infant_no_walk$number_contacts)


df_infants <- df_infants %>% 
  mutate(walk_status = if_else(infant_mobility != 4, 1, 2)) 
#Negative Binomial Regression to find the Odds
negBinWalk<- (glm.nb(number_contacts ~ walk_status, data=df_infants, link=log))
#Display Ratio
exp(negBinWalk$coefficients[2])



#HH Mobility


#Bicycle Usage

infant_bike <- subset(df_infants,bicycle==1)
dim(infant_bike)[1]
mean(infant_bike$number_contacts)
sd(infant_bike$number_contacts)

#No bike usage
infant_no_bike<- subset(df_infants,bicycle!=1)
dim(infant_no_bike)[1]
mean(infant_no_bike$number_contacts)
sd(infant_no_bike$number_contacts)

#Negative Binomial Regression to find the Odds
negBinBike<- (glm.nb(number_contacts ~ bicycle, data=df_infants, link=log))
#Display Ratio
exp(negBinBike$coefficients[2])




#Motorcycle Usage

infant_Mbike <- subset(df_infants,motorbike==1)
dim(infant_Mbike)[1]
mean(infant_Mbike$number_contacts)
sd(infant_Mbike$number_contacts)

#No bike usage
infant_no_Mbike<- subset(df_infants,motorbike=1)
dim(infant_no_Mbike)[1]
mean(infant_no_Mbike$number_contacts)
sd(infant_no_Mbike$number_contacts)

#Negative Binomial Regression to find the Odds
negBinMBike<- (glm.nb(number_contacts ~ motorbike, data=df_infants, link=log))
#Display Ratio
exp(negBinMBike$coefficients[2])


#Car Usage

infant_Car <- subset(df_infants,car==1)
dim(infant_Car)[1]
mean(infant_Car$number_contacts)
sd(infant_Car$number_contacts)

#No car usage
infant_no_Car<- subset(df_infants,car!=1)
dim(infant_no_Car)[1]
mean(infant_no_Car$number_contacts)
sd(infant_no_Car$number_contacts)

#Negative Binomial Regression to find the Odds
negBinCar<- (glm.nb(number_contacts ~ car, data=df_infants, link=log))
#Display Ratio
exp(negBinCar$coefficients[2])


#Does HH Walk?

infant_Walking <- subset(df_infants,walk==1)
dim(infant_Walking)[1]
mean(infant_Walking$number_contacts)
sd(infant_Walking$number_contacts)

#No Walking
infant_no_Walking<- subset(df_infants,walk!=1)
dim(infant_no_Walking)[1]
mean(infant_no_Walking$number_contacts)
sd(infant_no_Walking$number_contacts)

#Negative Binomial Regression to find the Odds
negBinWalking<- (glm.nb(number_contacts ~ walk, data=df_infants, link=log))
#Display Ratio
exp(negBinWalking$coefficients[2])


#Public Transport Usage

infant_PT <- subset(df_infants,pub_transport==1)
dim(infant_PT)[1]
mean(infant_PT$number_contacts)
sd(infant_PT$number_contacts)

#No bike usage
infant_no_PT<- subset(df_infants,pub_transport!=1)
dim(infant_no_PT)[1]
mean(infant_no_PT$number_contacts)
sd(infant_no_PT$number_contacts)

#Negative Binomial Regression to find the Odds
negBinPT<- (glm.nb(number_contacts ~ pub_transport, data=df_infants, link=log))
#Display Ratio
exp(negBinPT$coefficients[2])



#Day Care Attendance

#Attends daycare
infant_Daycare <- subset(df_infants,nursery!=1)
dim(infant_Daycare)[1]
mean(infant_Daycare$number_contacts)
sd(infant_Daycare$number_contacts)

#No daycare attendance
infant_no_Daycare<- subset(df_infants,nursery==1)
dim(infant_no_Daycare)[1]
mean(infant_no_Daycare$number_contacts)
sd(infant_no_Daycare$number_contacts)

df_infants <- df_infants %>% 
  mutate(daycare = if_else(nursery ==1 , 0, 1)) 

#Negative Binomial Regression to find the Odds
negBinDaycare<- (glm.nb(number_contacts ~ daycare, data=df_infants, link=log))
#Display Ratio
exp(negBinDaycare$coefficients[2])


#Pneumo Status

#Has pneumo
infant_pneumo <- subset(df_infants,infant_pneumo_status==1)
dim(infant_pneumo)[1]
mean(infant_pneumo$number_contacts)
sd(infant_pneumo$number_contacts)

#Does not have pneumo
infant_no_pneumo<- subset(df_infants,infant_pneumo_status==0)
dim(infant_no_pneumo)[1]
mean(infant_no_pneumo$number_contacts)
sd(infant_no_pneumo$number_contacts)



#Negative Binomial Regression to find the Odds
negBinPneumo<- (glm.nb(number_contacts ~ infant_pneumo_status, data=df_infants, link=log))
#Display Ratio
exp(negBinPneumo$coefficients[2])




