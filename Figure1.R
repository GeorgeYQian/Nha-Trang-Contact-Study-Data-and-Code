
#################################

#This file reproduces Figure 1 

#################################

#Obtain required packages
require(tidyverse)

#Load data
df_infants <- read.csv("df_infants_public.csv")
df_contacts <- read.csv("df_contacts_public.csv")

#Define age grouping
brks<-0:99
df_contacts$ctagegp <- cut(df_contacts$contact_age,brks, right=F)

#Sort contacts by age and contact time
df_contacts_grouped=df_contacts%>%group_by(ctagegp,contact_time) %>% summarise(N = length(contact_age))

#Make contact time continuous
contact_time_continuous<-as.numeric(df_contacts_grouped$contact_time)
ctagegp_continuous<-as.numeric(df_contacts_grouped$ctagegp)

#Plot using a histogram
ggplot(data=df_contacts_grouped, aes(x=ctagegp_continuous,y=N,fill=factor(contact_time) )) +
  geom_bar(stat="identity")+
  xlab("Contact Age (years)") + ylab("Number of Contacts") + 
  scale_fill_brewer(palette = "Accent", name = "Contact Time", labels = c("Short (<5 minutes)", "Medium (5-60 minutes)", "Long (>60 minutes)")) +  theme(legend.position = "none")+ 
  theme_classic() +
  theme(legend.position = c(0.8, 0.8)) 
ggsave("Figures/Figure1.pdf", unit = "cm", width=15, height=9)
