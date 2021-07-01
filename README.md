This repository provides the source code and related data files used in:

    "Pneumococcal exposure routes for infants, a nested cross-sectional survey in Nha Trang, Vietnam."
    [Full journal info and link to be provided]


### Data files

1. df_contacts_public.csv

    This is a csv file that contains relevant information about the reported contacts of the infants enrolled in the contact study. The columns contain:

    - a serial number for each contact (contact_serial_no), a serial number for each infant (infant_serial_no), 
    - the age, in years, of each contact (contact_age) 
    - the contact time (contact_time). Specifically, the contact_time column states whether the contact was reported to be of short (under 5 minutes, contact_time=1), medium (between 5 minutes to 1 hour, contact_time=2) or long (over 1 hour, contact_time=3) duration 

2. df_infants_public.csv

    This csv file contains relevant information about the infants themselves. The columns contain: 

    - a serial number for each infant (infant_serial_no), which matches that found in df_contacts_public.csv
    - the gender of the infant: either male (1) or female (2). This is found under the column heading: gender
    - the number of siblings in the infant's household (under the column heading: siblings)
    - The number of people residing in the infant's household (under the column heading: household_size)
    - The employment status of the infant's primary caregiver (under the column heading: caregiver_employment): 1 if the caregiver does not work, 2 if he/she works less than 2         days a week, and 3 if he/she works more than 2 days/week    
    - The most advanced education level achieved by any member of the infant's household (under the column heading: household_education): 1 if primary school has not been              completed, 2 if primary school has been completed, 3 if junior high school has been completed and 4 if any education after junior high has been attained
    - The most advanced level of mobility that the infant has achieved (under the column heading: infant_mobility): 1 if he/she cannot yet sit, 2 if sitting has been achieved, 3       if he/she can crawl and 4 if he/she can walk
    - the age of the infant in months (infant_age)
    - a number – either 0 or 1 – in infant_pneumo_status, signifying whether the infant is a pneumococcal carrier (1) or is not (0).
    - an ID number representing the commune to which the infant belongs (commune_id)
    

3.	CarriageDataCarlaHalfYear.csv

    This is a csv file that converts Table 6.1 from Carla Talarico’s thesis (page 179, available online: https://deepblue.lib.umich.edu/handle/2027.42/64706) into line data. This is required for reproducing Figure 2 in the manuscript.


### Code

1. NhaTrangContactStudyCode.r

    This is the main R file. Running this produces the p-spline fit to age, calculates PEI values for each infant, and then, using JAGS, performs the Bayesian logistic regression described in the paper. It also reproduces Figures 2 and 3 in the manuscript.

2. spline_model.r

    This is an R file that defines the p-spline (described by Eilers and Marx, 1996) required for the spline fit of pneumococcal carriage data to age.



### RStudio and package versions

The above code was written with R v4.0.3 in mind.

Furthermore, the following package versions were used:

-	binom v1.1-1
-	ggplot2 v3.3.2
-	reshape2 v1.4.4
-	rjags v4-10
-	runjags v2.0.4-6
-	splines v4.0.3
-	tidyverse v1.3.0
