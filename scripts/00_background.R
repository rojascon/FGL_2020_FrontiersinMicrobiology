#################################################################################
#
#               Anaerobic Microbial Communities in a Stratified Sulfidic Lake
#                      
#              Rojas et al 2021.Organic electron donors and terminal electron 
#       acceptors structure anaerobic microbial communities and interactions in 
#                     a permanently stratified sulfidic lake
#
#                               By: Connie Rojas
#                               Created: 4 Aug 2020
#                            Last updated: 8 Feb 2021
################################################################################

##CODE FOR: configuring R workspace and printing R version and package versions
#for reader

################################################################################
#             1.  Configure the workspace for subsequent R project scripts                 
################################################################################

#set conditions for R session
rm(list=ls());
options(scipen=999);
options(stringsAsFactors = FALSE);
options(show.error.locations = TRUE); # show line numbers on error

#load necessary packages
library(pacman);
pacman::p_load("car","MASS","dplyr","tidyr","reshape2","vegan","ggplot2",
               "lme4","lmtest","multcomp","FSA");


################################################################################
#             2. Communicate the R version and package versions to reader                 
################################################################################
print("This code was developed with R version 3.6.2");

print("The packages used and their versions were: multcomp_1.4-15| lmtest_0.9-38| 
lme4_1.1-26| ggplot2_3.3.3| vegan_2.5-7| reshape2_1.4.4| tidyr_1.1.2| dplyr_1.0.3| 
MASS_7.3-53| car_3.0-10| pacman_0.5.1| FSA_0.8.32");

