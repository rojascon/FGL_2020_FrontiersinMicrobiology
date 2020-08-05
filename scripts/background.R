rm(list=ls());                          # clear Console Window
options(show.error.locations = TRUE); # show line numbers on error
options(scipen=999);                   #turns off scientific notation
library(package=pacman)              # load necessary packages
pacman::p_load("FSA","ggplot2","tidyr","vegan")