## save to csv file
library(readr)
load("data/CLEOPATRA_2A.rda")
dat <- CLEOPATRA_2A
dat$arm <- factor(ifelse(dat$arm == "control", "control", "experimental"))
readr::write_csv(dat, "data/CLEOPATRA_2A.csv")
