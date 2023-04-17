rm(list=ls())

setwd("Scripts")

library(devtools)

# Installing updated version of M3C
install("M3C", force = TRUE)
print(packageVersion("M3C"))

# Installing latest version of fake
untar("fake_1.4.0.tar.gz")
install.packages("fake")
print(packageVersion("fake"))

# Installing latest version of sharp
untar("sharp_1.4.0.tar.gz")
install("sharp", force = TRUE)
print(packageVersion("sharp"))
