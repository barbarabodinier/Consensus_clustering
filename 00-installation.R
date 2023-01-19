setwd("Packages")

library(devtools)

# Installing updated version of M3C
install("M3C", force = TRUE)
print(packageVersion("M3C"))

# Installing latest version of fake
install.packages("fake")

# Installing latest version of sharp
untar("sharp_1.3.0.9000.tar.gz")
install("sharp", force = TRUE)
print(packageVersion("sharp"))
