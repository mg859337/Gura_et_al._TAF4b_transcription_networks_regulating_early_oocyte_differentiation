# Getting numbers for downsampling RNA-seq samples

# insert mapped read numbers of E14.5 Wts and Defs
x <- c(103311892, 86368224, 86353695, 101325148, 97394611, 84043745, 90035186, 82526034, 90347649, 88333599, 81302765, 127676050)

# is it a normal distribution? A uniform distribution?
hist(x)

#looks more uniform than anything else
y <- min(x)
z <- max(x)


# Get four random numbers within the parameters

k <- runif(4, min = y, max = z)
hist(k)

setwd("~/E14.5_F_M/")

write.csv(k, "read_depths_randomlygenerated.csv")
