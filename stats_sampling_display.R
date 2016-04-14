###
# Script for displaying stats sampling tests from C++ code
##########################################################

# Recovery_After_Replication.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Recovery_After_Replication.txt", header = FALSE)
hist( df$V1 , col="blue", main="Recovery_After_Replication") 

# Poisson.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Poisson.txt", header = FALSE)
hist( df$V1 , col="blue", main="Poisson") 

# Z.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Z.txt", header = FALSE)
hist( df$V1 , col="blue", main="Z") 


# Binom Dy.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Binomaial_Dying.txt", header = FALSE)
hist( df$V1 , col="blue", main="Binomial_Dying") 


# Binom NB.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Binomaial_Newborn.txt", header = FALSE)
hist( df$V1 , col="blue", main="Binomial_Newborn") 

# Binom Mut.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Binomaial_Mut.txt", header = FALSE)
hist( df$V1 , col="blue", main="Binomial_Mut") 

# Uniform MR.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Uniform_MR.txt", header = FALSE)
hist( df$V1 , col="blue", main="Uniform_MR") 

# Update PR Pqarams.txt
df <- read.table("~/Documents/TEM/src/core/Tests/Update_Proliferation_Rate.txt", header = FALSE)
hist( df$V1 , col="blue", main="Update_Proliferation_Rate") 





