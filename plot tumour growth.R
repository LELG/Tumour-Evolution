## Read files to plot tumour growth

path <-"~/Documents/TEM/src/core/Tester/TumourFunctions/Mitosis_Model/"

files <- list.files(path = "~/Documents/TEM/src/core/Tester/TumourFunctions/Mitosis_Model", pattern = NULL, all.files = FALSE,
           full.names = FALSE, recursive = FALSE,
           ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE)

df <- read.table("~/Documents/TEM/src/core/Tester/TumourFunctions/Mitosis_Model/te_file_0.txt", header = FALSE)
h <- numeric(0)
s <- numeric(0)

j=1
for(i in files  )
{
  
  val <- paste(path,i,sep="") 
  print(val)
  df <- read.table(val, header = FALSE)
  h[j] <- max(df$V2)
  s[j] <- max(df$V1)
  j=j+1
}

mh <-max(h)/8760
ms <- max(s)
plot(df$V2/8760, df$V1,xlim = c(0, mh),ylim = c(100,ms), col="blue", type  ="l", xlab="years", ylab="Tumour Size", main="Mitosis")

cl <- rainbow(length(files))
j=1
for(i in files  )
{
  
  val <- paste(path,i,sep="")
  df <- read.table(val, header = FALSE)
  lines(df$V2/8760, df$V1,xlim =c(0, mh),ylim =c(0,ms), col=cl[j], type  ="l")
  j=j+1
}

par(mfrow=c(2,1))
df <- read.table("~/Documents/TEM/src/core/Tester/TumourFunctions/Mitosis_Model/te_file_31.txt", header = FALSE)
plot(df$V2/8760, df$V1, col="blue", type  ="l", xlab="years", ylab="Tumour Size", main="Mitosis PR=0.03, DR= 0.03 MR=1 × 10^-6, Laplace")
plot(df$V2/8760, df$V3, col="blue", xlab="years", ylab="Tumour Size", main="Mitosis Clones PR=0.03, DR= 0.03 MR=1 × 10^-, Laplace")


