

mydata <- mydata[,-c(1)]
cor.mydata <- cor(mydata)
dim(cor.mydata)

n <- dim(cor.mydata)[1]

cor.lin <- vector(mode = 'numeric', length = (n*(n-1)/2))
length(cor.lin)

count <- 1
for(i in 1:(n-1)){
  for(j in (i+1):n){
    cor.lin[count] <- cor.mydata[i, j]
    count <- count+1
  }
}
max(cor.lin)
hist(cor.lin)

