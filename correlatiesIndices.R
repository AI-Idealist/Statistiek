
#--------------------------------------Bereken Z-scores
y=matrix(data=c(115,85,100,460,340,400), nrow=3, ncol=2)
y
sum(y)
rn <-c("IQ","Inkomen")
colnames(y) <- rn
rownames(y) <- c("A","B","C")
y<- as.table(y)
y
scale(y, center = TRUE, scale = TRUE)

#--------------------------------------Bereken CramersV
y=matrix(data=c(4,6,14,14,35,12,22,21,32), nrow=3, ncol=3)
colnames(y) <- c("CDA", "PVDA", "VVD")
rownames(y) <- c("NRC","VK","TG")
y
cramersV(y)

#-------------------------------------Bereken phi
#in libary psych
y=matrix(data=c(4,6,14,35), nrow=2, ncol=2)
colnames(y) <- c("CDA", "PVDA")
rownames(y) <- c("NRC","VK")
y
phi(y)

#--------------------------------------Bereken GoodmanenKruskal Tau
y=matrix(data=c(32,77,56,43,62,61), nrow=3, ncol=2)
colnames(y) <- c("Vrouw", "Man")
rownames(y) <- c("NRC","VK","TG")
y

GKtau <- function(y){

#bereken E2 fouten met info over onafhankelijke)
y_prc = (y / colSums(y))
err = y * (1-y_prc)
E2 = sum (rowSums(err))
E2

#bereken E1 fouten zonder info over de onafhankelijke
rtot = rowSums(y)
rtot
r_perc = (rtot/sum(rtot))
r_perc
err1 = rtot*(1-r_perc)
err1
E1 = sum(err1)
E1

tau = (E1 - E2) / E1
return(tau)
}

GoodManKruskaltau = GKtau(y)
GoodManKruskaltau

#-----------------------------------------------bereken lambda
y=matrix(data=c(5,2,3,4,6,2,1,3,4), nrow=3, ncol=3)
colnames(y) <- c("ONM", "GTST","ATWT")
rownames(y) <- c("breezer","cole","bier")
y


calc.lambda <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  SumRmax <- sum(apply(x, 1, max))
  SumCmax <- sum(apply(x, 2, max))
  MaxCSum <- max(colSums(x))
  MaxRSum <- max(rowSums(x))
  n <- sum(x)

  L.CR <- (SumRmax - MaxCSum) / (n - MaxCSum)
  L.RC <- (SumCmax - max(rowSums(x))) / (n - MaxRSum)
  L.S <- (SumRmax + SumCmax - MaxCSum - MaxRSum) /
        ((2 * n) - MaxCSum - MaxRSum)

  Llist <- list(L.CR, L.RC, L.S)
  names(Llist) <- c("L.CR", "L.RC", "L.S")

  Llist
}

lambda = calc.lambda(y)
lambda


# ---------------------------------------------------Calculate Goodman-Kruskal gamma
# x = table
y=matrix(data=c(177,787,7582,69,1150,738,529,290,465), nrow=3, ncol=3)
colnames(y) <- c("nooit", "soms","vaak")
rownames(y) <- c("vaak","soms","nooi")
y

concordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND > c)
  # for each matrix[r, c]
  mat.lr <- function(r, c)
  { 
    lr <- x[(r.x > r) & (c.x > c)]
    sum(lr)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.lr, r = r.x, c = c.x))
}

# Calculate DIScordant Pairs in a table
# cycle through x[r, c] and multiply by
# sum(x elements below and to the left of x[r, c])
# x = table
discordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND < c)
  # for each matrix[r, c]
  mat.ll <- function(r, c)
  { 
    ll <- x[(r.x > r) & (c.x < c)]
    sum(ll)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.ll, r = r.x, c = c.x))
}

# Calculate Pairs tied on Rows
# cycle through each row of x and multiply by
# sum(x elements to the right of x[r, c])
# x = table
ties.row <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  total.pairs <- 0

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  for (r in 1:rows)
  {
    for (c in 1:(cols - 1))
    {
      total.pairs <- total.pairs + (x[r, c] * sum(x[r, (c + 1):cols]))
    }
  }

  total.pairs
}

# Calculate Pairs tied on Columns
# cycle through each col of x and multiply by
# sum(x elements below x[r, c])
# x = table
ties.col <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  total.pairs <- 0

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  for (c in 1:cols)
  {
    for (r in 1:(rows - 1))
    {
      total.pairs <- total.pairs + (x[r, c] * sum(x[(r + 1):rows, c]))
    }
  }

  total.pairs
}


calc.gamma <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)

  gamma <- (c - d) / (c + d)

  gamma
}

gamm= calc.gamma(y)
gamm


# ------------------------------------------------------------------Calculate Kendall's Tau-b
# x = table

y=matrix(data=c(2,2,10,5,5,5,10,5,2), nrow=3, ncol=3)
colnames(y) <- c("*", "**","***")
rownames(y) <- c("***","**","*")
y


calc.KTb <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)

  # An alternative computation is:
  #Tr <- ties.row(x)
  #Tc <- ties.col(x)
  #KTb <- (c - d) / sqrt((c + d + Tc) * (c + d + Tr))

  # The "preferred" computation is:
  n <- sum(x)
  SumR <- rowSums(x)
  SumC <- colSums(x)

  KTb <- (2 * (c - d)) / sqrt(((n ^ 2) - (sum(SumR ^ 2))) * ((n ^ 2) - (sum(SumC ^ 2))))

  KTb
}

KTB = calc.KTb(y)
KTB


#-------------------------------------------------------- Calculate Kendall-Stuart Tau-c
# x = table
calc.KSTc <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)
  m <- min(dim(x))
  n <- sum(x)

  KSTc <- (m * 2 * (c - d)) / ((n ^ 2) * (m - 1))

  KSTc
}

#------------------------------------------------------------- Calculate Somers' d
# Return 3 values:
# 1. Sd C~R
# 2. Sd R~C
# 3. Sd Symmetric
# x = table

y=matrix(data=c(923,756,1022,361,1201,7417), nrow=3, ncol=2)
colnames(y) <- c("jong", "oud")
rownames(y) <- c("vaak","soms","nooi")
y


calc.Sd <- function(x)
{
  #-x <- matrix(as.numeric(x), dim(x))
  c <- concordant(x)
  d <- discordant(x)
  n <- sum(x)
  SumR <- rowSums(x)
  SumC <- colSums(x)

  Sd.CR <- (2 * (c - d)) / ((n ^ 2) - (sum(SumR ^ 2)))
  Sd.RC <- (2 * (c - d)) / ((n ^ 2) - (sum(SumC ^ 2)))
  Sd.S <- (2 * (c - d)) / ((n ^ 2) - (((sum(SumR ^ 2)) + (sum(SumC ^ 2))) / 2))

  Sdlist <- list(Sd.CR, Sd.RC, Sd.S)
  names(Sdlist) <- c("Sd.CR", "Sd.RC", "Sd.S")

  Sdlist
}

sd = calc.Sd(y)
sd


#-------------------------------------------- spearman rho

y=matrix(data=c(200,151,129,100,88,75,70,55), nrow=8, ncol=1)
colnames(y) <- c("vrienden")
rownames(y) <- c("thoams", "Sanne","david","tijn", "Anne", "Merel", "Emma", "Thuis")
x=matrix(data=c(190,220,160,90,120,0,60,30), nrow=8, ncol=1)
colnames(y) <- c("tijd")
rownames(y) <- c("thoams", "Sanne","david","tijn", "Anne", "Merel", "Emma", "Thuis")

library(Hmisc);
y1 =rank(y)
x1 =rank(x)
rcorr(x1,y1,type="spearman");


## Ties
y=matrix(data=c(270,240,210,180,160,120,90,60,30,0), nrow=10, ncol=1)
colnames(y) <- c("vrienden")
rownames(y) <- c("thoams", "Sanne","david","tijn", "Anne", "Merel", "Emma", "Thuis","abel", "connie")
x=matrix(data=c(140,140,115,130,110,110,100,104,90,110), nrow=10, ncol=1)
colnames(y) <- c("tijd")
rownames(y) <- c("thoams", "Sanne","david","tijn", "Anne", "Merel", "Emma", "Thuis","abel","connie")
y1 =rank(-y)
x1 =rank(-x)
y1
x1

rcorr(x1,y1,type="spearman");


#-----------------------------------------bereken correlatiecooefcient

y=matrix(data=c(2, 4, 1.5, 2, 3, 7, 3, 8, 8, 6. ), nrow=5, ncol=2)
colnames(y) <- c("x", "y")
rownames(y) <- c("A","B","C","D")
y
cov(y)
cor(y)


#----------------------------------------lineaire regressie
y = c(30, 70, 80) 
x = c(10,30,50)

eruption.lm = lm(y ~ x) 
eruption.lm
summary(eruption.lm)
library(psych);
describe(y)
#--y= y-mean(y)
#--y
plot(x,y)
library(ggplot2)
#--qplot(y,x)
abline(eruption.lm)
abline(h = mean(y), v = , col = "red")

#--------------------------------- Cohens Kappa
rater1 = c(1,2,3,4,5,6,7,8,9) # rater one's ratings
rater2 = c(1,3,1,6,1,5,5,6,7) # rater one's ratings
cohen.kappa(x=cbind(rater1,rater2))
 
#data matrix taken from Cohen
cohen <- matrix(c(
3, 0, 0,
1, 3, 1,
0, 0, 2),ncol=3,byrow=TRUE)

cohen.kappa(cohen)
 
#cohen.weights  weight differences
cohen.weights <- matrix(c(
0,1,3,
1,0,6,
3,6,0),ncol=3)
 
 
cohen.kappa(cohen,cohen.weights,n.obs=200)



#--------------------------------- Cronbach alpha
#--- werkt niet ????

library(psych)
data(expsy)
cronbach(expsy[,1:10])  ## not good because item 2 is reversed (1 is high and 4 is low)
cronbach(cbind(expsy[,c(1,3:10)],-1*expsy[,2]))  ## better


m <- matrix(c(1:10, 11:20), nrow = 10, ncol = 2)
m
apply(m, 1, mean)
m
apply(m, 1:2, function(x) x/2)


dnorm(4, mean = 10, sd = 6, log = FALSE)
dnorm(4, mean = 0, sd = 6, log = FALSE)
0.8*dnorm(4, mean = 10, sd = 6, log = FALSE)/( 0.8*dnorm(4, mean = 10, sd = 6, log = FALSE) +  0.2*dnorm(4, mean = 0, sd = 6, log = FALSE))
