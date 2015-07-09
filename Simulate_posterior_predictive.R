# Name: Simulate_posterior_predictive.R
# Desc: Simple script to use create dummy data
# Date: 9/07/2015
# Auth: u.niazi@imperial.ac.uk

####
# load data
dfDat = read.csv(file.choose(), header=T)

dfDat.bk = dfDat
dfDat = dfDat.bk[,c(2,4)]
dfDat = na.omit(dfDat)
fGroups = dfDat$Status
fGroups = factor(fGroups, levels = c('Remote', 'Recent'))

ivDat.rec = dfDat$TNFaTEFF[fGroups == 'Recent']
ivDat.rem = dfDat$TNFaTEFF[fGroups == 'Remote']
ivDat.com = dfDat$TNFaTEFF

col = rainbow(2)
c = col[as.numeric(fGroups)]

i = order(ivDat.com)
d = ivDat.com[i]
c = c[i]

barplot(d, col=c)

# use existing recent data to produce future data
v.dat = var(ivDat.rec)
m.dat = mean(ivDat.rec)
# calculate degrees of freedom
v = sum((ivDat.rec - m.dat)^2)
v = v/length(ivDat.rec)
# calculate inverse chisquar, which is variance
x = rchisq(10000, v)
theta.var = (v*v.dat)/x

# calculate possible values of mean
#s = sample(theta.var, size = 10000, replace = T)
theta.mean = rep(NA, length=10000)
for (i in 1:10000) theta.mean[i] = rnorm(1, mean=m.dat, sample(theta.var, 1))
#theta.mean = rnorm(10000, mean = m.dat, sd = s)

mu.s = sample(theta.mean, size = 1000)
var.s = sample(theta.var, size=1000)
# create matrix for liklihood calculation
lik = matrix(NA, 1000, 1000)
for (m in 1:1000){
  for (v in 1:1000){
    lik[m,v] = prod (dnorm(ivDat.rec, mu.s[m], sqrt(var.s[v]) ) )
  }
}
rownames(lik) = round(mu.s)
colnames(lik) = round(var.s)

which(lik == max(lik), arr.ind = T)
# use the col of matrix with the maximum liklihood for the variance to get the mean vector 
# rejection sampling from prior mean
post.mean = sample(mu.s, size = 10000, replace = T, prob = lik[,555])
post.var = sample(var.s, size=10000, replace = T, prob= lik[857,])

# create new data points by sampling from the data
x.new = rep(NA, 100)
for (i in 1:100) {
  s = NULL
  while(TRUE){ 
    s = rnorm(1, sample(post.mean, 1), sqrt(sample(post.var, size = 1)))
    if (s >=0) break;
  }
  x.new[i] = s
}
