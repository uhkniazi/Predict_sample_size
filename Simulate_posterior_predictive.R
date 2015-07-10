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

# ivDat.rec = dfDat$TNFaTEFF[dfDat$Status == 'Recent']
# ivDat.rem = dfDat$TNFaTEFF[fGroups == 'Remote']
ivDat.com = dfDat$TNFaTEFF
ivDat.rec = dfDat$TNFaTEFF[dfDat$Status == 'Recent']

# bar plot of 2 level data with cutoff for outliers
col = rainbow(2)
c = col[as.numeric(fGroups)]
# order the data for barplot
i = order(ivDat.com)
d = ivDat.com[i]
c = c[i]
barplot(d, col=c)

# calculate the outlier values, using maximum liklihood ratios
lik.dn = dnorm(d, mean = mean(d), sd(d))
lik.dn = lik.dn/max(lik.dn)
outlier = which(lik.dn <= 0.2)

ivDat.com.no.outlier = d[-outlier]
cutoff = mean(ivDat.com.no.outlier) + 2.5*sd(ivDat.com.no.outlier)

# generate groups based on outlier and plot again
fGroups = as.character(dfDat$Status)
i = which(ivDat.com > cutoff)
fGroups[i] = 'Outlier'
fGroups = factor(fGroups)

# create barplot to show outliers
col = rainbow(3)
c = col[as.numeric(fGroups)]

i = order(ivDat.com)
d = ivDat.com[i]
c = c[i]

barplot(d, col=c)
legend('topleft', legend = as.character(levels(fGroups)), fill = col)
abline(h = cutoff)


# use existing recent data to produce future data
#ivDat.rec = ivDat.com

# calculate hyperparameters using old data
v.dat = var(ivDat.rec)
m.dat = mean(ivDat.rec)
# calculate degrees of freedom
v = sum((ivDat.rec - m.dat)^2)
v = v/length(ivDat.rec)
n0 = length(ivDat.rec)
# calculate inverse chisquar, which is variance
x = rchisq(10000, v)
theta.var = (v*v.dat)/x

# calculate possible values of mean
#s = sample(theta.var, size = 10000, replace = T)
theta.mean = rep(NA, length=10000)
for (i in 1:10000) {
  sd.s = sample(theta.var, 1)#/n0
  sd.s = sqrt(sd.s)
  theta.mean[i] = rnorm(1, mean=m.dat, sd.s)
}

# take a sample from this prior distribution, to calculate liklihoods
mu.s = sample(theta.mean, size = 1000)
var.s = sample(theta.var, size=1000)
# create matrix for likelihood calculation for remote data without outliers
lik = matrix(NA, 1000, 1000)
for (m in 1:1000){
  for (v in 1:1000){
    lik[m,v] = prod (dnorm(ivDat.rec[ivDat.rec <= cutoff], mu.s[m], sqrt(var.s[v]) ) )
  }
}
rownames(lik) = round(mu.s)
colnames(lik) = round(var.s)

l = which(lik == max(lik), arr.ind = T)
# use the col of matrix with the maximum liklihood for the variance to get the mean vector 
# rejection sampling from prior mean
post.mean = sample(mu.s, size = 10000, replace = T, prob = lik[,l[1,'col']])
post.var = sample(var.s, size=10000, replace = T, prob= lik[l[1,'row'],])

post.mean.cont = post.mean
post.var.cont = post.var

# repeat for outliers i.e. cases
# create matrix for likelihood calculation for remote data without outliers
lik = matrix(NA, 1000, 1000)
for (m in 1:1000){
  for (v in 1:1000){
    lik[m,v] = prod (dnorm(ivDat.rec[ivDat.rec > cutoff], mu.s[m], sqrt(var.s[v]) ) )
  }
}
rownames(lik) = round(mu.s)
colnames(lik) = round(var.s)

l = which(lik == max(lik), arr.ind = T)
# use the col of matrix with the maximum liklihood for the variance to get the mean vector 
# rejection sampling from prior mean
post.mean = sample(mu.s, size = 10000, replace = T, prob = lik[,l[1,'col']])
post.var = sample(var.s, size=10000, replace = T, prob= lik[l[1,'row'],])

post.mean.case = post.mean
post.var.case = post.var

# create new data points by sampling from the data
iSize = 50
x.cont = rep(NA, iSize)
for (i in 1:iSize) {
  s = NULL
  while(TRUE){ 
    s = rnorm(1, sample(post.mean.cont, 1), sqrt(sample(post.var.cont, size = 1)))
    if (s >=0) break;
  }
  x.cont[i] = s
}

x.case = rep(NA, iSize)
for (i in 1:iSize) {
  s = NULL
  while(TRUE){ 
    s = rnorm(1, sample(post.mean.case, 1), sqrt(sample(post.var.case, size = 1)))
    if (s >=0) break;
  }
  x.case[i] = s
}

# plot the old and new data together
p.old = par(mfrow=c(2,1))
col = rainbow(3)
c = col[as.numeric(fGroups)]

i = order(ivDat.com)
d = ivDat.com[i]
c = c[i]

barplot(d, col=c)
legend('topleft', legend = as.character(levels(fGroups)), fill = col)
abline(h = cutoff)

fGroups = gl(2, k = iSize, labels = c('Recent', 'Outlier'))
col = c('green', 'red')
c = col[as.numeric(fGroups)]
ivDat.sim = c(x.cont, x.case)
i = order(ivDat.sim)
d = ivDat.sim[i]
c = c[i]
barplot(d, col=c)
legend('topleft', legend = as.character(levels(fGroups)), fill = col)
abline(h = cutoff)

delta = abs(mean(x.case) - mean(x.cont))
sd.sim = sd(c(x.case, x.cont))
sig.l = 0.01
iSamSize = 50

pow = rep(NA, length=iSamSize)

for (i in 1:iSamSize) pow[i] = power.t.test(n = i, delta, sd.sim, sig.l, type = 'two.sample')$power

plot(pow, main='Power as function of number of samples', type='l')

