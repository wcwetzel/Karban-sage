# Analysis of Karban's damage experiment
# 25 May 2012

library(gstat)
library(sp)
library(lme4)
library(bbmle)
library(spdep)
#library(rgl)
library(arm)
library(rethinking)
library(RANN)
library(pgirmess)

d = read.csv("/Users/will/Documents/Analysis for colleagues/Karban - Damage Neighbors/Sage99plants2011.csv")

# Center X AND Y!!!
d$y = d$y - mean(d$y)
d$x = d$x - mean(d$x)

coordinates(d) = c('x', 'y')
d$x = coordinates(d)[,1]
d$y = coordinates(d)[,2]

# Spatial blocks with 3 blocks
# block 1 is lower right, 2 is upper right, 3 is left
d$spblock[d$receiver > 73] = 1 
d$spblock[is.na(d$spblock)][d$x[is.na(d$spblock)] > 0] = 2
d$spblock[is.na(d$spblock)][d$x[is.na(d$spblock)] < 0] = 3

plot(y ~ x, data=d, type='n')
text(y ~ x, data=d, labels=d$spblock)

# Spatial blocks with 2 blocks
# block 1 is right side (lower herb), 2 is left side (high herb)
d$spblock2[d$receiver > 73] = 1 
d$spblock2[is.na(d$spblock2)][d$x[is.na(d$spblock2)] > 0] = 1
d$spblock2[is.na(d$spblock2)][d$x[is.na(d$spblock2)] < 0] = 2

plot(y ~ x, data=d, type='n')
text(y~x, data=d, labels=d$spblock2)

hist(d$damage[d$spblock2==1], xlim=c(0,100), border='blue')
hist(d$damage[d$spblock2==2], xlim=c(0,100), add=TRUE, lty=2, 
	breaks=9, border='red')


d$healthy = 100 - d$damage
d$damage.prop = d$damage/100

d$asdamage = asin(sqrt(d$damage.prop))
hist(d$asdamage)
hist(log(d$damage.prop))
hist(log(d$damage))


# some plotting to explore spatial patterns
plot(y ~ x, data=d, cex=log(d$damage), asp=1)
plot(y ~ x, data=d, cex=d$damage/50, asp=1)

plot(y ~ x, data=d, cex=d$r+1, asp=1)

plot3d(d$damage ~ d$x + d$y, type='s', size=0.5)
plot3d(d$damage ~ d$x + d$y, type='h', add=TRUE)

plot3d(d$damage ~ d$x + d$r, type='s', size=0.5)
plot3d(d$damage ~ d$x + d$r, type='h', add=TRUE)

plot3d(d$r ~ d$x + d$y, type='s', size=0.5)
plot3d(d$r ~ d$x + d$y, type='h', add=TRUE)

plot(damage ~ x, data=d, pch=as.character(d$spblock2))
plot(damage ~ r, data=d, pch=as.character(d$spblock2))
plot(asdamage ~ r, data=d)
plot(damage ~ spblock2, data=d)

spplot(d, 'damage', do.log=TRUE)
bubble(d, 'damage')

# distance matrix
dists = as.matrix(dist(coordinates(d)))

# autocovariate of damage
d$autodam = numeric(65)
for(i in 1:nrow(d)){
	d$autodam[i] = sum( d$damage[-i] / (-0.5 + dists[i,-i]) )
}


plot(damage ~ autodam, data=d)
corr(cbind(d$damage, d$autodam))
summary(lm(damage ~ autodam, data=d))


# Maybe I don't care about spatial autocorrelation.
# Maybe there is no spatial autocorrelation once I take out the
	# effect of the x dimension?
# So maybe I just need to include x in a model w/o autocorr.

k1 = lm(log(damage) ~ I(x^2) + I(y^2) + I(x*y) + x + y, data=d)

k2 = glm(cbind(d$damage, d$healthy) ~ 1, family='binomial')

k3 = glm(cbind(d$damage, d$healthy) ~ 
	I(x^2) + I(y^2) + I(x*y) + x + y, data=d, family='binomial')

k4 = glmer(cbind(d$damage, d$healthy) ~ 
	I(x^2) + I(y^2) + I(x*y) + x + y + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

k5 = glmer(cbind(d$damage, d$healthy) ~ 
	I(x^2) + I(y^2) + I(x*y) + x + y + r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

k6 = glmer(cbind(d$damage, d$healthy) ~ 
	x + r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

k7 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

k8 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (x|receiver), data=d, 
	family='binomial', REML=FALSE)


k9 = glmer(cbind(d$damage, d$healthy) ~ 
	(x|receiver), data=d, 
	family='binomial', REML=FALSE)

k10 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (0+x|receiver), data=d, 
	family='binomial', REML=FALSE)
	
k11 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

k12 = glmer(cbind(d$damage, d$healthy) ~ 
	(1|receiver), data=d, 
	family='binomial', REML=FALSE)

k13 = glmer(cbind(d$damage, d$healthy) ~
	(1|receiver) + (x), data=d,
	family='binomial', REML=FALSE)

local = localG(d$damage, rlistW)

k14 = glmer(cbind(d$damage, d$healthy) ~
	as.vector(local) + (1|receiver), data=d,
	family='binomial', REML=FALSE)

k15 = glmer(cbind(d$damage, d$healthy) ~
	r + as.vector(local) + (1|receiver), data=d,
	family='binomial', REML=FALSE)


k16 = glmer(cbind(d$damage, d$healthy) ~ 
	(1 + x|receiver), data=d, 
	family='binomial', REML=FALSE)

k17 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1 + x|receiver), data=d, 
	family='binomial', REML=FALSE)


k18 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k19 = glmer(cbind(d$damage, d$healthy) ~ 
	x + r + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k20 = glmer(cbind(d$damage, d$healthy) ~ 
	(x|receiver) + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k21 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (x|receiver) + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)


k24 = glmer(cbind(d$damage, d$healthy) ~ 
	(1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k25 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k26 = glmer(cbind(d$damage, d$healthy) ~ 
	(1 + autodam + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k27 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1 + autodam + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k28 = glmer(cbind(d$damage, d$healthy) ~ 
	autodam + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

k29 = glmer(cbind(d$damage, d$healthy) ~ 
	r + autodam + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)





AICctab(k6, k7, k11, k12, nobs=65)
anova(k6, k7)

extractDIC(k6) - extractDIC(k6)
extractDIC(k7) - extractDIC(k6)

anova(k14, k15)

AICctab(k8, k9, nobs=65)
anova(k8, k9)

AICctab(k11, k12, nobs=65)
anova(k11, k12)


AICctab(k16, k17, nobs=65)
anova(k16, k17)

AICctab(k18, k19, nobs=65)
AICtab(k18, k19)
anova(k18, k19)

AICctab(k20, k21, nobs=65)
anova(k20, k21)


AICctab(k24, k25, nobs=65)
anova(k24, k25)

AICctab(k26, k27, nobs=65)
anova(k26,k27)

AICctab(k28, k29, nobs=65)
anova(k28,k29)






AICctab(k4,k5, k6, k7, k8, k9, k10, nobs=nrow(d))

plot(damage ~ r, data=d[-4,])
plot(resid(k5) ~ r, data=d)

plot(y ~ x, data=d, cex= (resid(k11)+1)/ (max(resid(k11) + 1)) *3 )



#############

# random and fixed effect of x, McElreath's suggestion
k22 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1 + x|receiver), data=d, 
	family='binomial', REML=FALSE)
k23 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + (1 + x|receiver), data=d, 
	family='binomial', REML=FALSE)
k23.5 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + r:x + (1 + x|receiver), data=d, 
	family='binomial', REML=FALSE)
extractDIC(k22) - extractDIC(k23.5)

# fixed effect of x
k30 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
k31 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
k31.5 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + r:x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

# random slopes for r and random and fixed x
k32 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1 + x + r|receiver), data=d, 
	family='binomial', REML=FALSE)
k33 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + (1 + x + r|receiver), data=d, 
	family='binomial', REML=FALSE)

# random slopes for r and fixed effect for x
k34 = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)
k35 = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + (1 + r|receiver), data=d, 
	family='binomial', REML=FALSE)

# no effect of x
k36 = glmer(cbind(d$damage, d$healthy) ~ 
	(1|receiver), data=d, 
	family='binomial', REML=FALSE)
k37 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)


AICctab(k22, k23, k23.5, nobs=65)
BICtab(k22, k23, nobs=65)
anova(k22,k23)
anova(k22, k23.5)
extractDIC(k22)
extractDIC(k23)


AICctab(k30, k31, k31.5, nobs=65)
BICtab(k22, k23, nobs=65)
anova(k30,k31)
extractDIC(k30)
extractDIC(k31)


AICctab(k32, k33, nobs=65)
anova(k32,k33)
extractDIC(k32)
extractDIC(k33)

AICctab(k34, k35, nobs=65)
anova(k34,k35)
extractDIC(k34)
extractDIC(k35)

BICtab(k32,k33, k22,k23, nobs=65)
extractDIC(k33)
extractDIC(k23)

AICctab(k36, k37, nobs=65)
anova(k36,k37)
extractDIC(k36)
extractDIC(k37)

AICctab(k22, k23, k23.5,k30, k31, k31.5, nobs=65)

#---------- final models with x instead of sp block --------------------------#

# fixed effect of x, no random effect of x
x0 = glmer(cbind(d$damage, d$healthy) ~ 
	(1|receiver), data=d, 
	family='binomial', REML=FALSE)
xx = glmer(cbind(d$damage, d$healthy) ~ 
	x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
xr = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
xxr = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
xxrint = glmer(cbind(d$damage, d$healthy) ~ 
	r + x + r:x + (1|receiver), data=d, 
	family='binomial', REML=FALSE)


AICctab(x0, xx, xr, xxr, xxrint, weights=TRUE, nobs=65)




#------------ models with 3 spatial blocks (spblock) ------------#

d$spblock = as.factor(d$spblock)
d$receiver = as.factor(d$receiver)

sp0 = glmer(cbind(d$damage, d$healthy) ~ 
	(1| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0, sp1, nobs=65)
anova(sp0, sp1)

sp0r = glmer(cbind(d$damage, d$healthy) ~ 
	(1 + r| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1r = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1 + r| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0r, sp1r, nobs=65)
anova(sp0r, sp1r)

sp0int = glmer(cbind(d$damage, d$healthy) ~ 
	(1| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)
sp0.5int = glmer(cbind(d$damage, d$healthy) ~ 
	spblock + (1| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1int = glmer(cbind(d$damage, d$healthy) ~ 
	r + spblock + r:spblock + (1| spblock / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0int, sp1int, sp0.5int, nobs=65)
anova(sp0.5int, sp1int)



#------ explore effects of spatial blocks  -----------------#

plot(damage.prop ~ r, data=d, pch=as.character(d$spblock),
	col=d$spblock)
lines(loess.smooth(d$r[d$spblock==1], d$damage.prop[d$spblock==1]))
lines(loess.smooth(d$r[d$spblock==2], d$damage.prop[d$spblock==2]), 
	col='red')
lines(loess.smooth(d$r[d$spblock==3], d$damage.prop[d$spblock==3]),
	col='green')
abline(lm(damage.prop ~ r, data=d, subset=d$spblock==1), lty=2)
abline(lm(damage.prop ~ r, data=d, subset=d$spblock==2), lty=2,
	col='red')
abline(lm(damage.prop ~ r, data=d, subset=d$spblock==3), lty=2,
	col='green')

summary(lm(damage.prop ~ r, data=d, subset=d$spblock==1))

#------------ models with 2 spatial blocks (spblock) ------------#

d$spblock2 = as.factor(d$spblock2)
d$receiver = as.factor(d$receiver)

sp0 = glmer(cbind(d$damage, d$healthy) ~ 
	(1| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1 = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0, sp1, nobs=65)
anova(sp0, sp1)

sp0r = glmer(cbind(d$damage, d$healthy) ~ 
	(1 + r| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1r = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1 + r| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0r, sp1r, nobs=65)
anova(sp0r, sp1r)

sp0int = glmer(cbind(d$damage, d$healthy) ~ 
	(1| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)
sp0.5int = glmer(cbind(d$damage, d$healthy) ~ 
	spblock2 + (1| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)
sp1int = glmer(cbind(d$damage, d$healthy) ~ 
	r + spblock2 + r:spblock2 + (1| spblock2 / receiver), data=d, 
	family='binomial', REML=FALSE)

AICctab(sp0int, sp1int, sp0.5int, nobs=65)
anova(sp0.5int, sp1int)

sp0intx = glmer(cbind(d$damage, d$healthy) ~ 
	(1|receiver), data=d, 
	family='binomial', REML=FALSE)
sp0.5intx = glmer(cbind(d$damage, d$healthy) ~ 
	spblock2 + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
spnoint2x = glmer(cbind(d$damage, d$healthy) ~ 
	r + spblock2 + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
spnoint1x = glmer(cbind(d$damage, d$healthy) ~ 
	r + (1|receiver), data=d, 
	family='binomial', REML=FALSE)
sp1intx = glmer(cbind(d$damage, d$healthy) ~ 
	r + spblock2 + r:spblock2 + (1|receiver), data=d, 
	family='binomial', REML=FALSE)

sp1intxmeans = glmer(cbind(d$damage, d$healthy) ~ 
	r + spblock2 + r:spblock2 + (1|receiver) - 1, data=d, 
	family='binomial', REML=FALSE)
	
AICctab(sp0intx, sp1intx, sp0.5intx, spnoint1x, spnoint2x, nobs=65, weights=TRUE)
AICtab(sp0intx, sp1intx, sp0.5intx, spnoint1x, spnoint2x, weights=TRUE)
extractDIC(sp1intx) - extractDIC(sp1intx)
extractDIC(sp0.5intx) - extractDIC(sp1intx)
extractDIC(spnoint2x) - extractDIC(sp1intx)
extractDIC(spnoint1x) - extractDIC(sp1intx)
extractDIC(sp0intx) - extractDIC(sp1intx)




AICctab(sp0intx, sp1intx, sp0.5intx, spnoint1x, spnoint2x, x0, xx, xr, xxr, xxrint, nobs=65, weights=TRUE)

AICctab(sp0intx, sp1intx, sp0.5intx, sp0int, sp1int, sp0.5int, nobs=65, weights=TRUE)
BICtab(sp0intx, sp1intx, sp0.5intx, sp0int, sp1int, sp0.5int, nobs=65)

anova(sp0.5intx, sp1intx)
anova(sp0intx, sp1intx, sp0.5intx, spnoint1x, spnoint2x, nobs=65)


extractDIC(sp0.5intx) - extractDIC(sp1intx)

AICctab(sp0intx, sp1intx, sp0.5intx, spnoint1x, spnoint2x, k22, k23, k23.5,
	k30, k31, k31.5, nobs=65, weights=TRUE)



# simulate counts using ran eff estimate from sp1intx
# does it match the observed distribution?
par(mfrow=c(3,2))
hist(d$damage, xlim=c(0,100))
hist(rbinom(65, 100, mean(d$damage/100)), xlim=c(0,100))
herbsim = rbinom(65, 100, 
	plogis( rnorm(65, mean = qlogis(mean(d$damage/100)), sd = 1.0404) ) )
hist(herbsim, xlim=c(0,100))
herbsim = rbinom(65, 100, 
	plogis( rnorm(65, mean = qlogis(mean(d$damage/100)), sd = 1.0404) ) )
hist(herbsim, xlim=c(0,100))
herbsim = rbinom(65, 100, 
	plogis( rnorm(65, mean = qlogis(mean(d$damage/100)), sd = 1.0404) ) )
hist(herbsim, xlim=c(0,100))
herbsim = rbinom(65, 100, 
	plogis( rnorm(65, mean = qlogis(mean(d$damage/100)), sd = 1.0404) ) )
hist(herbsim, xlim=c(0,100))
# looks close enough, mode could be closer to zero.


#------ explore effects of 2 spatial blocks  -----------------#

plot(damage.prop ~ r, data=d, pch=as.character(d$spblock2),
	col=d$spblock2)
lines(loess.smooth(d$r[d$spblock2==1], d$damage.prop[d$spblock2==1]))
lines(loess.smooth(d$r[d$spblock2==2], d$damage.prop[d$spblock2==2]), 
	col='red')
abline(lm(damage.prop ~ r, data=d, subset=d$spblock2==1), lty=2)
abline(lm(damage.prop ~ r, data=d, subset=d$spblock2==2), lty=2,
	col='red')


summary(lm(damage.prop ~ r, data=d, subset=d$spblock2==1))
summary(lm(damage.prop ~ r, data=d, subset=d$spblock2==2))

#----------- plot prediction of spblock2 interaction x model-------#

post.sp1intx = sample.naive.posterior(sp1intx)

newr = seq(-0.6, 0.8, length=100)
#newblock = sort(rep(1:2, 100))

pred.sp1intx.block1 = sapply(newr, 
	function(z) plogis(mean(post.sp1intx[,1] + post.sp1intx[,2] * z)))
ci.sp1intx.block1 = sapply(newr,
	function(z) plogis(HPDI(post.sp1intx[,1] + post.sp1intx[,2] * z)))

pred.sp1intx.block2 = sapply(newr, 
	function(z) plogis(mean(post.sp1intx[,1] + post.sp1intx[,2] * z +
	post.sp1intx[,3] + post.sp1intx[,4] * z)))
ci.sp1intx.block2 = sapply(newr,
	function(z) plogis(HPDI(post.sp1intx[,1] + post.sp1intx[,2] * z +
	post.sp1intx[,3] + post.sp1intx[,4] * z)))

plot(damage.prop ~ r, data=d, pch=20,col=d$spblock2, 
	ylab='Rate of herbivory on receiver',
	xlab='Relatedness between volatile donor and receiver',
	las=1)
lines(newr, pred.sp1intx.block1)
lines(newr, ci.sp1intx.block1[1,], lty=2)
lines(newr, ci.sp1intx.block1[2,], lty=2)
lines(newr, pred.sp1intx.block2, col='red')
lines(newr, ci.sp1intx.block2[1,], lty=2, col='red')
lines(newr, ci.sp1intx.block2[2,], lty=2, col='red')
legend(-0.15, 1, pch=20, lty=1, col=c(1,2),
	legend=c('East block, low herbivory', 'West block, high herbivory'))


#----------- exploratory variogram analysis ----------------#

hscat(d$damage ~ 1, d, (0:9) * 5)
hscat(log(d$damage) ~ 1, d, (0:9) * 5)
plot(variogram(d$damage ~ 1, d, cloud=FALSE))
plot(variogram(log(d$damage) ~ 1, d, cloud=FALSE))


plot(d$x, d$y, cex=resid(k11)+1)
plot(resid(k11)+1 ~ d$x)

#------ neighbors and weight lists -----------#
# turn x-y coordinates into a square distance matrix
dists = as.matrix(dist(cbind(d$x,d$y), diag=TRUE, upper=TRUE))

### inverse distance weights ###
curve(1/(1+x), 0, max(dists))
rW = 1/(1+dists)
diag(rW) = 0 # ack, I think I've always had diag of weights matrix
# as 1, but I think it needs to be zero.

#plot(nelistW, coordinates(d))
rlistW = mat2listw(rW, row.names=d$receiver)
summary(rlistW)

### negative exponential ###
# convert distance matrix into weights matrix
# with negative exp decline in weight with distance
curve(exp(-x * 0.5), 0, max(dists))
nexpW = exp(-dists * 0.5)
diag(nexpW) = 0
# mat2listw converts a square spatial weights matrix into 
# a spdep weights list object
nelistW = mat2listw(nexpW, row.names=d$receiver)
summary(nelistW)

### nearest neighbors ###
nb2listw

#------------- Moran's I -----------------------------#
# sim data with no spat autocorr
uncorr = rbinom(65, 100, mean(d$damage/100))
moran.test(d$damage, listw=nelistW)
moranI.damage.rW = moran.test(d$damage, listw=rlistW)
moranI.damage.rW$statistic

moran.test(uncorr, listw=nelistW)
moran.test(uncorr, listw=rlistW)

# bootstrap it
nboots = 999
moranIboot = numeric(nboots)
for(i in 1:nboots){
	bootdata = sample(d$damage, length(d$damage), replace=TRUE)
	mt = moran.test(bootdata, listw = rlistW)
	moranIboot[i] = mt$statistic
}

hist(moranIboot, xlim=c(min(c(moranIboot, moranI.damage.rW$statistic)),
	max(c(moranIboot, moranI.damage.rW$statistic))))
abline(v=moranI.damage.rW$statistic)
ci = HPDI(moranIboot)
abline(v=ci, lty=2, col='grey')

#------------------ Geary's C ------------------------#
gearyC.damage.rW = geary.test(d$damage, listw=rlistW)
geary.test(d$damage, listw=nelistW)
geary.test(uncorr, listw=nelistW)

# bootstrap it
nboots = 99
gearyCboot = numeric(nboots)
for(i in 1:nboots){
	bootdata = sample(d$damage, length(d$damage), replace=TRUE)
	mt = geary.test(bootdata, listw = rlistW)
	gearyCboot[i] = mt$statistic
}

hist(gearyCboot, xlim=c(min(c(gearyCboot, gearyC.damage.rW$statistic)),
	max(c(gearyCboot, gearyC.damage.rW$statistic))))
abline(v=gearyC.damage.rW$statistic)
ci = HPDI(gearyCboot)
abline(v=ci, lty=2, col='grey')

#------------------ Geary's C  log(d$damage) ------------------------#
gearyC.log.damage.rW = geary.test(log(d$damage), listw=rlistW)
geary.test(log(d$damage), listw=nelistW)
geary.test(log(uncorr), listw=rlistW)

# bootstrap it
nboots = 99
gearyCboot = numeric(nboots)
for(i in 1:nboots){
	bootdata = sample(log(d$damage), length(d$damage), replace=TRUE)
	mt = geary.test(bootdata, listw = rlistW)
	gearyCboot[i] = mt$statistic
}

hist(gearyCboot, xlim=c(min(c(gearyCboot, gearyC.log.damage.rW$statistic)),
	max(c(gearyCboot, gearyC.log.damage.rW$statistic))))
abline(v=gearyC.log.damage.rW$statistic)
ci = HPDI(gearyCboot)
abline(v=ci, lty=2, col='grey')


#------------------ Global G ------------------------#
globalG.damage.rW = globalG.test(d$damage, listw=rlistW)
geary.test(d$damage, listw=nelistW)
geary.test(uncorr, listw=nelistW)

# bootstrap it
nboots = 99
globalGboot = numeric(nboots)
for(i in 1:nboots){
	bootdata = sample(d$damage, length(d$damage), replace=TRUE)
	mt = globalG.test(bootdata, listw = rlistW)
	globalGboot[i] = mt$statistic
}

hist(globalGboot, xlim=c(min(c(globalGboot, globalG.damage.rW$statistic)),
	max(c(globalGboot, globalG.damage.rW$statistic))))
abline(v=globalG.damage.rW$statistic)
ci = HPDI(globalGboot)
abline(v=ci, lty=2, col='grey')

#------------------ Global G log(d$damage) ------------------------#
globalG.damage.rW = globalG.test(log(d$damage), listw=rlistW)
globalG.test(log(d$damage), listw=nelistW)

# bootstrap it
nboots = 99
globalGboot = numeric(nboots)
for(i in 1:nboots){
	bootdata = sample(log(d$damage), length(d$damage), replace=TRUE)
	mt = globalG.test(bootdata, listw = rlistW)
	globalGboot[i] = mt$statistic
}

hist(globalGboot, xlim=c(min(c(globalGboot, globalG.damage.rW$statistic)),
	max(c(globalGboot, globalG.damage.rW$statistic))))
abline(v=globalG.damage.rW$statistic)
ci = HPDI(globalGboot)
abline(v=ci, lty=2, col='grey')

#------------- set up neighbor lists -------------------#
# nearest neighbor for each plant
nblist = knn2nb(knearneigh(coordinates(d)))
# max, min, mean nearest neighbor relationships
max(unlist(nbdists(nblist, coordinates(d))))
min(unlist(nbdists(nblist, coordinates(d))))
mean(unlist(nbdists(nblist, coordinates(d))))

# neighbors by distance
# 10m = smallest distance so that all plants have a neighbor
nb10<-dnearneigh(coordinates(d),0,15,row.names=d$receiver)
plot(nb10, coordinates(d))


#---------------- Correlograms -------------------------#

cor.nb10 = sp.correlogram(nblist, d$damage, order=10, method='I')

correl = correlog(coordinates(d), d$damage)
# log transformed
correllog = correlog(coordinates(d), log(d$damage))

plot(correl[,1], correl[,2], pch=(c(19,19, rep(21,11))), type='b')
abline(h=0, lty=2, col='grey')
# log transofmed
plot(correllog[1:10,1], correllog[1:10,2], pch=(c(19,19, rep(21,8))), type='b',
	xlab='Distance', ylab='Moran I', las=1)
abline(h=0, lty=2, col='grey')
#autocorrelation to about 11 m

# with residuals #
correl.glmm.r.res = correlog(coordinates(d), residuals(k11, type='pearson'))
plot(correl.glmm.r.res[,1], correl.glmm.r.res[,2], pch=(c(19,19, rep(21,11))), type='b')
abline(h=0, lty=2, col='grey')



#-------------- Spline correlogram ---------------------#
library(ncf)
spline.correlog.dam = spline.correlog(x = d$x, y = d$y, z = d$damage, xmax=50)
plot.spline.correlog(spline.correlog.dam)
summary(spline.correlog.dam)

spline.correlog.damprop = spline.correlog(x = d$x, y = d$y, z = d$damage.prop, xmax=50)
plot.spline.correlog(spline.correlog.damprop)
summary(spline.correlog.damprop)

spline.correlog.logdam = spline.correlog(x = d$x, y = d$y, z = log(d$damage), xmax=50)
plot.spline.correlog(spline.correlog.logdam)
summary(spline.correlog.logdam)

spline.correlog.glmm.intercept.res = spline.correlog(x = d$x, y = d$y, 
	z = residuals(k12, type='pearson'), xmax=50)
plot.spline.correlog(spline.correlog.glmm.intercept.res)

spline.correlog.glmm.r.res = spline.correlog(x = d$x, y = d$y, 
	z = residuals(k11, type='pearson'), xmax=50)
plot.spline.correlog(spline.correlog.glmm.r.res)



#--------------- local G ------------------------------#

local = localG(d$damage, rlistW)
plot(coordinates(d), cex=local-min(local), asp=1)





#-------------- Autocovariate regression ----------------#

library(spdep)

autocov = autocov_dist(d$damage, nbs=10, coordinates(d))

plot(autocov ~ d$autodam)
plot(d$damage ~ autocov)
abline(lm(d$damage ~ autocov))

ac0 = glm(cbind(damage, healthy) ~ autocov, data=d, family='binomial')
ac1 = glm(cbind(damage, healthy) ~ r + autocov, data=d, family='binomial')

nac0 = glm(cbind(damage, healthy) ~ 1, data=d, family='binomial')
nac1 = glm(cbind(damage, healthy) ~ r, data=d, family='binomial')

acglmm0 = glmer(cbind(damage, healthy) ~ autocov + (1|receiver), data=d, family='binomial')
acglmm1 = glmer(cbind(damage, healthy) ~ r + autocov + (1|receiver), data=d, family='binomial')


AICctab(ac0, ac1, nobs=65)
AICctab(nac0, nac1, nobs=65)
AICctab(acglmm0, acglmm1, nobs=65)

anova(ac0, ac1)

#----------- Transform and use autoregressive model ------------#

z = log(d$damage.prop)

plot(z ~ d$r)


lm1 = lm(z ~ r, data=d)
plot(residuals(lm1) ~ d$x)
lm.morantest(lm1, rlistW)

m0 = spautolm(z ~ 1, data=d, family='SAR', listw=rlistW)
m1 = spautolm(z ~ r, data=d, family='SAR', listw=rlistW)
AICtab(m0,m1)

m0 = spautolm(z ~ 1, data=d, family='CAR', listw=rlistW)
m1 = spautolm(z ~ r, data=d, family='CAR', listw=rlistW)
AICtab(m0,m1)





#-------------------- GAM on transformed data -----------------#
library(mgcv)
gam0 = gam(z ~ 1)
gam1 = gam(z ~ r, data=d)

AICtab(gam0, gam1)
anova(gam0, gam1)

sg0xy = gam(z ~ 1 + s(x,y), data=d)
sg1xy = gam(z ~ r + s(x,y), data=d)
sg0x = gam(z ~ 1 + s(x), data=d)
sg1x = gam(z ~ r + s(x), data=d)

AICtab(sg0xy, sg1xy, sg0x, sg1x)


#------------------- GAM binomial ----------------------------#

sg0xy = gam(cbind(damage, healthy) ~ 1 + s(x,y), data=d, family=binomial)
sg1xy = gam(cbind(damage, healthy) ~ r + s(x,y), data=d, family=binomial)
sg0x = gam(cbind(damage, healthy) ~ 1 + s(x), data=d, family=binomial)
sg1x = gam(cbind(damage, healthy) ~ r + s(x), data=d, family=binomial)
AICtab(sg0xy, sg1xy, sg0x, sg1x)

#--------------------- lme with spat autocorr -------------------------------#
library(MASS)
library(nlme)
sp1 = corSpatial(form = ~ x + y, type = 'gaussian')
scor = Initialize(sp1, as(d, "data.frame")[, c('x', 'y')], nugget=FALSE)

m0 = lme(z ~ 1, random = ~1|receiver, data=as.data.frame(d), 
	correlation=scor, method='ML')
m1 = lme(z ~ r, random = ~1|receiver, data=as.data.frame(d), 
	correlation=scor, method='ML')

AICtab(m0,m1)


#--------------------- glmmPQL -------------------------------#
library(MASS)
library(nlme)
sp1 = corSpatial(1, form = ~ x + y, type = 'gaussian')
scor = Initialize(sp1, as(d, "data.frame")[, c('x', 'y')], nugget=FALSE)

ddd = as.data.frame(d)

m1PQL = glmmPQL(cbind(damage, healthy) ~ r, random = ~1 | receiver, family=binomial,
	data=as.data.frame(d), correlation=scor)
summary(m1PQL)



#------------------- check for spat autocorr within the 2 blocks --------#


### first for block 1 ###
# turn x-y coordinates into a square distance matrix
dists1 = as.matrix(dist(cbind(d$x[d$spblock2==1],d$y[d$spblock2==1]), diag=TRUE, upper=TRUE))

## inverse distance weights ##
rW1 = 1/(1+dists1)
diag(rW1) = 0
rlistW1 = mat2listw(rW1, row.names=d$receiver[d$spblock2==1])
summary(rlistW1)

## negative exponential ##
nexpW1 = exp(-dists1 * 0.5)
diag(nexpW1) = 0
nelistW1 = mat2listw(nexpW1, row.names=d$receiver[d$spblock2==1])
summary(nelistW1)

## Moran's I ##
# sim data with no spat autocorr
uncorr1 = rbinom(35, 100, mean(d$damage[d$spblock2==1]/100))
moran.test(d$damage[d$spblock2==1], listw=nelistW1)
moran.test(log(d$damage[d$spblock2==1]), listw=nelistW1)
moran.test(d$damage[d$spblock2==1], listw=rlistW1)
moran.test(log(d$damage[d$spblock2==1]), listw=rlistW1)
moran.test(uncorr1, listw=nelistW1)
moran.test(uncorr1, listw=rlistW1)


### second for block 2 ###
dists2 = as.matrix(dist(cbind(d$x[d$spblock2==2],d$y[d$spblock2==2]), diag=TRUE, upper=TRUE))

## inverse distance weights ##
rW2 = 1/(1+dists2)
diag(rW2) = 0
rlistW2 = mat2listw(rW2, row.names=d$receiver[d$spblock2==2])
summary(rlistW1)

## negative exponential ##
nexpW2 = exp(-dists2 * 0.5)
diag(nexpW2) = 0
nelistW2 = mat2listw(nexpW2, row.names=d$receiver[d$spblock2==2])
summary(nelistW1)

## Moran's I ##
# sim data with no spat autocorr
uncorr2 = rbinom(30, 100, mean(d$damage[d$spblock2==2]/100))
moran.test(d$damage[d$spblock2==2], listw=nelistW2)
moran.test(log(d$damage[d$spblock2==2]), listw=nelistW2)
moran.test(d$damage[d$spblock2==2], listw=rlistW2)
moran.test(log(d$damage[d$spblock2==2]), listw=rlistW2)

moran.test(uncorr1, listw=nelistW1)
moran.test(uncorr1, listw=rlistW1)




