library(ncf)

x = runif(60)
y=runif(60)

mu = rep(log(22), 60)

plot(x,y)

dists = as.matrix(dist(cbind(x,y)))
Sigma = exp(10 * -dists^2)
diag(Sigma) = 1
plot(Sigma ~ dists)

lambda = exp(mvrnorm(n=1, mu=mu, Sigma = Sigma))

plot(x,y, cex=lambda/max(lambda) * 2)

z = rpois(length(mu), lambda)

plot(x,y, cex=z/max(z) * 2)


spline.corr = spline.correlog(x=x, y=y, z=z)
plot.spline.correlog(spline.corr)


