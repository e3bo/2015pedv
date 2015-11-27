                                        # Adapted from example in DiceKrigin::km

#' ## Check for consistency of nugget estimates

d <- 3          # problem dimension
n <- 75                 # size of the experimental design
design <- matrix(runif(n*d), n, d)

covtype <- "matern5_2"          
theta <- c(0.3, 0.5, 1)         # the parameters to be found by estimation
sigma <- 2
nugget <- 3  # choose a numeric value if you want to estimate nugget 
nugget.estim <- TRUE # choose TRUE if you want to estimate it

n.simu <- 30            # number of simulations
sigma2.estimate <- nugget.estimate <- mu.estimate <- matrix(0, n.simu, 1)
coef.estimate <- matrix(0, n.simu, length(theta))

model <- km(~1, design=data.frame(design), response=rep(0,n), covtype=covtype, 
            coef.trend=0, coef.cov=theta, coef.var=sigma^2, nugget=nugget)
y <- simulate(model, nsim=n.simu)

for (i in 1:n.simu) {
        # parameter estimation: tune the optimizer by changing optim.method, control
        model.estimate <- km(~1, design=data.frame(design), response=data.frame(y=y[i,]), 
        covtype=covtype, optim.method="BFGS", control=list(pop.size=50, trace=FALSE), 
        nugget.estim=nugget.estim) 
        
        # store results
        coef.estimate[i,] <- covparam2vect(model.estimate@covariance)
        sigma2.estimate[i] <- model.estimate@covariance@sd2
        mu.estimate[i] <- model.estimate@trend.coef
        if (nugget.estim) nugget.estimate[i] <- model.estimate@covariance@nugget
}

# comparison true values / estimation
cat("\nResults with ", n, "design points, 
    obtained with ", n.simu, "simulations\n\n",
    "Median of covar. coef. estimates: ", apply(coef.estimate, 2, median), "\n",
    "Median of trend  coef. estimates: ", median(mu.estimate), "\n", 
    "Mean of the var. coef. estimates: ", mean(sigma2.estimate), "\n")
if (nugget.estim) cat("\nMean of the nugget effect estimates: ", 
                      mean(nugget.estimate), "\n")

# one figure for this specific example - to be adapted
split.screen(c(2,1))        # split display into two screens
split.screen(c(1,3), screen = 2) # now split the bottom half into 3

screen(1)
boxplot(coef.estimate[,1], coef.estimate[,2], coef.estimate[,3], 
        names=c("theta1", "theta2", "theta3"))
abline(h=theta, col="red")
fig.title <- paste("Empirical law of the parameter estimates 
                    (n=", n , ", n.simu=", n.simu, ")", sep="")
title(fig.title)

screen(3)
boxplot(mu.estimate, xlab="mu")
abline(h=0, col="red")

screen(4)
boxplot(sigma2.estimate, xlab="sigma2")
abline(h=sigma^2, col="red")

screen(5)
boxplot(nugget.estimate, xlab="nugget")
abline(h=nugget, col="red")

close.screen(all = TRUE)  
     
#' ## Now check that nugget parameter can estimate noise variance

noise.var <- 4
model <- km(~1, design=data.frame(design), response=rep(0,n), covtype=covtype, 
            coef.trend=0, coef.cov=theta, coef.var=sigma^2, nugget=0)
y <- simulate(model, nsim=n.simu*noise.var)

y.tilde <- y + rnorm(y, sd=sqrt(noise.var))

for (i in 1:(n.simu*noise.var)) {
        # parameter estimation: tune the optimizer by changing optim.method, control
        model.estimate <- km(~1, design=data.frame(design), response=data.frame(y=y.tilde[i,]), 
        covtype=covtype, optim.method="BFGS", control=list(pop.size=50, trace=FALSE), 
        nugget.estim=nugget.estim) 
        
        # store results
        coef.estimate[i,] <- covparam2vect(model.estimate@covariance)
        sigma2.estimate[i] <- model.estimate@covariance@sd2
        mu.estimate[i] <- model.estimate@trend.coef
        if (nugget.estim) nugget.estimate[i] <- model.estimate@covariance@nugget

        
}

# comparison true values / estimation
cat("\nResults with ", n, "design points, 
    obtained with ", n.simu, "simulations\n\n",
    "Median of covar. coef. estimates: ", apply(coef.estimate, 2, median), "\n",
    "Median of trend  coef. estimates: ", median(mu.estimate), "\n", 
    "Mean of the var. coef. estimates: ", mean(sigma2.estimate), "\n")
if (nugget.estim) cat("\nMean of the nugget effect estimates: ", 
                      mean(nugget.estimate), "\n")

# one figure for this specific example - to be adapted
split.screen(c(2,1))        # split display into two screens
split.screen(c(1,3), screen = 2) # now split the bottom half into 3

screen(1)
boxplot(coef.estimate[,1], coef.estimate[,2], coef.estimate[,3], 
        names=c("theta1", "theta2", "theta3"))
abline(h=theta, col="red")
fig.title <- paste("Empirical law of the parameter estimates 
                    (n=", n , ", n.simu=", n.simu, ")", sep="")
title(fig.title)

screen(3)
boxplot(mu.estimate, xlab="mu")
abline(h=0, col="red")

screen(4)
boxplot(sigma2.estimate, xlab="sigma2")
abline(h=sigma^2, col="red")

screen(5)
boxplot(nugget.estimate, xlab="noise.variance")
abline(h=noise.var, col="red")

close.screen(all = TRUE)  

#' ## Now check predictions of y given y.tilde

n <- 200                 # size of the experimental design
design <- matrix(runif(n*d), n, d)


noise.var <- 4
model <- km(~1, design=data.frame(design), response=rep(0,n), covtype=covtype, 
            coef.trend=0, coef.cov=theta, coef.var=sigma^2, nugget=0)
n.simu <- 200
y <- simulate(model, nsim=n.simu)
y.tilde <- y + rnorm(y, sd=sqrt(noise.var))
y.hat <- matrix(NA, nrow=nrow(y.tilde), ncol=ncol(y.tilde))

for (i in 1:(n.simu)) {
    pred.model <- km(~1, design=data.frame(design), response=data.frame(y.tilde[i,]), covtype=covtype, 
                     coef.trend=0, coef.cov=theta, coef.var=sigma^2, noise.var=rep(noise.var, len=ncol(y.tilde)))
    y.hat[i, ] <- predict(pred.model, newdata=data.frame(design), type='SK')$mean
}
res <- y.hat - y.tilde
mean(res)
var(as.numeric(res))
qqnorm(res)
