## -----------------------------------------------------------------------------
library(ctsmTMB)
model <- ctsmTMB$new()

## -----------------------------------------------------------------------------
print(model)

## -----------------------------------------------------------------------------
model$addSystem(dx ~ theta * (mu - x) * dt + sigma_x * dw)

## ----eval=FALSE---------------------------------------------------------------
# model$addSystem(dx ~ theta * (mu - x) * dt + sigma_x1 * dw1 + x * sigma_x2 * dw2 + dw3)

## -----------------------------------------------------------------------------
model$addObs(y ~ x)

## -----------------------------------------------------------------------------
model$setVariance(y ~ sigma_y^2 * u)

## -----------------------------------------------------------------------------
model$addInput(u)

## -----------------------------------------------------------------------------
model$setParameter(
  theta   = c(initial = 5,    lower = 0,    upper = 20),
  mu      = c(initial = 0,    lower = -10,  upper = 10),
  sigma_x = c(initial = 1e-1, lower = 1e-5, upper = 5),
  sigma_y = c(initial = 1e-1, lower = 1e-5, upper = 5)
)

## -----------------------------------------------------------------------------
model$setParameter(
  sigma_y  = 0.05
)

## -----------------------------------------------------------------------------
print(model)

## -----------------------------------------------------------------------------
initial.state <- list(mean=1, cov=1e-1)
model$setInitialState(initial.state=initial.state)

## ----fig.height=5,fig.width=9, out.width="100%", fig.align='center'-----------
library(ggplot2)

# Set simulation settings
set.seed(11)
true.pars <- c(theta=10, mu=1, sigma_x=1, sigma_y=0.05)
dt.sim <- 1e-3
t.end <- 5
t.sim <- seq(0, t.end, by=dt.sim)
df.sim <- data.frame(t=t.sim, u=1, y=NA)

# Perform simulation
sim <- model$simulate(data=df.sim, 
                      pars=true.pars, 
                      n.sims=1, 
                      silent=T, 
                      initial.state=initial.state)
x <- sim$states$x$i0$x1

# Extract observations from simulation and add noise
iobs <- seq(1,length(t.sim), by=10)
t.obs <- t.sim[iobs]
y = x[iobs] + true.pars["sigma_y"] * rnorm(length(iobs))

# Create data-frame
data <- data.frame(
  t = t.obs,
  u = 1,
  y = y
)

# Plot the simulation and observed data
ggplot() +
  geom_line(aes(x=t.sim,y=x,color="Simulation")) +
  geom_point(aes(x=t.obs,y=y,fill="Observations")) +
  ctsmTMB:::getggplot2theme() + labs(x="t", y="x",color="",fill="")

## -----------------------------------------------------------------------------
fit <- model$estimate(data)

## ----eval=FALSE,echo=FALSE----------------------------------------------------
# # fit2 <- model$estimate(data, method="lkf")
# # fit3 <- model$estimate(data, method="laplace")

## -----------------------------------------------------------------------------
names(fit)

## -----------------------------------------------------------------------------
fit$convergence

## -----------------------------------------------------------------------------
fit$nll

## -----------------------------------------------------------------------------
fit$nll.gradient

## -----------------------------------------------------------------------------
fit$nll.hessian

## -----------------------------------------------------------------------------
print(fit)

## ----eval=FALSE, echo=FALSE---------------------------------------------------
# # print(fit2)
# # print(fit3)

## -----------------------------------------------------------------------------
fit$par.fixed

## -----------------------------------------------------------------------------
fit$sd.fixed

## -----------------------------------------------------------------------------
fit$cov.fixed

## -----------------------------------------------------------------------------
solve(fit$nll.hessian)

## ----fig.height=5,fig.width=9, out.width="100%", fig.align='center', eval=F,include=F----
# # Extract time, and prior/posterior mean and standard deviation estimates
# t         = fit$states$mean$posterior$t
# xprior    = fit$states$mean$prior$x
# xprior_sd = fit$states$sd$prior$x
# xpost     = fit$states$mean$posterior$x
# xpost_sd  = fit$states$sd$posterior$x
# 
# # Create vector c(xprior[1], xpost[1], xprior[2], xpost[2],...)
# t2 <- c(rbind(t,t))
# xprior_post <- c(rbind(xprior, xpost))
# xprior_post_sd <- c(rbind(xprior_sd,xpost_sd))
# 
# # Plot
# ggplot() +
#   geom_line(aes(x=t2,y=xprior_post,color="State Estimates (Prior-Posterior)"),lwd=1) +
#   geom_point(aes(x=data$t,data$y, color="Observations")) +
#   labs(x = "Time", y = "", color="") +
#   ctsmTMB:::getggplot2theme()

## ----eval=FALSE, include=FALSE------------------------------------------------
#   # geom_line(aes(x=t,y=xpost,color="State Estimates (Posterior)"),lwd=1) +
#   # geom_line(aes(x=t,y=xprior,color="State Estimates (Prior)"),lwd=1) +
#   # geom_ribbon(aes(x=t,ymin=xpost-2*xpost_sd,ymax=xpost+2*xpost_sd),fill="grey",alpha=0.5) +

## ----fig.height=5,fig.width=7, out.width="100%", fig.align='center'-----------
plot(fit, type="states",state.type="prior",against="y")

## ----fig.height=5,fig.width=7, out.width="100%", fig.align='center'-----------
plot(fit, type="states",state.type="posterior",against="y")

## -----------------------------------------------------------------------------
names(fit$residuals)

## ----fig.height=9,fig.width=9, out.width="100%", fig.align='center'-----------
plot(fit)

## ----fig.height=5,fig.width=9, out.width="100%", fig.align='center'-----------
a <- fit$par.fixed["theta"] - 3*fit$sd.fixed["theta"]
b <- fit$par.fixed["theta"] + 3*fit$sd.fixed["theta"]
prof <- profile(fit, list("theta"=seq(a,b,length.out=50)), silent=TRUE)
plot(prof)

## ----fig.height=5,fig.width=9, out.width="100%", fig.align='center'-----------
# a <- fit$par.fixed["mu"] - 8*fit$sd.fixed["mu"]
# b <- fit$par.fixed["mu"] + 8*fit$sd.fixed["mu"]
# prof <- profile(fit, list("mu"=seq(a,b,length.out=50)), silent=TRUE)
# plot(prof)

## ----fig.height=5,fig.width=9, out.width="100%", fig.align='center'-----------
# a <- fit$par.fixed["theta"] - 5*fit$sd.fixed["theta"]
# b <- fit$par.fixed["theta"] + 5*fit$sd.fixed["theta"]
# c <- fit$par.fixed["mu"] - 10*fit$sd.fixed["mu"]
# d <- fit$par.fixed["mu"] + 10*fit$sd.fixed["mu"]
# prof <- profile(fit, list("theta"=seq(a,b,length.out=50),
#                           "mu"=seq(c,d,length.out=50)), silent=TRUE)
# plot(prof)

## -----------------------------------------------------------------------------
model$setAlgebraics(theta ~ exp(logtheta))

## -----------------------------------------------------------------------------
model$setParameter(logtheta = log(c(initial=5, lower=0, upper=20)))

