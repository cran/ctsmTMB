## ----include=FALSE------------------------------------------------------------
library(ctsmTMB)
library(ggplot2)

## ----eval=FALSE---------------------------------------------------------------
# model$simulate(data,
#                pars = NULL,
#                use.cpp = FALSE,
#                method = "ekf",
#                ode.solver = "rk4",
#                ode.timestep = diff(data$t),
#                simulation.timestep = diff(data$t),
#                k.ahead = nrow(data)-1,
#                return.k.ahead = 0:k.ahead,
#                n.sims = 100,
#                initial.state = self$getInitialState(),
#                estimate.initial.state = private$estimate.initial,
#                silent = FALSE)

## -----------------------------------------------------------------------------
model = ctsmTMB$new()
model$addSystem(dx ~ theta * (t*u^2-cos(t*u) - x) * dt + sigma_x*dw)
model$addObs(y ~ x)
model$setVariance(y ~ sigma_y^2)
model$addInput(u)
model$setParameter(
  theta   = c(initial = 2, lower = 0,    upper = 100),
  sigma_x = c(initial = 0.2, lower = 1e-5, upper = 5),
  sigma_y = c(initial = 5e-2)
)
model$setInitialState(list(1, 1e-1*diag(1)))

## ----include=TRUE-------------------------------------------------------------
# Set simulation settings
set.seed(20)
true.pars <- c(theta=20, sigma_x=1, sigma_y=5e-2)
dt.sim <- 1e-3
t.sim <- seq(0, 1, by=dt.sim)
u.sim <- cumsum(rnorm(length(t.sim),sd=0.1))
df.sim <- data.frame(t=t.sim, y=NA, u=u.sim)

# Simulate data
sim <- model$simulate(data=df.sim, 
                      pars=true.pars, 
                      n.sims=1,
                      silent=T)

# Grab first simulation trajectory
x <- sim$states$x$i0$x1

# Extract observations from simulation and add noise
iobs <- seq(1,length(t.sim), by=10)
t.obs <- t.sim[iobs]
u.obs <- u.sim[iobs]
y = x[iobs] + true.pars["sigma_y"] * rnorm(length(iobs))

# Create data-frame
.data <- data.frame(
  t = t.obs,
  u = u.obs,
  y = y
)

## -----------------------------------------------------------------------------
sim <- model$simulate(data=.data, 
                      pars=c(20,1,0.05), 
                      n.sims=100,
                      silent=T)

## ----fig.height=5,fig.width=9,out.width="100%", fig.align='center'------------
# Get the first (and only in this case) k-step simulation data.frame
X <- sim$states$x$i0

# Grab all the simulations (the first five columns are indices, time, etc.)
Y <- X[,-c(1:5)]

# Grab prediction time column
t <- X[,"t.j"]

# Plot
matplot(t,Y,type="l", ylim=c(-4,4))

## -----------------------------------------------------------------------------
sim <- model$simulate(data=.data, 
                      pars=c(20,3,0.05), 
                      n.sims=100,
                      silent=T)

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
# Get the first (and only in this case) k-step simulation data.frame
X <- sim$states$x$i0

# Grab all the simulations (the first five columns are indices, time, etc.)
Y <- X[,-c(1:5)]

# Grab prediction time column
t <- X[,"t.j"]

# Plot
matplot(t,Y,type="l",ylim=c(-4,4))

## -----------------------------------------------------------------------------
sim <- model$simulate(data=.data, 
                      pars=c(50,1,0.05), 
                      n.sims=100,
                      silent=T)

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
# Get the first (and only in this case) k-step simulation data.frame
X <- sim$states$x$i0

# Grab all the simulations (the first five columns are indices, time, etc.)
Y <- X[,-c(1:5)]

# Grab prediction time column
t <- X[,"t.j"]

# Plot
matplot(t,Y,type="l", ylim=c(-4,4))

