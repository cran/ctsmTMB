## ----include=FALSE------------------------------------------------------------
library(ctsmTMB)

## -----------------------------------------------------------------------------
# Simulate data using Euler Maruyama
set.seed(10)
theta=10; mu=1; sigma_x=1; sigma_y=1e-1
# 
dt.sim = 1e-3
t.sim = seq(0,1,by=dt.sim)
dw = rnorm(length(t.sim)-1,sd=sqrt(dt.sim))
x = 3
for(i in 1:(length(t.sim)-1)) {
  x[i+1] = x[i] + theta*(mu-x[i])*dt.sim + sigma_x*dw[i]
}

# Extract observations and add noise
dt.obs = 1e-2
t.obs = seq(0,1,by=dt.obs)
y = x[t.sim %in% t.obs] + sigma_y * rnorm(length(t.obs))

# Create data
.data = data.frame(
  t = t.obs,
  y = y
)

## -----------------------------------------------------------------------------
# Create model object
obj = ctsmTMB$new()

# Add system equations
obj$addSystem(
  dx ~ theta * (mu-x) * dt + sigma_x*dw
)

# Add observation equations
obj$addObs(
  y ~ x
)

# Set observation equation variances
obj$setVariance(
  y ~ sigma_y^2
)

# Specify algebraic relations
obj$setAlgebraics(
  theta ~ exp(logtheta),
  sigma_x ~ exp(logsigma_x),
  sigma_y ~ exp(logsigma_y)
)

# Specify parameter initial values and lower/upper bounds in estimation
obj$setParameter(
  logtheta   = log(c(initial = 5,    lower = 0,    upper = 20)),
  mu         = c(    initial = 0,    lower = -10,  upper = 10),
  logsigma_x = log(c(initial = 1e-1, lower = 1e-5, upper = 5)),
  logsigma_y = log(c(initial = 1e-1, lower = 1e-5, upper = 5))
)

# Set initial state mean and covariance
obj$setInitialState(list(x[1], 1e-1*diag(1)))

## -----------------------------------------------------------------------------
fit = obj$estimate(.data)

## ----eval=FALSE---------------------------------------------------------------
# nll = TMB::MakeADFun(...)
# opt = stats::nlminb(start=nll$par, objective=nll$fn, grad=nll$gr, hessian=nll$he)

## -----------------------------------------------------------------------------
nll = obj$likelihood(.data)

## -----------------------------------------------------------------------------
nll$par

## -----------------------------------------------------------------------------
nll$fn(nll$par)

## -----------------------------------------------------------------------------
nll$gr(nll$par)

## -----------------------------------------------------------------------------
nll$he(nll$par)

## -----------------------------------------------------------------------------
pars = obj$getParameters()
print(pars)

## -----------------------------------------------------------------------------
opt = stats::optim(par=nll$par, 
                   fn=nll$fn, 
                   gr=nll$gr, 
                   method="L-BFGS-B", 
                   lower=pars$lower, 
                   upper=pars$upper)

## -----------------------------------------------------------------------------
# Estimated parameters
data.frame(external=opt$par, internal=fit$par.fixed)

# Neg. Log-Likelihood
data.frame(external=opt$value, internal=fit$nll)

# Gradient components
data.frame(external=t(nll$gr(opt$par)), internal=t(nll$gr(fit$par.fixed)))

