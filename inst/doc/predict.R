## ----include=FALSE------------------------------------------------------------
library(ctsmTMB)
library(ggplot2)

## ----eval=FALSE---------------------------------------------------------------
# model$predict(data,
#               pars = NULL,
#               use.cpp = FALSE,
#               method = "ekf",
#               ode.solver = "euler",
#               ode.timestep = diff(data$t),
#               k.ahead = Inf,
#               return.k.ahead = NULL,
#               return.covariance = TRUE,
#               initial.state = self$getInitialState(),
#               estimate.initial.state = private$estimate.initial,
#               silent = FALSE
# )

## ----eval=FALSE---------------------------------------------------------------
# k.ahead * (nrow(data) - k.ahead)

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

# Set simulation settings
set.seed(20)
true.pars <- c(theta=20, sigma_x=1, sigma_y=5e-2)
dt.sim <- 1e-3
t.sim <- seq(0, 1, by=dt.sim)
u.sim <- cumsum(rnorm(length(t.sim),sd=0.1))
df.sim <- data.frame(t=t.sim, y=NA, u=u.sim)

# Perform simulation
sim <- model$simulate(data=df.sim, 
                      pars=true.pars, 
                      n.sims=1,
                      silent=T)
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

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
ggplot() + 
  geom_point(aes(x=t.obs, y=y, color="y(t)")) +
  geom_line(aes(x=t.obs, y=u.obs, color="u(t)")) +
  ctsmTMB:::getggplot2theme() + 
  theme(legend.text = ggplot2::element_text(size=15)) +
  labs(color="",x="Time",y="")

## ----message=FALSE------------------------------------------------------------
pred = model$predict(.data, k.ahead=nrow(.data)-1, pars=c(1, 1, 0.05))
pred1 = model$predict(.data, k.ahead=nrow(.data)-1, pars=c(10, 1, 0.05))
pred2 = model$predict(.data, k.ahead=nrow(.data)-1, pars=c(50, 1, 0.05))
pred3 = model$predict(.data, k.ahead=nrow(.data)-1, pars=c(100, 1, 0.05))

## -----------------------------------------------------------------------------
head(pred$states) 

## -----------------------------------------------------------------------------
head(pred$observations)

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
t <- pred$states$t.j
latex.str <- lapply(
  c(sprintf("theta[%s]", c(1,10,50,100)),"y[t[k]]"),
  str2expression
)
ggplot() +
  geom_line(aes(x=t, y=pred$states$x, color="label1")) +
  geom_line(aes(x=t, y=pred1$states$x, color="label2")) +
  geom_line(aes(x=t, y=pred2$states$x, color="label3")) +
  geom_line(aes(x=t, y=pred3$states$x, color="label4")) +
  # geom_line(aes(x=t, y=pred4$states$x, color="label5")) +
  # geom_line(aes(x=t, y=pred5$states$x, color="label6")) +
  geom_point(aes(x=t.obs, y=y, color="y(t)")) +
  scale_color_discrete(labels=latex.str) +
  ctsmTMB:::getggplot2theme() + 
  theme(legend.text = ggplot2::element_text(size=15)) +
  labs(color="",x="Time",y="")

## -----------------------------------------------------------------------------
fit = model$estimate(.data)
print(fit)

## -----------------------------------------------------------------------------
pred.horizon <- 25
pred = model$predict(.data, k.ahead=pred.horizon)

## -----------------------------------------------------------------------------
pred.H = pred$states[pred$states$k.ahead==pred.horizon,]

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
ggplot() +
  geom_line(aes(x=pred.H$t.j, y=pred.H$x,color="25-Step Predictions")) +
  geom_ribbon(aes(x=pred.H$t.j,ymin=pred.H$x-2*sqrt(pred.H$var.x),ymax=pred.H$x+2*sqrt(pred.H$var.x)),fill="grey",alpha=0.5) +
  geom_point(aes(x=t.obs,y,color="Observations")) +
  labs(color="",x="Time",y="") +
  ctsmTMB:::getggplot2theme()

## -----------------------------------------------------------------------------
rmse = c()
k.ahead = 1:pred.horizon
for(i in k.ahead){
  xy = data.frame(
    x = pred$states[pred$states$k.ahead==i,"x"],
    y = pred$observations[pred$observations$k.ahead==i,"y.data"]
  )
  rmse[i] = sqrt(mean((xy[["x"]] - xy[["y"]])^2))
}

## ----echo=FALSE, fig.height=5,fig.width=9,out.width="100%", fig.align='center'----
ggplot() +
  geom_line(aes(k.ahead, rmse), color="steelblue") + 
  geom_point(aes(k.ahead, rmse), color="red") +
  labs(
    title = "Root-Mean Square Errors for Different Prediction Horizons",
    x = "Prediction Steps",
    y = "Root-Mean-Square Errors"
  ) +
  ctsmTMB:::getggplot2theme()

