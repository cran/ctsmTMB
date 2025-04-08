## ----include=FALSE------------------------------------------------------------
library(ctsmTMB)
library(ggplot2)

## ----eval=FALSE---------------------------------------------------------------
# model$estimate(
#   data,
#   method = "ekf",
#   ode.solver = "euler",
#   ode.timestep = diff(data$t),
#   loss = "quadratic",
#   loss_c = NULL,
#   control = list(trace=1,iter.max=1e5,eval.max=1e5),
#   use.hessian = FALSE,
#   laplace.residuals = FALSE,
#   unconstrained.optim = FALSE,
#   estimate.initial.state = FALSE,
#   silent = FALSE
# )

## ----echo=FALSE, fig.height=5,fig.width=7,out.width="100%", fig.align='center', fig.cap="Loss Functions", warning=FALSE,message=FALSE----
c <- 5
f1 <- function(r) r^2
f2 <- function(r) ifelse(r <= c, r^2, c*(2*r-c))
f3 <- function(r) ifelse(r <= c, r^2, c^2)
sigmoid <- function(r_sqr) 1/(1+exp(-5*(sqrt(r_sqr)-c)))
huber.loss <- function(r_sqr) {
  s <- sigmoid(r_sqr)
  r_sqr * (1-s) + c * (2*sqrt(r_sqr)-c)*s
}
tukey.loss <- function(r_sqr) {
  s <- sigmoid(r_sqr)
  r_sqr * (1-s) + c^2*s
}

r <- seq(0,25,by=1e-2)
ggplot() +
  geom_line(aes(x=r,y=f3(r),col="Tukey"),linewidth=2) +
  geom_line(aes(x=r,y=f2(r),col="Huber"),linewidth=2) +
  geom_line(aes(x=r,y=f1(r),col="Quadratic",)) +
  geom_line(aes(x=r,y=huber.loss(r^2),col="Huber Smooth"),linewidth=0.5) +
  geom_line(aes(x=r,y=tukey.loss(r^2),col="Tukey Smooth"),linewidth=0.5) +
  geom_line(aes(x=c(c,c),y=c(0,9*c^2)),linetype="dashed",col="black") +
  geom_text(aes(x=c,y=c^2), hjust=2,vjust=-1,label="c",size=5) +
  ctsmTMB:::getggplot2theme() +
  scale_x_continuous(limits=c(0,3*c)) + scale_y_continuous(limits=c(0,9*c^2)) +
  scale_color_discrete(breaks=c("Quadratic","Huber","Huber Smooth","Tukey","Tukey Smooth")) +
  labs(x="r",y="Likelihood",color="") 

## -----------------------------------------------------------------------------
m <- 1
qchisq(0.95,df=m)

