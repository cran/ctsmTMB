## ----include=FALSE------------------------------------------------------------
library(ctsmTMB)

## -----------------------------------------------------------------------------
# Create model object
obj = ctsmTMB$new()

# Add system equations
obj$addSystem(
  dx ~ theta * (mu-x) * dt + sigma_x*dw
)

## -----------------------------------------------------------------------------
obj$addObs(
  log(y) ~ x, obsnames = "log_y"
)

## ----error=TRUE---------------------------------------------------------------
try({
obj$setVariance(
  y ~ sigma_y^2
)
})

## -----------------------------------------------------------------------------
obj$setVariance(
  log_y ~ sigma_y^2
)

## ----eval=FALSE---------------------------------------------------------------
# obj$addObs(
#   log(y) ~ x,
#   y ~ x,
#   y^2+z^3 ~ x,
#   obsnames = c("log_y", NA, "y2_plus_z3")
# )

