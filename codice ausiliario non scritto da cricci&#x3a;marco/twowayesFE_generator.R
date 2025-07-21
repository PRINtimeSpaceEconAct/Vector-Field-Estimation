nu = 100
T = 4 

## ID
idi = seq(1:nu)%x%rep(1,T) ## INDEX BASED ON INDIVIDUALS

idt = rep(1,nu)%x%seq(1:T) ## INDEX BASED ON TIME OCCASIONS



## INDIVIDUAL EFFECTS
alpha = rnorm(nu)

Alpha = alpha%x%rep(1,T)

gamma = rnorm(nu)

GGamma = rep(1,T)%x%gamma


## COVARIATES
X = Alpha + GGamma + rnorm(nu*T)

## DEP. VAR.
y = Alpha + GGamma + b1*X + rnorm(nu*T)
