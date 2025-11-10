model{
### PRIORS 
## occupancy intercept estimate of community (community mean)
beta0.mean ~ dnorm(0, 0.05)
beta0.tau ~ dgamma(0.1, 0.1)
beta0.sigma <- sqrt(1 / beta0.tau)

## detection intercept estimate of community (community mean)
alpha0.mean ~ dnorm(0, 0.05)
alpha0.tau ~ dgamma(0.1, 0.1)
alpha0.sigma <- sqrt(1 / alpha0.tau)

## Data augmentation parameter
# < empty > 

## Continuous site covariates on detection - Fixed effects
# < empty > 

## Continuous site covariates on detection - Independent effects
# < empty > 

## Continuous site covariates on detection - with random effects
# Covariate: Duration|Species

alpha.ranef.cont.Duration.mean ~ dnorm(0, 0.05)
alpha.ranef.cont.Duration.tau ~ dgamma(0.1, 0.1)
alpha.ranef.cont.Duration.sigma <- sqrt(1 / alpha.ranef.cont.Duration.tau)


# Covariate: JulianDate_Deployment|Species

alpha.ranef.cont.JulianDate_Deployment.mean ~ dnorm(0, 0.05)
alpha.ranef.cont.JulianDate_Deployment.tau ~ dgamma(0.1, 0.1)
alpha.ranef.cont.JulianDate_Deployment.sigma <- sqrt(1 / alpha.ranef.cont.JulianDate_Deployment.tau)


## Categorical site covariates on detection - Fixed effect
# < empty > 

## Categorical site covariates on detection - with random effects
# < empty > 

## Continuous observation-level covariates on detection - Fixed effects
# Covariate: effort
alpha.obs.fixed.cont.effort ~ dnorm(0, 0.05)

## Continuous observation-level covariates on detection - with random effects
# < empty > 

## Categorical observation-level covariates on detection - Fixed effect
# < empty > 

## Categorical observation-level covariates on detection - with random effects
# < empty > 

## Continuous site covariates on Occupancy - Fixed effects
# < empty > 

## Continuous site covariates on Occupancy - Independent effects
# < empty > 

## Continuous site covariates on occupancy - with random effects
# Covariate: d2hs|Species

beta.ranef.cont.d2hs.mean ~ dnorm(0, 0.05)
beta.ranef.cont.d2hs.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.d2hs.sigma <- sqrt(1 / beta.ranef.cont.d2hs.tau)


# Covariate: elevation|Species

beta.ranef.cont.elevation.mean ~ dnorm(0, 0.05)
beta.ranef.cont.elevation.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.elevation.sigma <- sqrt(1 / beta.ranef.cont.elevation.tau)


# Covariate: slope|Species

beta.ranef.cont.slope.mean ~ dnorm(0, 0.05)
beta.ranef.cont.slope.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.slope.sigma <- sqrt(1 / beta.ranef.cont.slope.tau)


# Covariate: canopy|Species

beta.ranef.cont.canopy.mean ~ dnorm(0, 0.05)
beta.ranef.cont.canopy.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.canopy.sigma <- sqrt(1 / beta.ranef.cont.canopy.tau)


# Covariate: evi|Species

beta.ranef.cont.evi.mean ~ dnorm(0, 0.05)
beta.ranef.cont.evi.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.evi.sigma <- sqrt(1 / beta.ranef.cont.evi.tau)


# Covariate: precipitation|Species

beta.ranef.cont.precipitation.mean ~ dnorm(0, 0.05)
beta.ranef.cont.precipitation.tau ~ dgamma(0.1, 0.1)
beta.ranef.cont.precipitation.sigma <- sqrt(1 / beta.ranef.cont.precipitation.tau)


## Categorical site covariates on Occupancy - Fixed effects
# < empty > 

## Categorical site covariates on occupancy - with random effects
# < empty > 

# Species-station random effect on detection probability
# < empty > 

## Draws of random effects other than species


### MODEL LOOPS 

# species loop
for (i in 1:M){
##  Draw species-specific random effect parameters from community distributions
# intercepts:
beta0[i] ~ dnorm(beta0.mean, beta0.tau)
alpha0[i] ~ dnorm(alpha0.mean, alpha0.tau)


# continuous detection covariate with random effects: Duration|Species
alpha.ranef.cont.Duration[i] ~ dnorm(alpha.ranef.cont.Duration.mean, alpha.ranef.cont.Duration.tau)
# continuous detection covariate with random effects: JulianDate_Deployment|Species
alpha.ranef.cont.JulianDate_Deployment[i] ~ dnorm(alpha.ranef.cont.JulianDate_Deployment.mean, alpha.ranef.cont.JulianDate_Deployment.tau)


# categorical detection covariates: no random effect of species


# continuous observation-level detection covariates: no random effect of species


# categorical observation covariates: no random effect of species


# continuous occupancy covariate with random effects: d2hs|Species
beta.ranef.cont.d2hs[i] ~ dnorm(beta.ranef.cont.d2hs.mean, beta.ranef.cont.d2hs.tau)
# continuous occupancy covariate with random effects: elevation|Species
beta.ranef.cont.elevation[i] ~ dnorm(beta.ranef.cont.elevation.mean, beta.ranef.cont.elevation.tau)
# continuous occupancy covariate with random effects: slope|Species
beta.ranef.cont.slope[i] ~ dnorm(beta.ranef.cont.slope.mean, beta.ranef.cont.slope.tau)
# continuous occupancy covariate with random effects: canopy|Species
beta.ranef.cont.canopy[i] ~ dnorm(beta.ranef.cont.canopy.mean, beta.ranef.cont.canopy.tau)
# continuous occupancy covariate with random effects: evi|Species
beta.ranef.cont.evi[i] ~ dnorm(beta.ranef.cont.evi.mean, beta.ranef.cont.evi.tau)
# continuous occupancy covariate with random effects: precipitation|Species
beta.ranef.cont.precipitation[i] ~ dnorm(beta.ranef.cont.precipitation.mean, beta.ranef.cont.precipitation.tau)


# categorical occupancy covariates: no random effect of species

# station loop
for (j in 1:J){

# Occupancy probability formula

logit.psi[i,j] <- beta0[i] + beta.ranef.cont.d2hs[i] * d2hs[j] + beta.ranef.cont.elevation[i] * elevation[j] + beta.ranef.cont.slope[i] * slope[j] + beta.ranef.cont.canopy[i] * canopy[j] + beta.ranef.cont.evi[i] * evi[j] + beta.ranef.cont.precipitation[i] * precipitation[j]
psi[i,j] <- exp(logit.psi[i,j]) / (exp(logit.psi[i,j]) + 1)
z[i,j] ~ dbern(psi[i, j])

# No random effect of species and station on detection probability
# occasion loop
for (k in 1:maxocc){
# Detection probability formula
logit.p[i,j,k] <- alpha0[i] + alpha.obs.fixed.cont.effort * effort[j, k] + alpha.ranef.cont.Duration[i] * Duration[j] + alpha.ranef.cont.JulianDate_Deployment[i] * JulianDate_Deployment[j]

# convert p to real scale
p[i,j,k] <- exp(logit.p[i,j,k]) / (1+exp(logit.p[i,j,k]))

# Ensure occasions without effort have p = 0


p.eff[i,j,k] <- z[i,j] * p[i,j,k] * effort_binary[j,k]
y[i,j,k] ~ dbern(p.eff[i,j,k])

### generate new data from model under consideration
new.y[i,j,k] ~ dbern(p.eff[i,j,k])
}   # close occasion loop

### calculate Freeman-Tukey residuals for real and new data
res[i,j] <- (sqrt(sum(y[i,j, 1:maxocc])) - sqrt(sum(p.eff[i,j, 1:maxocc])))^2
new.res[i,j] <- (sqrt(sum(new.y[i,j, 1:maxocc])) - sqrt(sum(p.eff[i,j, 1:maxocc])))^2
}   # close station loop

### sum residuals over stations
R2[i] <- sum(res[i, 1:J])
new.R2[i] <- sum(new.res[i, 1:J])

### species-level Bayesian p-value
Bpvalue_species[i] <- R2[i] > new.R2[i]


### Number of occupied stations for each species
NStationsOccupied[i] <- sum(z[i, 1:J])

### species is part of community?
speciesInCommunity[i] <- 1 - equals(NStationsOccupied[i],0)

### Fraction of stations occupied
fractionStationsOccupied[i] <- NStationsOccupied[i] /J

### Does species occur at all (at any station)
occ[i] <- 1 - equals(fractionStationsOccupied[i], 0)
##############
######### added code for mean psi and p per species
############

mu.psi[i] <- mean(psi[i, 1:J])
mu.p[i] <- mean(p.eff[i,1:J,])

#######################
#######################


}    # close species loop

###sum residuals over observed species
R3 <- sum(R2[1:M])
new.R3 <- sum(new.R2[1:M])
Bpvalue <- mean(Bpvalue_species[1:M]) ########## EDITED to give corrected BPV from average across each estimate

### total number of species
Nspecies <- sum(speciesInCommunity[1:M])

### Species richness at every location
for (j in 1:J){
Nspecies_station[j] <- sum(z[1:M,j])
}

}