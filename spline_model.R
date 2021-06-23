model{
  for (i in 1:n){
    y[i] ~ dbinom(p[i], N[i])
    logit(p[i]) <- mu[i]
    mu[i] = X[i, ] %*% beta[1:K]
  }
  
  beta[1:K] ~ dmnorm(beta.0[1:K,1] + beta.00, lambda*Q[1:K,1:K])
  
  gamma <- beta - beta.00
  
  lambda ~ dgamma(0.001, 0.001)
  
  beta.00 ~ dnorm(0, 1e-6)
  
  for (j in 1:m){
    y.rep[j] ~ dbinom(p.rep[j], 100)
    logit(p.rep[j]) <- mu.rep[j]
    mu.rep[j] = X.pred[j, ] %*% gamma[1:K] + beta.00
  }
}