library(gtools)
library(purrr)
library(ranger)
library(dplyr)

hard_threshold <- function(x, delta){
  ifelse(abs(x) < delta, 0, x)
}

# Create lag matrix
lagmat <- function(y, z, k){
  z <- as.numeric(z)
  
  mat <- matrix(NA, nrow = length(y), ncol = k + 2)
  
  mat[,1] <- y
  mat[,2] <- z
  for (ii in 1:k) {
    mat[,ii + 2] <- lag(y, n = ii)
  }
  df <- data.frame(mat[-(1:k),])
  colnames(df) <- c(c("y", "z"), sapply(1:k, function(i){paste("ym", toString(i), sep = '')}))
  
  df$y <- as.factor(df$y)
  for (ii in 1:k) {
    df[, ii + 2] <- as.factor(df[, ii + 2]) 
  } 
  
  return(df)
}

ar <- c(.15, .10, .10, .05)
ma <- c(.50, .40, .30, .20)

L <- 10000

eta <- 0.25
beta <- 1
rho <- 1
gamma <- 0
delta <- 0.5

set.seed(100)
z <- arima.sim(list(ar=ar, ma=ma), n=L)
xp <- arima.sim(list(ar=ar, ma=ma), n=L)
xpp <- arima.sim(list(ar=ar, ma=ma), n=L)

x <- sqrt(1-eta)*xpp + sqrt(eta)*z

p <- inv.logit(gamma * hard_threshold(x, delta) + beta * z + rho * xp)
y <- as.numeric(rbernoulli(n = L, p = p))

k <- 5
marginal_data <- lagmat(y, z, k)
marginal_model <- ranger(y ~ ., 
                         data = marginal_data, 
                         probability = TRUE,
                         min.node.size = 10,
                         importance = 'impurity')
phat <- marginal_model$predictions[,2]

plot(tail(p, length(p) - k), phat)
cor(tail(p, length(p) - k), phat)

marginal_z_model <- ranger(y ~ z, 
                           data = data.frame(y = as.factor(y), z = z), 
                           probability = TRUE)
phat_start <- marginal_z_model$predictions[1:k, 2]

phat_full <- c(phat_start, phat)

ysim <- as.numeric(rbernoulli(n = L, p = phat_full))

plot(p[200:400], type = "l")
lines(phat_full[200:400], col = "red")

