# The objective function ----
oracle <- function(idx, x){
  # Objective function
  # Args:
    # idx=-1: simple quadratic well-conditioned 0.01*x_1^2 + 0.005*x_2^2, L/mu=2
    # idx=0: simple quadratic well-conditioned 0.5*x_1^2 + x_2^2, L/mu=2
    # idx=1: simple quadratic ill-conditioned 0.005*x_1^2 + x_2^2, L/mu=200
    # idx=2: quadratic
    # idx=3: log-sum-exp
    # idx=4: absolute value |x|
  
  if(idx == -1){
    obj_value <- 1*10^(-2)*x[1]^2 + 5*10^(-3)*x[2]^2
    grad <- c(1*10^(-2)*2*x[1],5*10^(-3)*2*x[2])
  }
  
  if(idx == 0){
    obj_value <- 0.5*x[1]^2 + x[2]^2
    grad <- c(x[1],2*x[2])
  }
  
  if(idx == 1){
    obj_value <- 5*10^(-3)*x[1]^2 + x[2]^2
    grad <- c(5*10^(-3)*2*x[1],2*x[2])
  }
  
  if(idx == 2){
    obj_value <- 0.5 * t(x)%*%A%*%x + t(b)%*%x
    grad <- A%*%x + b
  }
  
  if(idx == 3){
    rho <- 20
    temp_vec <- exp((t(A)%*%x - b)/rho)
    obj_value <- rho*log(sum(temp_vec))
    grad <- A%*%temp_vec/sum(temp_vec)
  }
  
  if(idx == 4){
    obj_value <- abs(x)
    grad <- sign(x)
  }
  
  return(list(obj=obj_value, grad = grad))
}

# Optimization algorithms----
GD <- function(idx, x0, s, MaxIter){
  # The gradient descent algorithm
  obj_track <- rep(0,MaxIter+1)
  gradSq_track <- rep(0,MaxIter+1)
  obj_track[1] <- oracle(idx,x0)$obj
  gradSq_track[1] <- sum(oracle(idx,x0)$grad^2)
  x_old <- x0
  
  for(k in 1:MaxIter){
    x_new <- x_old - s * oracle(idx,x_old)$grad
    obj_track[k+1] <- oracle(idx,x_new)$obj
    gradSq_track[k+1] <- sum(oracle(idx,x_new)$grad^2)
    x_old <- x_new
  }
  return(list(x = x_new, obj_track = obj_track, gradSq_track = gradSq_track))
}

extSC <- function(idx, x0, s, mu, c0=1, c1, c2, MaxIter){
  # the general algorithm for strongly convex obj (parametrized by c0, c1, and c2)
  q <- mu * s
  x_old <- x0 # previous iteration
  x <- x0 - 2*s/(1+sqrt(q))*oracle(idx,x0)$grad # current iteration
  
  if(idx <= 1){ #for simple quadratic (x is 2-dimension), also record the trace of x
    x_track <- rbind(x_old,x)
  }
  
  obj_track <- rep(0,MaxIter+2)
  gradSq_track <- rep(0,MaxIter+2)
  obj_track[1] <- oracle(idx,x0)$obj
  gradSq_track[1] <- sum(oracle(idx,x0)$grad^2)
  obj_track[2] <- oracle(idx,x)$obj
  gradSq_track[2] <- sum(oracle(idx,x)$grad^2)
  
  for(k in 1:MaxIter){
    x_new <- x - c0*s*oracle(idx,x)$grad + (1-c1*sqrt(q))*(x-x_old) - (c2*sqrt(c0)-c0/2)*s*(oracle(idx,x)$grad - oracle(idx,x_old)$grad)
    x_old <- x
    x <- x_new
    
    obj_track[k+2] <- oracle(idx,x)$obj
    gradSq_track[k+2] <- sum(oracle(idx,x)$grad^2)
    if((idx == 0) | (idx == 1)){
      x_track <- rbind(x_track,x)
    }
  }
  
  if((idx == 0) | (idx == 1)){
    return(list(x = x, obj_track = obj_track, gradSq_track = gradSq_track, x_track = x_track))
  }else{
    return(list(x = x, obj_track = obj_track, gradSq_track = gradSq_track))
  }
}

NAGSC <- function(idx, x0, s, mu, c0=1, MaxIter){
  # Nesterov's accelerated gradient method (for strongly convex obj)
  q <- mu * s
  x_old <- x0 # previous iteration
  x <- x0 - 2*s/(1+sqrt(q))*oracle(idx,x0)$grad # current iteration
  
  if(idx <= 1){ #for simple quadratic (x is 2-dimension), also record the trace of x
    x_track <- rbind(x_old,x)
  }
  
  obj_track <- rep(0,MaxIter+2)
  gradSq_track <- rep(0,MaxIter+2)
  obj_track[1] <- oracle(idx,x0)$obj
  gradSq_track[1] <- sum(oracle(idx,x0)$grad^2)
  obj_track[2] <- oracle(idx,x)$obj
  gradSq_track[2] <- sum(oracle(idx,x)$grad^2)
  
  for(k in 1:MaxIter){
    x_new <- x - c0*s*oracle(idx,x)$grad + (1-sqrt(q))/(1+sqrt(q))*(x-x_old) - (1-sqrt(q))/(1+sqrt(q))*s*(oracle(idx,x)$grad - oracle(idx,x_old)$grad)
    x_old <- x
    x <- x_new
    
    obj_track[k+2] <- oracle(idx,x)$obj
    gradSq_track[k+2] <- sum(oracle(idx,x)$grad^2)
    if((idx == 0) | (idx == 1)){
      x_track <- rbind(x_track,x)
    }
  }
  
  if((idx == 0) | (idx == 1)){
    return(list(x = x, obj_track = obj_track, gradSq_track = gradSq_track, x_track = x_track))
  }else{
    return(list(x = x, obj_track = obj_track, gradSq_track = gradSq_track))
  }
}

extC <- function(idx, x0, s, r, gamma = 1, beta, MaxIter){
  # the general algorithm for convex obj (parametrized by r, gamma, and beta)
  x_old <- x0 # previous iteration
  x <- x0 - gamma*s*oracle(idx,x0)$grad # current iteration
  
  obj_track <- rep(0,MaxIter+2)
  gradSq_track <- rep(0,MaxIter+2)
  obj_track[1] <- oracle(idx,x0)$obj
  gradSq_track[1] <- sum(oracle(idx,x0)$grad^2)
  obj_track[2] <- oracle(idx,x)$obj
  gradSq_track[2] <- sum(oracle(idx,x)$grad^2)
  
  for(k in 1:MaxIter){
    x_new <- x - gamma*s*oracle(idx,x)$grad + k/(k+r+1)*(x-x_old) - k/(k+r+1)*beta*s*(oracle(idx,x)$grad - oracle(idx,x_old)$grad)
    x_old <- x
    x <- x_new
    
    obj_track[k+2] <- oracle(idx,x)$obj
    gradSq_track[k+2] <- sum(oracle(idx,x)$grad^2)
  }
  
  return(list(x = x, obj_track = obj_track, gradSq_track = gradSq_track))
}


