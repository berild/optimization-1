
fletcher_reeves <- function(arm,tol,point){
  start_time = Sys.time()
  f_k = f(arm$angle,arm$len,point)
  grad_f = calc_grad(arm$angle,arm$len,point)
  p_k = -grad_f
  k = 1
  while (mag(grad_f) > tol){
    a_k = wolfe_line_search(arm,point,p_k)
    if (a_k == 322){
      break
    }
    arm$angle = arm$angle + a_k*p_k
    new_grad_f = calc_grad(arm$angle,arm$len,point)
    beta = as.numeric((t(new_grad_f)%*%new_grad_f)/(t(grad_f)%*%grad_f))
    p_k = -new_grad_f + beta*p_k
    grad_f = new_grad_f 
    k = k +1
  }
  end_time = Sys.time()
  return(list(arm = arm, 
              k = k, 
              runtime = start_time - end_time,
              grad = mag(grad_f),
              f = f(arm$angle,arm$len,point)))
}


calc_grad <- function(angle,len,point){
  x = rep(0,length(angle))
  y = rep(0,length(angle))
  grad_f = rep(0,length(angle))
  for (k in seq(length(angle))){
    for (i in seq(k,length(angle))){
      x[k] = cos(sum(angle[1:i]))*len[i] + x[k]
      y[k] = sin(sum(angle[1:i]))*len[i] + y[k]
    }
    grad_f[k] = -y[k]*(x[1]-point[1]) + x[k]*(y[1]-point[2])
  }
  return(grad_f)
}


f <- function(angle,len,point){
  arm_point = c(0,0)
  for (i in seq(length(angle))){
    arm_point = c(len[i]*cos(sum(angle[1:i])),len[i]*sin(sum(angle[1:i]))) + arm_point
  }
  return(1/2*sum((arm_point-point)^2))
}

