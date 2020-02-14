gradient_descent <- function(arm,tol,point,a_k=1,c_1=1e-4,rho=0.5){
  grad_f = calc_grad(arm$angle,arm$len,point)
  norm_grad_f = mag(grad_f)
  while (norm_grad_f>tol){
    p_k = -grad_f/norm_grad_f
    a_k = backtracking_line_search(arm,point,grad_f,p_k,a_k=a_k,c_1=c_1,rho=rho)
    arm$angle = arm$angle + a_k*p_k
    grad_f = calc_grad(arm$angle,arm$len,point)
    norm_grad_f = mag(grad_f)
  }
  return(arm)
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

backtracking_line_search <- function(arm,point,grad_f,p_k,a_k=1,c_1=1e-4,rho=0.5){
  f_k = f(arm$angle,arm$len,point)
  new_f_k = f(arm$angle + a_k*p_k,arm$len,point)
  while (new_f_k > (f_k + c_1*a_k*as.numeric(grad_f%*%p_k))){
    a_k = rho*a_k
    new_f_k = f(arm$angle + a_k*p_k,arm$len,point)
  }
  return(a_k)
}

