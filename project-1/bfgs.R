
BFGS <- function(arm,tol,point,beta = 0.5){
  # inital point x_0, convergence tolerance  epsilon > 0 , inverse hessian approximation H_0
  H_k = beta*diag(nrow(arm))
  k = 1
  grad_f = calc_grad(arm$angle,arm$len,point)
  arm_point = point_arm(arm)
  first = T
  while (mag(arm_point[nrow(arm_point),] - point) > tol){
    p_k = - H_k%*%grad_f
    # line search for a_k
    a_k = 0.1
    s_k = a_k*p_k
    new_arm = arm$angle + s_k
    new_grad_f = calc_grad(new_arm,arm$len,point)
    y_k = new_grad_f - grad_f
    if (!first){
      rho_k = 1/as.numeric(t(y_k)%*%s_k)
      H_k = (diag(length(y_k))-rho_k*s_k%*%t(y_k))%*%H_k%*%(diag(length(y_k))-rho_k*y_k%*%t(s_k)) + rho_k*s_k%*%t(s_k)
      k = k + 1
      arm$angle = new_arm
      arm_point = point_arm(arm)
      browser()
      grad_f = new_grad_f
    }else{
      first = F
      H_k = as.numeric((t(y_k)%*%s_k)/(t(y_k)%*%y_k))*diag(length(y_k))
    }
  }
  return(arm)
}
#line_search <- function(angle,)


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
