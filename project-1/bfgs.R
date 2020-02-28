
BFGS <- function(arm,tol,point,beta = 0.5){
  # beta is a scalar for the initial hessian
  start_time = Sys.time()
  # initial hessian
  H_k = beta*diag(nrow(arm))
  k = 1
  grad_f = calc_grad(arm$angle,arm$len,point)
  first = T
  a_k = 1
  # magnitue of gradient as stopping criterion
  while (mag(grad_f) > tol){
    p_k = - H_k%*%grad_f
    # line search for a_k
    a_k = wolfe_line_search(arm,point,p_k)
    if (a_k == 322){
      break
    }
    s_k = a_k*p_k
    new_arm = arm$angle + s_k
    new_grad_f = calc_grad(new_arm,arm$len,point)
    y_k = new_grad_f - grad_f
    if (!first){
      rho_k = 1/as.numeric(t(y_k)%*%s_k)
      H_k = (diag(length(y_k))-rho_k*s_k%*%t(y_k))%*%H_k%*%(diag(length(y_k))-rho_k*y_k%*%t(s_k)) + rho_k*s_k%*%t(s_k)
      k = k + 1
      arm$angle = new_arm
      grad_f = new_grad_f
    }else{
      # only updating the hessian in the first step
      first = F
      H_k = as.numeric((t(y_k)%*%s_k)/(t(y_k)%*%y_k))*diag(length(y_k))
    }
  }
  end_time = Sys.time()
  return(list(arm = arm, 
              k = k, 
              runtime = as.numeric(end_time - start_time,units="secs"),
              grad = mag(grad_f),
              f = f(arm$angle,arm$len,point)))
}



