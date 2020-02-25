gradient_descent <- function(arm,tol,point,a_k=1,c_1=1e-4,rho=0.5){
  start_time = Sys.time()
  grad_f = calc_grad(arm$angle,arm$len,point)
  norm_grad_f = mag(grad_f)
  k = 1
  while (norm_grad_f>tol){
    p_k = -grad_f/norm_grad_f
    a_k = backtracking_line_search(arm,point,grad_f,p_k,a_k=a_k,c_1=c_1,rho=rho)
    arm$angle = arm$angle + a_k*p_k
    grad_f = calc_grad(arm$angle,arm$len,point)
    norm_grad_f = mag(grad_f)
    k = k+1
  }
  end_time = Sys.time()
  return(list(arm = arm, 
              k = k, 
              runtime = as.numeric(end_time - start_time,units="secs"),
              grad = mag(grad_f),
              f = f(arm$angle,arm$len,point)))
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

