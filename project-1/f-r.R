
fletcher_reeves <- function(arm,tol,point){
  start_time = Sys.time()
  f_k = f(arm$angle,arm$len,point)
  grad_f = calc_grad(arm$angle,arm$len,point)
  p_k = -grad_f
  k = 1
  while (mag(grad_f) > tol){
    # line search, strong wolfe conditions
    a_k = wolfe_line_search(arm,point,p_k,c_2 = 0.45)
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
              runtime = as.numeric(end_time - start_time,units="secs"),
              grad = mag(grad_f),
              f = f(arm$angle,arm$len,point)))
}
