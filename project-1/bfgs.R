
f <- function(angle,len,point){
  arm_point = c(0,0)
  for (i in seq(length(angle))){
    arm_point = c(len[i]*cos(sum(angle[1:i])),len[i]*sin(sum(angle[1:i]))) + arm_point
  }
  return(1/2*sum((arm_point-point)^2))
}

gradient_descent <- function(arm,point,tol,a_k=1,c_1=1e-4,rho=0.5){
  grad_f = calc_grad(arm$angle,arm$len,point)
  norm_grad_f = mag(grad_f)
  while (norm_grad_f>epsilon){
    p_k = -grad/norm_grad_f
    a_k = backtracking_line_search(arm,point,grad_f,p_k,a_k=a_k,c_1=c_1,rho=rho)
    arm$angle = arm$angle + a_k*p_k
    grad_f = calc_grad(arm$angle,arm$len,point)
    norm_grad_f = mag(grad_f)
  }
  return(arm)
}

BFGS <- function(arm,tol,point,beta = 0.5){
  # inital point x_0, convergence tolerance  epsilon > 0 , inverse hessian approximation H_0
  H_k = beta*diag(nrow(arm))
  k = 1
  grad_f = calc_grad(arm$angle,arm$len,point)
  first = T
  a_k = 1
  while (mag(grad_f) > tol){
    p_k = - H_k%*%grad_f
    # line search for a_k
    a_k = wolfe_line_search(arm,point,p_k)
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
      first = F
      H_k = as.numeric((t(y_k)%*%s_k)/(t(y_k)%*%y_k))*diag(length(y_k))
    }
  }
  return(arm)
}
#line_search <- function(angle,)

backtracking_line_search <- function(arm,point,grad_f,p_k,a_k=1,c_1=1e-4,rho=0.5){
  f_k = f(arm$angle,arm$len,point)
  new_f_k = f(arm$angle + a_k*p_k,arm$len,point)
  while (new_f_k > (f_k + c_1*a_k*as.numeric(grad_f%*%p_k))){
    a_k = rho*a_k
    new_f_k = f(arm$angle + a_k*p_k,arm$len,point)
  }
  return(a_k)
}

wolfe_line_search <- function(arm,point,p_k,c_1=1e-4,c_2=0.9,steps = 20){
  a_0 = a_i = 0
  a_m = 1
  phi_0 = phi_i = f(arm$angle,arm$len,point)
  d_phi_0 = phi_deriv(arm,point,p_k,a_0)
  for (i in seq(steps)){
    phi_m = f(arm$angle+a_m*p_k,arm$len,point)
    if ((phi_m > (phi_0 + c_1*a_m*d_phi_0)) | 
        ((phi_m >= phi_i) & (i > 0)) ){
      return(zoom(arm,point,p_k,phi_0,d_phi_0,a_i,a_m,c_1=c_1,c_2=c_2))
    }
    d_phi_m = phi_deriv(arm,point,p_k,a_m)
    if (abs(d_phi_m)<=(-c_2*d_phi_0)){
      return(a_m)
    }
    if (d_phi_m>=0){
      return(zoom(arm,point,p_k,phi_0,d_phi_0,a_m,a_i,c_1=c_1,c_2=c_2))
    }
    a_i = a_m
    a_m = 2*a_i
    phi_i = phi_m
  }
}

zoom <- function(arm,point,p_k,phi_0,d_phi_0,a_l,a_h,c_1 = 1e-4,c_2=0.9){
  while (T){
    a_j = (a_l + a_h)/2
    phi_j = f(arm$angle + a_j*p_k,arm$len,point)
    phi_l = f(arm$angle + a_l*p_k,arm$len,point)
    if ((phi_j > (phi_0 + c_1*a_j*d_phi_0)) | (phi_j>=phi_l)){
      a_h = a_j
    }else{
      d_phi_j = phi_deriv(arm,point,p_k,a_j)
      
      if (abs(d_phi_j)<=(-c_2*d_phi_0)){
        return(a_j)
      }
      if ((d_phi_j*(a_h-a_l))>=0){
        a_h = a_l
      }
      a_l = a_j
    }
  }
}

phi_deriv <- function(arm,point,p_k,a_k){
  return(as.numeric(calc_grad(arm$angle + a_k*p_k, arm$len, point)%*%p_k))
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
