library(parallel)
source("./project-1/general_functions.R")
source("./project-1/gd.R")
source("./project-1/bfgs.R")
source("./project-1/f-r.R")

find_inner_radius <-function(len){
  if (max(len) > sum(len[-which.max(len)])){
    return(max(len) - sum(len[-which.max(len)]))
  }else{
    return(0)
  }
}

generate_problem <- function(n){
  len = sample(seq(15),n,replace = T)
  angle = runif(n,min = 0, max = 2*pi)
  max_radius = sum(len)
  min_radius = find_inner_radius(len)
  tmp_r = runif(1,min=0, max = max_radius*1.2)
  tmp_angle = runif(1,min = 0, max = 2*pi)
  point = c (tmp_r*cos(tmp_angle),tmp_r*sin(tmp_angle))
  return(list(point = point, arm = data.frame(angle = angle, len = len)))
}


test_conv <- function(grad_f, eval_f,tol_1 = 1e-3, tol_2 = 1e-3){
  if ((grad_f < tol_1) & (eval_f < tol_2)){
    return(1)
  }else if (grad_f < tol_1){
    return(2)
  }else{
    return(0)
  }
}

run <- function(n_iter){
  n_arms = seq(3,10)
  res_bfgs = res_fr = res_gd = data.frame(
    steps = rep(NA,n_iter*length(n_arms)), 
    runtime = rep(NA,n_iter*length(n_arms)),
    conv = rep(NA,n_iter*length(n_arms)),
    n = rep(NA,n_iter*length(n_arms)),
    grad_f = rep(NA,n_iter*length(n_arms)),
    f = rep(NA,n_iter*length(n_arms))
  )
  k = 1
  pb <- txtProgressBar(min = 1, max = n_iter*length(n_arms), style = 3)
  for (n in n_arms){
    optim_list <- mclapply(seq(n_iter),function(x){
      prob = generate_problem(n)
      tmp_bfgs = BFGS(arm = prob$arm,tol = 0.001, point = prob$point)
      res_bfgs = c(tmp_bfgs$k,tmp_bfgs$runtime,test_conv(tmp_bfgs$grad,tmp_bfgs$f),n,tmp_bfgs$grad,tmp_bfgs$f)
      tmp_fr = fletcher_reeves(arm = prob$arm,tol = 0.001, point = prob$point)
      res_fr = c(tmp_fr$k,tmp_fr$runtime,test_conv(tmp_fr$grad,tmp_fr$f),n,tmp_fr$grad,tmp_fr$f)
      tmp_gd = gradient_descent(arm = prob$arm,tol = 0.001, point = prob$point)
      res_gd = c(tmp_gd$k,tmp_gd$runtime,test_conv(tmp_gd$grad,tmp_gd$f),n,tmp_gd$grad,tmp_gd$f)
      return(list(res_bfgs = res_bfgs, res_gd = res_gd, res_fr = res_fr))
    },mc.cores = detectCores())
    for (ele in optim_list){
      setTxtProgressBar(pb, k)
      res_bfgs[k,] = ele$res_bfgs
      res_fr[k,] = ele$res_fr
      res_gd[k,] = ele$res_gd
      k = k+1
    }
  }
  return(list(res_bfgs = res_bfgs, res_gd = res_gd, res_fr = res_fr))
}


res = run(10000)
save(res, file = "./project-1/optim-runs.Rdata")
get_stats <- function(res,tot_it){
  stats = data.frame(
    n = as.factor(c(sapply(seq(3,10),function(x){rep(x,3)}))),
    method = rep(c("GD","F-R","BFGS"),length(seq(3,10))),
    conv = rep(NA,3*length(seq(3,10))),
    it_q25 = rep(NA,3*length(seq(3,10))),
    it_q50 = rep(NA,3*length(seq(3,10))),
    it_q75 = rep(NA,3*length(seq(3,10))),
    run_q25 = rep(NA,3*length(seq(3,10))),
    run_q50 = rep(NA,3*length(seq(3,10))),
    run_q75 = rep(NA,3*length(seq(3,10)))
  )
  i = 1
  for (n in seq(3,10)){
    stats$conv[i] = sum(res$res_gd$conv[res$res_gd$n==n]!=0)/tot_it*100
    stats$conv[i+1] = sum(res$res_fr$conv[res$res_fr$n==n]!=0)/tot_it*100
    stats$conv[i+2] = sum(res$res_bfgs$conv[res$res_bfgs$n==n]!=0)/tot_it*100
    stats$it_q25[i] = quantile(res$res_gd$steps[res$res_gd$n==n],0.25)
    stats$it_q25[i+1] = quantile(res$res_fr$steps[res$res_fr$n==n],0.25)
    stats$it_q25[i+2] = quantile(res$res_bfgs$steps[res$res_bfgs$n==n],0.25)
    stats$it_q50[i] = quantile(res$res_gd$steps[res$res_gd$n==n],0.50)
    stats$it_q50[i+1] = quantile(res$res_fr$steps[res$res_fr$n==n],0.50)
    stats$it_q50[i+2] = quantile(res$res_bfgs$steps[res$res_bfgs$n==n],0.50)
    stats$it_q75[i] = quantile(res$res_gd$steps[res$res_gd$n==n],0.75)
    stats$it_q75[i+1] = quantile(res$res_fr$steps[res$res_fr$n==n],0.75)
    stats$it_q75[i+2] = quantile(res$res_bfgs$steps[res$res_bfgs$n==n],0.75)
    stats$run_q25[i] = quantile(res$res_gd$runtime[res$res_gd$n==n],0.25)
    stats$run_q25[i+1] = quantile(res$res_fr$runtime[res$res_fr$n==n],0.25)
    stats$run_q25[i+2] = quantile(res$res_bfgs$runtime[res$res_bfgs$n==n],0.25)
    stats$run_q50[i] = quantile(res$res_gd$runtime[res$res_gd$n==n],0.50)
    stats$run_q50[i+1] = quantile(res$res_fr$runtime[res$res_fr$n==n],0.50)
    stats$run_q50[i+2] = quantile(res$res_bfgs$runtime[res$res_bfgs$n==n],0.50)
    stats$run_q75[i] = quantile(res$res_gd$runtime[res$res_gd$n==n],0.75)
    stats$run_q75[i+1] = quantile(res$res_fr$runtime[res$res_fr$n==n],0.75)
    stats$run_q75[i+2] = quantile(res$res_bfgs$runtime[res$res_bfgs$n==n],0.75)
    i = i + 3
  }
  return(stats)
}
stats = get_stats(res,10000)
save(stats, file = "./project-1/optim-stats.Rdata")

