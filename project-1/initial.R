# This script is used to visualize the solution using the different methods
# To run the simulation test look at run_stats.R
# At the bottom of the page the data from the simulation is loaded and plotted

library(ggplot2) # required for the plots
# we have implemented our trails simulation in parallel since
# it would take to long if not.
library(parallel) # required for parallel computations
source("./general_functions.R") # functions for all methods
source("./gd.R") # gradient descent functions
source("./bfgs.R") # bfgs functions
source("./f-r.R") # fletcher reeves functions

# intialize a arm
init_arm <- function(len,angle){
  arm_polar = data.frame(len = len, angle = angle)
  arm_cart = point_arm(arm_polar)
  return(list(polar = arm_polar,cart = arm_cart))
}

# forward kinematic problem 
# find point from segment and angle
point_arm <- function(arm){
  x = y = rep(0,nrow(arm)+1)
  for (i in seq(nrow(arm))){
    x[i+1] = cos(sum(arm$angle[1:i]))*arm$len[i] + x[i]
    y[i+1] = sin(sum(arm$angle[1:i]))*arm$len[i] + y[i]
  }
  return(data.frame(x = x, y = y))
}

# find inner radius if configuration space, is zero if no inner restrictions
find_inner_radius <-function(len){
  if (max(len) > sum(len[-which.max(len)])){
    return(max(len) - sum(len[-which.max(len)]))
  }else{
    return(0)
  }
}

# to create a visually satisfying end-effector
create_arm_tip <- function(arm_cart,arm_polar){
  tot_angle = sum(arm_polar$angle)
  tot_len = sum(arm_polar$len)
  angles = c(tot_angle-pi/2,tot_angle + pi/2)
  seg = c(tot_len/24,tot_len/24)
  n = nrow(arm_cart)
  res = data.frame(
    x1 = c(arm_cart$x[n] - seg[1]/2*cos(tot_angle), seg[2]*cos(angles[1]) + arm_cart$x[n] - seg[1]/2*cos(tot_angle),  seg[1]*cos(tot_angle)+ seg[2]*cos(angles[1]) + arm_cart$x[n] - seg[1]/2*cos(tot_angle)),
    x2 = c(arm_cart$x[n] - seg[1]/2*cos(tot_angle), seg[2]*cos(angles[2]) + arm_cart$x[n] - seg[1]/2*cos(tot_angle),  seg[1]*cos(tot_angle)+ seg[2]*cos(angles[2]) + arm_cart$x[n] - seg[1]/2*cos(tot_angle)),
    y1 = c(arm_cart$y[n] - seg[1]/2*sin(tot_angle), seg[2]*sin(angles[1]) + arm_cart$y[n] - seg[1]/2*sin(tot_angle),  seg[1]*sin(tot_angle)+ seg[2]*sin(angles[1]) + arm_cart$y[n] - seg[1]/2*sin(tot_angle)),
    y2 = c(arm_cart$y[n] - seg[1]/2*sin(tot_angle), seg[2]*sin(angles[2]) + arm_cart$y[n] - seg[1]/2*sin(tot_angle),  seg[1]*sin(tot_angle)+ seg[2]*sin(angles[2]) + arm_cart$y[n] - seg[1]/2*sin(tot_angle))
  )
  return(res)
}

# function plotting the whole arm for one method
plot_arm <- function(arm,point){
  require(ggplot2)
  max_radius = sum(arm$len)
  min_radius = find_inner_radius(arm$len)
  circ = circleFun(inner_radius = min_radius, radius = max_radius)
  arm_point = point_arm(arm)
  arm_text = data.frame(label = seq(nrow(arm_point)-1),x = arm_point$x[-nrow(arm_point)],y = arm_point$y[-nrow(arm_point)])
  arm_tip = create_arm_tip(arm_point,arm)
  arm_point$x[nrow(arm_point)] = arm_tip$x1[1]
  arm_point$y[nrow(arm_point)] = arm_tip$y1[1]
  arm_fig <- ggplot() + 
    geom_line(data = circ$out_up, aes(x=x,y=y)) + 
    geom_line(data = circ$in_up, aes(x=x,y=y)) + 
    geom_line(data = circ$out_low, aes(x=x,y=y)) + 
    geom_line(data = circ$in_low, aes(x=x,y=y)) + 
    geom_ribbon(data = circ$fill_up, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) +
    geom_ribbon(data = circ$fill_low, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) + 
    geom_path(data = arm_point, aes(x=x,y=y),color = "firebrick",size = 1) + 
    geom_path(data = arm_tip, aes(x=x1,y=y1), color = "firebrick",size = 1) + 
    geom_path(data = arm_tip, aes(x=x2,y=y2), color = "firebrick",size = 1) + 
    geom_point(data= arm_point[-nrow(arm_point),], aes(x=x,y=y),color = "firebrick",size = 6,shape = 16) +
    geom_point(data = data.frame(x = point[1],y=point[2]), aes(x=x,y=y),size = 5, shape = 18,color = "dodgerblue") + 
    geom_text(data = arm_text, aes(x = x, y = y , label = label),color = "white",size=5) + 
    labs (x = "", y = "") + 
    theme_classic() + 
    theme(axis.line = element_blank())
  return(arm_fig)
}

# function plotting the whole arm for all methods
plot_arm_joint <- function(arm1,arm2,arm3,point){
  require(ggplot2)
  max_radius = sum(arm1$len)
  min_radius = find_inner_radius(arm1$len)
  circ = circleFun(inner_radius = min_radius, radius = max_radius)
  arm_point1 = point_arm(arm1)
  arm_point2 = point_arm(arm2)
  arm_point3 = point_arm(arm3)
  arm_text1 = data.frame(label = seq(nrow(arm_point1)-1),x = arm_point1$x[-nrow(arm_point1)],y = arm_point1$y[-nrow(arm_point1)])
  arm_text2 = data.frame(label = seq(nrow(arm_point2)-1),x = arm_point2$x[-nrow(arm_point2)],y = arm_point2$y[-nrow(arm_point2)])
  arm_text3 = data.frame(label = seq(nrow(arm_point3)-1),x = arm_point3$x[-nrow(arm_point3)],y = arm_point3$y[-nrow(arm_point3)])
  arm_tip1 = create_arm_tip(arm_point1,arm1)
  arm_tip2 = create_arm_tip(arm_point2,arm2)
  arm_tip3 = create_arm_tip(arm_point3,arm3)
  arm_point1$x[nrow(arm_point1)] = arm_tip1$x1[1]
  arm_point1$y[nrow(arm_point1)] = arm_tip1$y1[1]
  arm_point2$x[nrow(arm_point2)] = arm_tip2$x1[1]
  arm_point2$y[nrow(arm_point2)] = arm_tip2$y1[1]
  arm_point3$x[nrow(arm_point3)] = arm_tip3$x1[1]
  arm_point3$y[nrow(arm_point3)] = arm_tip3$y1[1]
  arm_fig <- ggplot() + 
    geom_line(data = circ$out_up, aes(x=x,y=y)) + 
    geom_line(data = circ$in_up, aes(x=x,y=y)) + 
    geom_line(data = circ$out_low, aes(x=x,y=y)) + 
    geom_line(data = circ$in_low, aes(x=x,y=y)) + 
    geom_ribbon(data = circ$fill_up, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) +
    geom_ribbon(data = circ$fill_low, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) + 
    geom_path(data = arm_point1, aes(x=x,y=y,color = "GD"),size = 1) + 
    geom_path(data = arm_point2, aes(x=x,y=y,color = "F-R"),size = 1) + 
    geom_path(data = arm_point3, aes(x=x,y=y,color = "BFGS"),size = 1) + 
    geom_path(data = arm_tip1, aes(x=x1,y=y1,color = "GD"),size = 1) + 
    geom_path(data = arm_tip2, aes(x=x1,y=y1,color = "F-R"),size = 1) + 
    geom_path(data = arm_tip3, aes(x=x1,y=y1,color = "BFGS"),size = 1) + 
    geom_path(data = arm_tip1, aes(x=x2,y=y2,color = "GD"),size = 1) + 
    geom_path(data = arm_tip2, aes(x=x2,y=y2,color = "F-R"),size = 1) + 
    geom_path(data = arm_tip3, aes(x=x2,y=y2,color = "BFGS"),size = 1) + 
    geom_point(data= arm_point1[-nrow(arm_point1),], aes(x=x,y=y,color = "GD"),size = 6,shape = 16) +
    geom_point(data= arm_point2[-nrow(arm_point2),], aes(x=x,y=y,color = "F-R"),size = 6,shape = 16) +
    geom_point(data= arm_point3[-nrow(arm_point3),], aes(x=x,y=y,color = "BFGS"),size = 6,shape = 16) +
    geom_point(data = data.frame(x = point[1],y=point[2]), aes(x=x,y=y),size = 5, shape = 18,color = "hotpink") + 
    geom_text(data = arm_text1, aes(x = x, y = y , label = label),color = "white",size=5) + 
    geom_text(data = arm_text2, aes(x = x, y = y , label = label),color = "white",size=5) + 
    geom_text(data = arm_text3, aes(x = x, y = y , label = label),color = "white",size=5) + 
    labs (x = "", y = "",color="") + 
    theme_classic() + 
    theme(axis.line = element_blank())
  return(arm_fig)
}

# creating the configuration space in figures
circleFun <- function(center=c(0,0), inner_radius = 0, radius=1, npoints=200, start=0, end=2){
  tt1 <- seq(start*pi, end/2*pi, length.out=npoints)
  tt2 <- seq(end/2*pi, end*pi, length.out=npoints)
  out_up <- data.frame(
    x = center[1] + radius * cos(tt1),
    y = center[2] + radius * sin(tt1)
  )
  out_low <- data.frame(
    x = center[1] + radius * cos(tt2),
    y = center[2] + radius * sin(tt2)
  )
  in_up <- data.frame(
    x = center[1] + inner_radius * cos(tt1),
    y = center[2] + inner_radius* sin(tt1)
  )
  in_low <- data.frame(
    x = center[1] + inner_radius * cos(tt2),
    y = center[2] + inner_radius * sin(tt2)
  )
  x_up = c(out_up$x[out_up$x>in_up$x[1]], in_up$x ,out_up$x[out_up$x<in_up$x[npoints]])
  x_low = c(out_low$x[out_low$x<in_low$x[1]], in_low$x ,out_low$x[out_low$x>in_low$x[npoints]])
  suppressWarnings({
    fill_up <- data.frame(x = x_up,
                          ymax = sin(acos(x_up/radius))*radius,
                          ymin = ifelse((x_up<in_up$x[1]) & (x_up>in_up$x[npoints]), sin(acos(x_up/inner_radius))*inner_radius, 0))
    fill_low <- data.frame(x = x_low,
                           ymax = ifelse((x_low>in_low$x[1]) & (x_low<in_low$x[npoints]), sin(acos(x_low/inner_radius)+pi)*inner_radius, 0),
                           ymin = sin(acos(x_low/radius)+pi)*radius)
  }) 
  return(list(out_up = out_up, out_low = out_low, in_up = in_up, in_low = in_low, fill_up=fill_up,fill_low=fill_low))
}


# Visualization test 1
point = c(-7,0)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(4,1,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point=point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(4,1,1), angle=c(pi/4,pi/6,pi/10)), tol = 0.001, point=point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(4,1,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)


# Visualization test 2
point = c(5,3)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(5,2,1), angle=c(pi/2,pi/6,pi/10)),tol = 0.001, point=point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(5,2,1), angle=c(pi/2,pi/6,pi/10)), tol = 0.001, point=point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(5,2,1), angle=c(pi/2,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)


# Visualization test 3
point = c(1,1)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)



# Visualization test 4
point = c(3,2)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)



point = c(0,0)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)



point = c(10,10)
# Fletcher-Reeves
res1 = fletcher_reeves(arm=data.frame(len = c(3,2,1,1,2,2,2,2), angle=pi*c(3,2,1,1,2,2,2,2)),tol = 0.001, point = point)
plot_arm(res1$arm,point = point)
# gradient descent
res2 = gradient_descent(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res2$arm,point = point)
# BFGS
res3 = BFGS(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res3$arm,point = point)
# all
plot_arm_joint(res2$arm,res1$arm,res3$arm,point = point)


# auto generate problems 
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

# classifying convergence of one problem
test_conv <- function(grad_f, eval_f,tol_1 = 1e-3, tol_2 = 1e-3){
  if ((grad_f < tol_1) & (eval_f < tol_2)){
    return(1)
  }else if (grad_f < tol_1){
    return(2)
  }else{
    return(0)
  }
}

# a simple simulation run 
# to run on server use the run_stats.R 
# DO NOT RUN: if you do not have 1 hour atleast to wait.
run <- function(){
  n_arms = c(3,5,7,9)
  n_iter = 1000
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
      res_gd[k,] = ele$res_fr
      k = k+1
    }
  }
  return(list(res_bfgs = res_bfgs, res_gd = res_gd, res_fr = res_fr))
}


res = run()

# generate a dataframe of statistics from a simulation
get_stats <- function(res,tot_it){
  stats = data.frame(
    n = as.factor(c(sapply(seq(3,9,2),function(x){rep(x,3)}))),
    method = rep(c("GD","F-R","BFGS"),4),
    conv = rep(NA,3*4),
    mean_it = rep(NA,3*4),
    sd_it = rep(NA,3*4),
    mean_run = rep(NA,3*4),
    sd_run = rep(NA,3*4)
  )
  i = 1
  for (n in seq(3,10,2)){
    stats$conv[i] = sum(res$res_gd$conv[res$res_gd$n==n]!=0)/tot_it
    stats$conv[i+1] = sum(res$res_fr$conv[res$res_fr$n==n]!=0)/tot_it
    stats$conv[i+2] = sum(res$res_bfgs$conv[res$res_bfgs$n==n]!=0)/tot_it
    stats$mean_it[i] = mean(res$res_gd$steps[res$res_gd$n==n])
    stats$mean_it[i+1] = mean(res$res_fr$steps[res$res_fr$n==n])
    stats$mean_it[i+2] = mean(res$res_bfgs$steps[res$res_bfgs$n==n])
    stats$sd_it[i] = sd(res$res_gd$steps[res$res_gd$n==n])
    stats$sd_it[i+1] = sd(res$res_fr$steps[res$res_fr$n==n])
    stats$sd_it[i+2] = sd(res$res_bfgs$steps[res$res_bfgs$n==n])
    stats$mean_run[i] = mean(res$res_gd$runtime[res$res_gd$n==n])
    stats$mean_run[i+1] = mean(res$res_fr$runtime[res$res_fr$n==n])
    stats$mean_run[i+2] = mean(res$res_bfgs$runtime[res$res_bfgs$n==n])
    stats$sd_run[i] = sd(res$res_gd$runtime[res$res_gd$n==n])
    stats$sd_run[i+1] = sd(res$res_fr$runtime[res$res_fr$n==n])
    stats$sd_run[i+2] = sd(res$res_bfgs$runtime[res$res_bfgs$n==n])
    i = i + 3
  }
  return(stats)
}

load(file = "./optim-stats-large.Rdata")
stats # print statistics from simulation
# create iteration plot in project paper
it_plot <- ggplot(stats,aes(x = n)) + 
  geom_ribbon(aes(ymin = log(it_q25), ymax = log(it_q75), group = method, fill=method, color = method),alpha = 0.2) + 
  geom_line(aes(y = log(it_q50),group = method,color = method),size = 1) + 
  theme_bw()+ 
  labs(x = "# segments", y = "log iterations",fill="",color="")
it_plot

# create runtime plot in project paper
run_plot <- ggplot(stats,aes(x = n)) + 
  geom_ribbon(aes(x = n, ymin = log(run_q25), ymax = log(run_q75), group = method,fill = method, color = method),alpha = 0.2) + 
  geom_line(aes(x = n, y = log(run_q50),group = method, color = method),size = 1) +
  labs(x = "# segments", y = "log runtime", color = "", fill="") + 
  theme_bw()
run_plot


library(ggpubr) #required for join plot of iteration and runtime
ptot <- ggarrange(it_plot, run_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ptot

