library(ggplot2)
# gradient descent
# conjugate gradient
# newtons method
# bfgs
# nelder-mead

init_arm <- function(len,angle){
  arm_polar = data.frame(len = len, angle = angle)
  arm_cart = point_arm(arm_polar)
  return(list(polar = arm_polar,cart = arm_cart))
}



point_arm <- function(arm){
  x = y = rep(0,nrow(arm)+1)
  for (i in seq(nrow(arm))){
    x[i+1] = cos(sum(arm$angle[1:i]))*arm$len[i] + x[i]
    y[i+1] = sin(sum(arm$angle[1:i]))*arm$len[i] + y[i]
  }
  return(data.frame(x = x, y = y))
}

find_inner_radius <-function(len){
  if (max(len) > sum(len[-which.max(len)])){
    return(max(len) - sum(len[-which.max(len)]))
  }else{
    return(0)
  }
}

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
    #geom_point(data = arm_point[nrow(arm_point),], aes(x=x,y=y),color = "firebrick",size = 5,shape = 15) + 
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


mag <- function(x){
  sqrt(sum(x^2))  
}

point = c(3,2)
res = fletcher_reeves(arm=data.frame(len = c(3,2,2), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point=point)
plot_arm(res$arm,point = point)
res = gradient_descent(arm=data.frame(len = c(3,2,2), angle=c(pi/4,pi/6,pi/10)), tol = 0.001, point=point)
plot_arm(res$arm,point = point)
res = BFGS(arm=data.frame(len = c(3,2,2), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)

point = c(1,1)
res = fletcher_reeves(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = gradient_descent(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = BFGS(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)

point = c(3,2)
res = fletcher_reeves(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = gradient_descent(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = BFGS(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/3)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)

point = c(0,0)
res = fletcher_reeves(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = gradient_descent(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)
res = BFGS(arm=data.frame(len = c(3,2,1,1), angle=c(pi/4,pi/6,pi/10,pi/15)),tol = 0.001, point = point)
plot_arm(res$arm,point = point)

generate_problem <- function(n){
  len = sample(seq(15),n,replace = T)
  angle = runif(n,min = 0, max = 2*pi)
  max_radius = sum(len)
  min_radius = find_inner_radius(len)
  tmp_r = runif(1,min=min_radius, max = max_radius)
  tmp_angle = runif(1,min = 0, max = 2*pi)
  point = c (tmp_r*cos(tmp_angle),tmp_r*sin(tmp_angle))
  return(list(point = point, arm = data.frame(angle = angle, len = len)))
}

run <- function(n){
  prob = generate_problem(n)
  res = fletcher_reeves(arm = prob$arm,tol = 0.001, point = prob$point)
  return(plot_arm(res$arm,point = prob$point))
}

run(3)
