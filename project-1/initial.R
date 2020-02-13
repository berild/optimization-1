library(ggplot2)
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

plot_arm <- function(arm,point){
  require(ggplot2)
  max_radius = sum(arm$len)
  min_radius = find_inner_radius(arm$len)
  circ = circleFun(inner_radius = min_radius, radius = max_radius)
  arm_point=point_arm(arm)
  arm_fig <- ggplot() + 
    geom_path(data = arm_point, aes(x=x,y=y),color = "firebrick",size = 1) + 
    geom_point(data = arm_point[nrow(arm_point),], aes(x=x,y=y),color = "firebrick",size = 5,shape=15) + 
    geom_point(data= arm_point[-nrow(arm_point),], aes(x=x,y=y),color = "firebrick",size = 3) +
    geom_point(data = data.frame(x = point[1],y=point[2]), aes(x=x,y=y),size = 5, shape = 18,color = "dodgerblue") + 
    geom_line(data = circ$out_up, aes(x=x,y=y)) + 
    geom_line(data = circ$in_up, aes(x=x,y=y)) + 
    geom_line(data = circ$out_low, aes(x=x,y=y)) + 
    geom_line(data = circ$in_low, aes(x=x,y=y)) + 
    geom_ribbon(data = circ$fill_up, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) +
    geom_ribbon(data = circ$fill_low, aes(x=x,ymin=ymin,ymax=ymax),alpha = 0.2) + 
    labs (x = "", y = "") + 
    theme_bw()
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

point = c(1,1)
arm = BFGS(arm=data.frame(len = c(1,4,1), angle=c(pi/4,pi/6,pi/10)),tol = 0.001, point = point)
plot_arm(arm,point = point)





