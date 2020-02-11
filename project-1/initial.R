library(ggplot2)
init_arm <- function(len,angle){
  arm_polar = data.frame(len = len, angle = angle)
  arm_cart = point_arm(arm_polar)
  return(list(polar = arm_polar,cart = arm_cart))
}

point_arm <- function(arm){
  x = y = rep(0,nrow(arm)+1)
  for (i in seq(nrow(arm))){
    x[i+1] = cos(arm$angle[i])*arm$len[i] + x[i]
    y[i+1] = sin(arm$angle[i])*arm$len[i] + y[i]
  }
  return(data.frame(x = x, y = y))
}

find_inner_radius <-function(len){
  if (max(len) > sum(len[-which.max(len)])){
    return(max(len) > sum(len[-which.max(len)]))
  }else{
    return(0)
  }
}

plot_arm <- function(arm){
  require(ggplot2)
  max_radius = sum(arm$polar$len)
  min_radius = find_inner_radius(arm$polar$len)
  circ = circleFun(inner_radius = min_radius, radius = max_radius)
  arm_fig <- ggplot() + 
    geom_line(data = arm$cart, aes(x=x,y=y),color = "firebrick",size = 1) + 
    geom_point(data= arm$cart, aes(x=x,y=y),color = "firebrick",size = 3) +
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

circleFun <- function(center=c(0,0), inner_radius = 0, radius=1, npoints=500, start=0, end=2){
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
  fill_up <- data.frame(x = out_up$x,
                        ymax = out_up$y,
                        ymin = ifelse((out_up$x>in_up$x[1]) | (out_up$x<in_up$x[npoints]),0,sin(acos(out_up$x/inner_radius))*inner_radius))
  fill_low <- data.frame(x = out_low$x,
                         ymin = out_low$y,
                         ymax = ifelse((out_low$x<in_low$x[1]) | (out_low$x>in_low$x[npoints]),0,sin(acos(out_low$x/inner_radius)+pi)*inner_radius))
  return(list(out_up = out_up, out_low = out_low, in_up = in_up, in_low = in_low, fill_up=fill_up,fill_low=fill_low))
}

arm_optim <- function(init,p){
  
}

arm = init_arm(c(2,4,1),c(pi/4,pi/6,pi/10))
plot_arm(arm)
