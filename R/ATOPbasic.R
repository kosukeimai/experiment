#### functions for calculating bounds

f1=function(x,y,u){u-(u-x)/y }
f2=function(x,y,l){l+(x-l)/y }

fLB= function (x1,y1,x0,y0,l,u){
  ifelse(y1>0,  max(l,f1(x1,y1,u)),l)-ifelse(y0>0,min(u,f2(x0,y0,l)),u)
}

fUB= function (x1,y1,x0,y0,l,u){
  ifelse(y1>0,  min(u,f2(x1,y1,l)),u)-ifelse(y0>0,max(l,f1(x0,y0,u)),l)
}

###  functions for calculating CI


CalC= function (x,hat.LB,hat.UB,sigma.LB,sigma.UB,alpha){
  pnorm(x+(hat.UB-hat.LB)/max(sigma.LB,sigma.UB))-pnorm(-x)-(1-alpha)
}