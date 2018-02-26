##### basic functions for calculating the CADE

Difflist=function(a,b){
  J=length(a)
  c=a
  for (j in 1:J){
    c[[j]]=a[[j]]-b[[j]]
  }
  return(c)
}

Meanlist=function(a){
  J=length(a)
  s=0
  for(j in 1:J){
    s=s+mean(a[[j]])
  }
  return(s/J)
}

Productlist=function(a,b){
  J=length(a)
  c=a
  for (j in 1:J){
    c[[j]]=a[[j]]*b[[j]]
  }
  return(c)
}
