#ffsvd <- function(x, nu=1, nv=1, method="best")
#{
#  method <- match.arg(tolower(method), c("best", "exact", "fastmap", "randsvd"))
#  
#  
#  
#  if (method=="best"){
#    
#    
#    
#    method <- 
#  }
#  
#  
#  if (method=="fastmap")
#    
#  else if (method=="randsvd")
#    
#  else if (method=="exact")
#    return(svd(x, nu=nu, nv=nv, LINPACK=FALSE))
#}


#ffpca <- function(x, retx=TRUE, center=TRUE, scale=FALSE, method="best", tol=NULL, ...)
#{
#  method <- match.arg(tolower(method), c("best", "fastmap", "randsvd", "exact"))
#  
#  
#  
#  if (method=="exact")
#    return(prcomp(x=x, retx=retx, center=center, scale.=scale, tol=tol))
#  
#  
#  
#}
