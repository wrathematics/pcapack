# Copyright 2013, Schmidt and Ostrouchov


# R code
rand.svd_R <- function(A, k, q, compute.u, compute.vt)
{
  q <- as.integer(q)
  n <- ncol(A)
  
  ### Stage A
  Omega <- matrix(rnorm(n*2L*k), nrow=n, ncol=2L*k)
  Y <- A %*% Omega
  Q <- qr.Q(qr(Y))
  tA <- t(A)
  
  for (i in 1:q)
  {
    Y <- tA %*% Q
    Q <- qr.Q(qr(Y))
    Y <- A %*% Q
    Q <- qr.Q(qr(Y))
  }
  
  
  ### Stage B
  B <- t(Q) %*% A
  
  if (!compute.u)
    nu <- 0
  else
    nu <- min(nrow(B), ncol(B))
  
  if (!compute.vt)
    nv <- 0
  else
    nv <- min(nrow(B), ncol(B))
  
  svd.B <- La.svd(x=B, nu=nu, nv=nv)
  
  
  # Produce u/vt as desired
  if (compute.u)
  {
    u <- svd.B$u
    u <- Q %*% u
    
    d <- svd.B$d
    
    d <- d[1L:k]
    u <- u[, 1L:k]
  }
  
  if (compute.vt)
  {
    vt <- svd.B$vt[1L:k, ]
  }
  
  # wrangle return
  if (compute.u)
  {
    if (compute.vt)
      svd <- list(d=d, u=u, vt=vt)
    else
      svd <- list(d=d, u=u)
  }
  else
  {
    if (compute.vt)
      svd <- list(d=d, vt=vt)
    else
      svd <- list(d=d)
  }
  
  return( svd )
}



# Unified method
rsvd <- function(x, k=1, q=3, ret.u=TRUE, ret.vt=TRUE)
{
  ### Cheap checks first
  assert.type(ret.u, "logical")
  assert.type(ret.vt, "logical")
  
  assert.natnum(k)
  k <- as.integer(k)
  if (k > nrow(x))
    stop("'k' must be no greater than nrow(x)")
  
  assert.natnum(q)
  q <- as.integer(q)
  
  assert.type(x, "numeric")
    ### TODO do this in one pass, idiot
#  if (any(is.na(A)))
#    stop("missing values not allowed")
#  if (any(!is.finite(A))) 
#    stop("infinite values not allowed")
  if (!is.matrix(x))
    dim(A) <- c(length(x), 1L)
  
  .Call("R_pcapack_randsvd", )
}


