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
    u <- Q %*% U
    
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



# Fortran code
rand.svd_F <- function(A, k, q, compute.u, compute.vt)
{
  method <- 'R'
  
  m <- as.integer(nrow(A))
  n <- as.integer(ncol(A))
  
  if (!is.double(A))
    storage.mode(A) <- "double"
  
  if (compute.u)
    compute.u <- 'V'
  else
    compute.u <- 'N'
  
  if (compute.vt)
    compute.vt <- 'V'
  else
    compute.vt <- 'N'
  
  ret <- .Call("R_randsvd", method, m, n, A, as.integer(k), as.integer(q), compute.u, compute.vt)
  
  return( ret )
}



# Unified method
rand.svd <- function(A, k=1, q=3, compute.u=TRUE, compute.vt=TRUE, method="Fortran")
{
  ### input parameter checking
  
  # check A
  if (!is.numeric(A)) 
    stop("'A' must be numeric")
  if (any(is.na(A)))
    stop("missing values not allowed")
  if (any(!is.finite(A))) 
    stop("infinite values not allowed")
  if (!is.matrix(A))
    dim(A) <- c(length(A), 1L)
  
  # check k
  if (!is.int(k) || k<0)
    stop("'k' must be a positive integer")
  
  k <- as.integer(k)
  
  # check q
  if (!is.int(q) || q<0)
    stop("'q' must be a positive integer")
  if (k > nrow(A))
    stop("'k' must be no greater than nrow(A)")
  
  q <- as.integer(q)
  
  # check compute.u and compute.vt
  if (!is.logical(compute.u))
    stop("'compute.u' must be logical")
  if (!is.logical(compute.vt))
    stop("'compute.vt' must be logical")
  
  # check method
  method <- match.arg(tolower(method), c("fortran", "r"))
  
  
  ### Call the appropriate function
  if (method == "fortran")
    ret <- rand.svd_F(A=A, k=k, q=q, compute.u=compute.u, compute.vt=compute.vt)
  else
    ret <- rand.svd_R(A=A, k=k, q=q, compute.u=compute.u, compute.vt=compute.vt)
  
  
  return( ret )
}


