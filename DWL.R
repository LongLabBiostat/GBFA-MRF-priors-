

dwl.new <- function(X=NULL,y=NULL,XX=NULL,Xy=NULL,yy=NULL)
{
  dwl <- new.env(parent=emptyenv())
  
  if ( !is.null(X) & !is.null(y) )
  {
    dwl$X <- X
    dwl$y <- as.vector(y)
    
    dwl$n <- nrow(X)
    dwl$p <- ncol(X)
    dwl$m <- min(dwl$n,dwl$p)
    
    dwl$XxMax <- min(3*dwl$n,dwl$p)
    dwl$Xx <- matrix(0,dwl$p,dwl$XxMax)
    dwl$XxIdx <- rep(0,dwl$p)
    dwl$XxCnt <- 0
    
    dwl$Xy <- as.vector(t(X) %*% y)
    dwl$yy <- sum(y^2)
  }
  else if ( !is.null(XX) & !is.null(Xy) )
  {
    dwl$p <- ncol(XX)
    dwl$m <- dwl$p
    
    dwl$XxMax <- dwl$p
    dwl$Xx <- as.matrix(XX)
    dwl$XxIdx <- 1:dwl$p
    dwl$XxCnt <- dwl$p
    
    dwl$Xy <- as.vector(Xy)
    dwl$yy <- as.numeric(yy)
  }
  else
  {
    message("X and y or XX and Xy must be provided.")
    stop()
  }
  
  dwl$lam <- rep(max(abs(dwl$Xy))*1.2,dwl$p)
  
  dwl$Idx <- rep(0,dwl$p)
  dwl$A <- numeric()
  dwl$nA <- 0
  
  dwl$B <- numeric()
  dwl$S <- numeric()
  
  dwl$C <- dwl$Xy
  dwl$iXXa <- matrix(0,dwl$m,dwl$m)
  
  dwl$coef = rep(0,dwl$p)
  dwl$sign = rep(0,dwl$p)
  
  dwl
}


dwl.lookahead <- function(dwl,tlam)
{
  ######################
  # Find the direction #
  ######################
  
  dwl$dlam <- tlam - dwl$lam
  
  # calculate dwl$dB/dalpha and dC/dalpha
  if ( dwl$nA > 0 )
  {
    dwl$dB <- as.vector(-dwl$iXXa[1:dwl$nA,1:dwl$nA] %*% (dwl$S * dwl$dlam[dwl$A]))
    dC <- as.vector(-dwl.getXXI(dwl,dwl$A) %*% dwl$dB)
  }
  else
    dC <- rep(0,dwl$p)
  
  
  ######################
  # How far can we go? #
  ######################
  
  # find breakpoint
  dwl$alpha <- 1
  dwl$type <- 0
  dwl$idx <- 0
  
  if ( dwl$nA > 0 )
  {
    pbp0 = -dwl$B/dwl$dB
    for ( l in 1:dwl$nA )
      if ( (dwl$B[l]+dwl$dB[l])*dwl$S[l] < 0 & pbp0[l] < dwl$alpha )
      {
        dwl$alpha <- pbp0[l]
        dwl$type <- 0
        dwl$idx <- dwl$A[l]
      }
  }
  
  pbp1 = (dwl$lam-dwl$C)/(dC-dwl$dlam)
  pbp2 = -(dwl$lam+dwl$C)/(dC+dwl$dlam)
  for ( k in 1:dwl$p )
    if ( dwl$Idx[k] == 0 )
    {
      if ( dwl$C[k]+dC[k] > dwl$lam[k]+dwl$dlam[k] & pbp1[k] < dwl$alpha )
      {
        dwl$alpha <- pbp1[k]
        dwl$type <- 1
        dwl$idx <- k
      }
      if ( dwl$C[k]+dC[k] < -dwl$lam[k]-dwl$dlam[k] & pbp2[k] < dwl$alpha )
      {
        dwl$alpha <- pbp2[k]
        dwl$type <- -1
        dwl$idx <- k
      }
    }
  
}


dwl.updateBC <- function(dwl)
{
  dwl$coef <- rep(0,dwl$p)
  if ( dwl$nA > 0 )
  {
    dwl$B <- as.vector(dwl$iXXa[1:dwl$nA,1:dwl$nA] %*% ( dwl$Xy[dwl$A] - dwl$S*dwl$lam[dwl$A] ))
    dwl$C <- dwl$Xy - as.vector(dwl.getXXI(dwl,dwl$A) %*% dwl$B)
    dwl$coef[dwl$A] <- dwl$B
  }
  else
  {
    dwl$B <- numeric()
    dwl$C <- dwl$Xy
  }
  dwl$sign <- sign(dwl$coef)
}


dwl.fit <- function(dwl,tlam)
{
  # Eat free lunch
  for ( i in 1:dwl$p )
  {
    if ( dwl$Idx[i] == 0 & dwl$lam[i] < tlam[i] )
      dwl$lam[i] <- tlam[i]
  }
  
  niter = 0
  repeat
  {
    niter = niter + 1
    
    ##############
    # Look ahead #
    ##############
    
    lh = dwl.lookahead(dwl,tlam)
    
    
    ############################
    # Add or remove a variable #
    ############################
    
    if ( dwl$alpha < 1 )
      if ( dwl$type == 0 )
        dwl.remove(dwl,dwl$idx)
    else
    {
      dwl$B <- dwl$B + dwl$alpha*dwl$dB
      dwl.add(dwl,dwl$idx,dwl$type)
    }
    
    
    ##########
    # Update #
    ##########
    
    dwl$lam <- dwl$lam + dwl$dlam*dwl$alpha
    dwl.updateBC(dwl)
    
    
    if ( dwl$alpha ==  1 )
      break
  }
  
  
  list(coef=dwl$coef,sign=dwl$sign,niter=niter)
}



dwl.add = function(dwl,k,sgn)
{
  b = dwl.getXXI(dwl,k)
  
  if ( dwl$nA > 0 )
  {
    a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
    del = drop(a %*% b[dwl$A])
    d = drop(b[k] - crossprod(del,b[dwl$A]))
    
    if ( d < 1e-8 )
    {
      # message("Warning: numerical instability")
      
      pos = which.max(del*sgn/dwl$B)
      dwl.remove(dwl,dwl$A[pos])
      
      if ( dwl$nA > 0 )
      {
        a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
        del = drop(a %*% b[dwl$A])
        d = drop(b[k] - crossprod(del,b[dwl$A]))
      }
    }
  }
  
  # Now add k
  if ( dwl$nA > 0 )
  {
    dwl$iXXa[1:dwl$nA,1:dwl$nA] <- a + del %*% t(del) / d
    dwl$iXXa[1:dwl$nA,dwl$nA+1] <- -del / d
    dwl$iXXa[dwl$nA+1,1:dwl$nA] <- -del / d
    dwl$iXXa[dwl$nA+1,dwl$nA+1] <- 1/d
  }
  else
  {
    dwl$iXXa[1] <- 1/b[k]
  }
  
  dwl$nA <- dwl$nA+1
  dwl$Idx[k] <- dwl$nA
  dwl$A <- c(dwl$A,k)
  dwl$S <- c(dwl$S,sgn)
}



dwl.remove = function(dwl,k)
{
  l = dwl$Idx[k]
  dwl$Idx[k] <- 0
  if ( l<dwl$nA )
    dwl$Idx[dwl$A[(l+1):dwl$nA]] <- dwl$Idx[dwl$A[(l+1):dwl$nA]] - 1
  dwl$A <- dwl$A[-l]
  dwl$S <- dwl$S[-l]
  
  if ( dwl$nA>1 )
  {
    a = dwl$iXXa[1:dwl$nA,1:dwl$nA]
    b = a[,l]
    dwl$iXXa[1:(dwl$nA-1),1:(dwl$nA-1)] <- a[-l,-l] - b[-l] %*% t(b[-l]) / b[l]
  }
  dwl$iXXa[,dwl$nA] <- 0
  dwl$iXXa[dwl$nA,] <- 0
  
  dwl$nA <- dwl$nA-1
}



dwl.getXXI <- function(dwl,I)
{
  for ( k in I )
    if ( dwl$XxIdx[k] == 0 )
    {
      dwl$XxCnt <- dwl$XxCnt + 1
      if ( dwl$XxCnt > dwl$XxMax )
      {
        oldmax = dwl$XxMax
        oldXx = dwl$Xx
        dwl$XxMax <- min(oldmax*2,dwl$p)
        dwl$Xx <- matrix(0,dwl$p,dwl$XxMax)
        dwl$Xx[,1:oldmax] <- oldXx
      }
      dwl$XxIdx[k] <- dwl$XxCnt
      dwl$Xx[,dwl$XxCnt] <- t(dwl$X) %*% dwl$X[,k]
    }
  
  as.matrix(dwl$Xx[,dwl$XxIdx[I]])
}


