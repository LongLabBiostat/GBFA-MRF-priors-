# This file contains the function GBFA_EM.
# written by Changgee Chang
# version: 20180711
#

#
# function GBFA_EM: runs EM for GBFA with Spike and Slab plus MRF prior
#
# Parameters:
# X: p by n data matrix where p is the number of covariates and n is the number of samples
# type: p by 1 vector of data types (see below)
# param: p by 1 vector of parameters nu, n, or r (see data types below)
# E: e by 2 matrix for edges where e is the number of edges.
#    E must be sorted by the first column and then the second column.
#    Both (j,k) and (k,j) must be added and therefore e is twice the actual number of edges.
# L: number of factors
# v0,v1: variances of gaussian spike and slab for W. For SSL, lam0=1/v0 and lam1=1/v1.
# delta: shrinkage parameter in MRF prior
# eta: smoothing parameter in MRF prior
# xi: shrinkage parameter for mu_Z
# fit.m: logical; whether to fit m or not.
# m.init: initial value for m; 0: zero, 1: trimmed mean(default), 2: median, 3: mean, p by 1 vector: manual
# scale: logical; find the scale g or not. If not, set to 1.
# sstype: spike and slab type; "lasso" or "gaussian"
# smoothing: "Ising" or "MRF"
# PXL: If TRUE, use the PXL-EM algorithm.
# W.init: initial value of W
#
# Data types:
# 0: Gaussian - nu_j is required as param
# 1: Binomial - n_j must be provided as param
# 2: Negative Binomial - r_j must be provided as param
# 3: Poisson


source("E:/Post_doc/BFA_01/MCMC/DWL.R")

GBFA_EM <- function(X,type,param,E,L,v0,v1,delta,eta,xi=0,fit.m=FALSE,m.init=1,scale=FALSE,sstype="lasso",smoothing="Ising",PXL=FALSE,W.init=NULL,eps=1e-3,maxIter=500)
{
  p = nrow(X)
  n = ncol(X)
  
  # edge
  e = nrow(E)
  Eidx = rep(1,p+1)
  for ( j in 1:p )
  {
    Eidx[j+1] = Eidx[j]
    while ( Eidx[j+1] <= e )
    {
      if ( E[Eidx[j+1],1] == j )
        Eidx[j+1] = Eidx[j+1] + 1
      else
        break
    }
  }
  
  # initialization
  psi = matrix(0,p,n)
  kappa = matrix(0,p,n)
  b = matrix(4,p,n)
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      psi[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      kappa[j,] = X[j,]-param[j]/2
      b[j,] = param[j]
    }
    if ( type[j] == 2 )
    {
      kappa[j,] = (X[j,]-param[j])/2
      b[j,] = param[j]+X[j,]
    }
    if ( type[j] == 3 )
    {
      kappa[j,] = X[j,]
    }
  }
  
  Poidx = which(type==3)
  if ( length(Poidx) > 0 )
    Poisson.exist = TRUE
  else
    Poisson.exist = FALSE
  
  
  # initialize Y
  Y = init(X,type,param)
  
  # initialize m
  if ( length(m.init) == p )
    m = m.init
  else if ( m.init == 0 )
    m = rep(0,p)
  else if ( m.init == 1 )
    m = apply(Y,1,function(x) mean(x,0.2))
  else if ( m.init == 2 )
    m = apply(Y,1,median)
  else if ( m.init == 3 )
    m = apply(Y,1,meam)
  else
    m = rep(0,p)
  
  # find scale
  if ( scale )
    g = sqrt(apply((Y-m)^2,1,mean))
  else
    g = rep(1,p)
  
  # W and Z init
  WZinit = (Y-m)/g
  if ( !is.null(W.init) )
  {
    W = W.init
    WW = t(W)%*%W
    Z = chol2inv(chol(WW)) %*% t(W) %*% WZinit
  }
  else
  {
    init_svd = svd(WZinit,L,L)
    W = init_svd$u %*% diag(init_svd$d[1:L],L) / sqrt(n)
    Z = t(init_svd$v) * sqrt(n)
#    W.varimax = varimax(W)
#    W = W %*% W.varimax$rotmat
#    Z = t(W.varimax$rotmat) %*% Z
  }
  SigmaZ = array(diag(1,L),c(L,L,n))
  
  GW = W*g
  GWZ = GW%*%Z
  mu = m + GWZ
  gWSigmaZ = array(GW%*%matrix(SigmaZ,L),c(p,L,n))
  g2WWSigmaZ = apply(gWSigmaZ*rep(GW,n),c(1,3),sum)

  
  rho = b/4
  rho[Poidx,] = exp(mu[Poidx,]+g2WWSigmaZ[Poidx,]/2)

  alpha = matrix(0,p,L)
  theta = matrix(0.5,p,L)
  
  # for SSL
  lam0 = 1/v0
  lam1 = 1/v1
  dlam = lam1-lam0
  llam0 = log(lam0)
  llam1 = log(lam1)
  dllam = llam1-llam0
  
  iv0 = 1/v0
  iv1 = 1/v1
  div = iv1-iv0
  lv0 = log(v0)
  lv1 = log(v1)
  dlv = lv1-lv0
  
  
  # storage
  eta = eta[order(eta)]
  if ( eta[1] != 0 )
    eta = c(0,eta)
  nEta = length(eta)
  
  Ws = array(0,c(p,L,nEta))
  Zs = array(0,c(L,n,nEta))
  rhos = array(0,c(p,n,nEta))
  thetas = array(0,c(p,L,nEta))
  LLKs = rep(0,nEta)
  BICs = rep(0,nEta)
  iters = rep(0,nEta)
  
  for ( etai in 1:nEta )
  {
    s_Z = rep(1,n)
    s_SigmaZ = rep(1,n)
    s_W = rep(1,p)
    iter = 0
    
    while ( iter < maxIter )
    {
      pW = W
      iter = iter + 1
      print(iter)
      
      # E-step for gamma
      for ( l in 1:L )
      {
        dir = rep(0,p)
        for ( j in 1:p )
          if ( Eidx[j] < Eidx[j+1] )
            if ( smoothing == "MRF" )
              dir[j] = eta[etai]*sum(theta[E[Eidx[j]:(Eidx[j+1]-1),2],l])
            else
              dir[j] = eta[etai]*sum(2*theta[E[Eidx[j]:(Eidx[j+1]-1),2],l]-1)
            if ( sstype == "Gaussian" )
              dir = dir - delta - dlv/2 - W[,l]^2*div/2 - alpha[,l]
            else
              dir = dir - delta + dllam - abs(W[,l])*dlam - alpha[,l]
            grad = theta[,l]*(1-theta[,l])*dir
            dg = sum(dir*grad)
            
            stheta = sum(theta[,l])
            if ( smoothing == "MRF" )
              Q = eta[etai]*sum(theta[E[,1],l]*theta[E[,2],l])/2
            else
              Q = eta[etai]*sum(2*theta[E[,1],l]*theta[E[,2],l]-theta[E[,1],l]-theta[E[,2],l]+1)/2
            if ( sstype == "Gaussian" )
              Q = Q - stheta*dlv/2 - sum(W[,l]^2*theta[,l])*div/2 - delta*stheta
            else
              Q = Q + stheta*dllam - sum(abs(W[,l])*theta[,l])*dlam - delta*stheta
            
            s = 1
            while (TRUE)
            {
              newalpha = alpha[,l] + s*dir
              newtheta = 1/(1+exp(-newalpha))
              newstheta = sum(newtheta)
              if ( smoothing == "MRF" )
                newQ = eta[etai]*sum(newtheta[E[,1]]*newtheta[E[,2]])/2
              else
                newQ = eta[etai]*sum(2*newtheta[E[,1]]*newtheta[E[,2]]-newtheta[E[,1]]-newtheta[E[,2]]+1)/2
              if ( sstype == "Gaussian" )
                newQ = newQ - newstheta*dlv/2 - sum(W[,l]^2*newtheta)*div/2 - delta*newstheta
              else
                newQ = newQ + newstheta*dllam - sum(abs(W[,l])*newtheta)*dlam - delta*newstheta
              if ( newQ-Q >= s*dg/2 )
              {
                alpha[,l] = newalpha
                theta[,l] = newtheta
                break
              }
              s = s / 2
            }
      }
      
      # E-step for Z
      if ( Poisson.exist )
      {
        # E-step for Sigma_Z
        cSigmaZ = array(apply(SigmaZ,3,chol),c(L,L,n))
        iSigmaZ = array(apply(cSigmaZ,3,chol2inv),c(L,L,n))
        Q = apply(apply(SigmaZ,3,diag),2,sum)/2 - apply(log(apply(cSigmaZ,3,diag)),2,sum)
        Q = Q + apply(rho[Poidx,],2,sum)
        Q = Q + apply(g2WWSigmaZ[-Poidx,]*rho[-Poidx,],2,sum)/2
        grad = (array(t(GW)%*%(GW[,rep(1:L,n)]*rho[,rep(1:n,each=L)]),c(L,L,n)) + rep(diag(1,L),n) - iSigmaZ)/2
        dir = -(grad+array(apply(grad,3,t),c(L,L,n)))/2
        dg = apply(dir*grad,3,sum)
        gWdir = array(GW%*%matrix(dir,L),c(p,L,n))
        g2WWdir = apply(gWdir*rep(GW,n),c(1,3),sum)
        for ( i in 1:n )
        {
          while (TRUE)
          {
            newSigmaZ = SigmaZ[,,i] + s_SigmaZ[i]*dir[,,i]
            noerr = TRUE
            tryCatch(chol(newSigmaZ),error=function(e){noerr<<-FALSE})
            if ( noerr )
              break
            s_SigmaZ[i] = s_SigmaZ[i] / 2
          }
          
          while (TRUE)
          {
            newSigmaZ = SigmaZ[,,i] + s_SigmaZ[i]*dir[,,i]
            cnewSigmaZ = chol(newSigmaZ)
            newQ = sum(diag(newSigmaZ))/2 - sum(log(diag(cnewSigmaZ)))
            newrho = rho[Poidx,i]*exp(s_SigmaZ[i]*g2WWdir[Poidx,i]/2)
            newQ = newQ + sum(newrho)
            newQ = newQ + sum((g2WWSigmaZ[-Poidx,i]+s_SigmaZ[i]*g2WWdir[-Poidx,i])*rho[-Poidx,i])/2
            if ( Q[i]-newQ >= -s_SigmaZ[i]*dg[i]/2 )
            {
              SigmaZ[,,i] = newSigmaZ
              rho[Poidx,i] = newrho
              gWSigmaZ[,,i] = gWSigmaZ[,,i] + s_SigmaZ[i]*gWdir[,,i]
              g2WWSigmaZ[,i] = g2WWSigmaZ[,i] + s_SigmaZ[i]*g2WWdir[,i]
              s_SigmaZ[i] = s_SigmaZ[i] * 1.2
              break
            }
            s_SigmaZ[i] = s_SigmaZ[i] / 2
          }
        }
        
  
        # E-step for mu_Z
        C = rho*(mu-psi)
        C2 = C*(mu-psi)/2
        C[Poidx,] = rho[Poidx,]
        C2[Poidx,] = rho[Poidx,]
        Z2 = apply(Z*Z,2,sum)
        Q = Z2/2 + apply(C2-kappa*mu,2,sum)
        grad = Z + t(GW)%*%(C-kappa)
        
        if ( xi == 0 )
        {
          H = array(t(GW)%*%(GW[,rep(1:L,n)]*rho[,rep(1:n,each=L)]),c(L,L,n)) + rep(diag(1,L),n)
          for ( i in 1:n )
          {
            dir = -chol2inv(chol(H[,,i]))%*%grad[,i]
            dg = crossprod(dir,grad[,i])
            
            dir2 = crossprod(dir,dir)
            dirZ = crossprod(dir,Z[,i])
            GWdir = GW%*%dir
    
            s = 1        
            while (TRUE)
            {
              newZ = Z[,i] + s*dir
              newmu = mu[,i] + s*GWdir
              newQ = Z2[i]/2 + s*dirZ + s^2*dir2/2 - sum(kappa[,i]*newmu)
              newQ = newQ + crossprod(rho[-Poidx,i],(newmu[-Poidx]-psi[-Poidx,i])^2)/2
              newrho = rho[Poidx,i]*exp(s*GWdir[Poidx])
              newQ = newQ + sum(newrho)
              if ( Q[i]-newQ >= -s*dg/2 )
              {
                Z[,i] = newZ
                mu[,i] = newmu
                rho[Poidx,i] = newrho
                break
              }
              s = s / 2
            }
          }
        }
        else
        {
          for ( i in 1:n )
          {
            while (TRUE)
            {
              tmp = Z[,i] - s_Z[i]*grad[,i]
              newZ = sign(tmp)*pmax(abs(tmp)-s_Z[i]*xi,0)
              Gt = (Z[,i]-newZ)/s_Z[i] # newZ = Z[,i] - t*Gt
              Gt2 = crossprod(Gt,Gt)
              if ( Gt2 == 0 )
                break
              gradGt = crossprod(grad[,i],Gt)
              newmu = mu[,i] - s_Z[i]*GW%*%Gt
              newrho = rho[Poidx,i]*exp(newmu[Poidx]-mu[Poidx,i])
              newC2 = rho[,i]*(newmu-psi[,i])^2/2
              newC2[Poidx] = newrho
              newQ = sum(newZ^2)/2 + sum(newC2-kappa[,i]*newmu)
              
              if ( newQ <= Q[i] - s_Z[i]*gradGt + s_Z[i]*Gt2/2 )
              {
                Z[,i] = newZ
                mu[,i] = newmu
                rho[Poidx,i] = newrho
                s_Z[i] = s_Z[i] * 1.2
                break
              }
              s_Z[i] = s_Z[i] / 2
            }
          }
        }
        GWZ = mu-m
      }
      else # Poisson doesn't exist
      {
        C = kappa + rho*(psi-m)
        GWC = t(GW)%*%C
        H = array(t(GW)%*%(GW[,rep(1:L,n)]*rho[,rep(1:n,each=L)]),c(L,L,n)) + rep(diag(1,L),n)
        cH = array(apply(H,3,chol),c(L,L,n))
        SigmaZ = array(apply(cH,3,chol2inv),c(L,L,n))
        for ( i in 1:n )
          if ( xi == 0 )
            Z[,i] = SigmaZ[,,i] %*% GWC[,i]
          else
          {
            dwl = dwl.new(XX=H[,,i],Xy=GWC[,i])
            dwl.fit(dwl,rep(xi,L))
            Z[,i] = dwl$coef
            rm(dwl)
          }
        GWZ = GW%*%Z
        mu = m + GWZ
        gWSigmaZ = array(GW%*%matrix(SigmaZ,L),c(p,L,n))
        g2WWSigmaZ = apply(gWSigmaZ*rep(GW,n),c(1,3),sum)
      }
  
          
      # E-step for rho
      mu_psi = mu - psi
      for ( j in 1:p )
      {
        if ( type[j] == 0 )
        {
          a_rho = (param[j]+n)/2
          b_rho = param[j]/2 + sum(mu_psi[j,]^2)/2 + sum(g2WWSigmaZ[j,])/2
          rho[j,] = a_rho/b_rho
        }
        else if ( type[j] != 3 )
        {
          varphi = sqrt(mu_psi[j,]^2 + g2WWSigmaZ[j,])
          rho[j,] = b[j,]*(1-varphi/2+varphi^2/6)/(4-2*varphi+varphi^2)
          vphidx = which(varphi > 1e-5)
          rho[j,vphidx] = b[j,vphidx]*(1-exp(-varphi[vphidx]))/2/varphi[vphidx]/(1+exp(-varphi[vphidx]))
        }
      }
  
          
      # M-step for W
      lam = lam0 + theta*dlam
      Q = apply(rho-kappa*mu,1,sum)
      grad = ((rho-kappa)%*%t(Z) + apply(array(matrix(gWSigmaZ,p)*rho[,rep(1:n,each=L)],c(p,L,n)),c(1,2),sum))*g
      if ( fit.m )
      {
        C = kappa + rho*psi
        dm = apply(rho-kappa,1,sum)
        A = matrix(0,L+1,L+1)
        sumC = apply(C,1,sum)
        for ( j in 1:p )
        {
          if ( type[j] != 3 )
          {
            A[1,1] = sum(rho[j,])
            A[-1,1] = g[j]*Z%*%rho[j,]
            A[1,-1] = A[-1,1]
            A[-1,-1] = g[j]^2*(apply(SigmaZ*rep(rho[j,],each=L^2),c(1,2),sum) + Z%*%(t(Z)*rho[j,]))
            B = c(sumC[j],g[j]*Z%*%C[j,])
            dwl = dwl.new(XX=A,Xy=B)
            dwl.fit(dwl,c(1e-10,lam[j,]))
            m[j] = dwl$coef[1]
            W[j,] = dwl$coef[-1]
            rm(dwl)
          }
          else
          {
            while (TRUE) 
            {
              newm = m[j] - s_W[j]*dm[j]
              tmp = W[j,] - s_W[j]*grad[j,]
              newW = sign(tmp)*pmax(abs(tmp)-s_W[j]*lam[j,],0)
              Gt = c(dm[j],(W[j,]-newW)/s_W[j]) # newW = W[j,] - t*Gt
              Gt2 = crossprod(Gt,Gt)
              if ( Gt2 == 0 )
                break
              gradGt = crossprod(c(dm[j],grad[j,]),Gt)
              gGtSigmaZ = g[j]*matrix(t(Gt[-1])%*%matrix(SigmaZ,L),L)
              g2WSigmaZGt = g[j]*t(Gt[-1])%*%gWSigmaZ[j,,]
              g2Gt2SigmaZ = g[j]*t(Gt[-1])%*%gGtSigmaZ
              newmu = mu[j,] - s_W[j]*(dm[j]+g[j]*t(Gt[-1])%*%Z)
              newg2WWSigmaZ = g2WWSigmaZ[j,] - 2*s_W[j]*g2WSigmaZGt + s_W[j]^2*g2Gt2SigmaZ
              newrho = exp(newmu+newg2WWSigmaZ/2)
              newQ = sum(newrho - kappa[j,]*newmu)
              if ( newQ <= Q[j] - s_W[j]*gradGt + s_W[j]*Gt2/2 )
              {
                m[j] = newm
                W[j,] = newW
                rho[j,] = newrho
                s_W[j] = s_W[j] * 1.2
                break
              }
              s_W[j] = s_W[j] / 2
            }
          }
        }
      }
      else
      {
        C = kappa + rho*(psi-m)
        for ( j in 1:p )
        {
          if ( type[j] != 3 )
          {
            A = g[j]^2*(apply(SigmaZ*rep(rho[j,],each=L^2),c(1,2),sum) + Z%*%(t(Z)*rho[j,]))
            B = g[j]*Z%*%C[j,]
            dwl = dwl.new(XX=A,Xy=B)
            dwl.fit(dwl,lam[j,])
            W[j,] = dwl$coef
            rm(dwl)
          }
          else
          {
            while (TRUE)
            {
              tmp = W[j,] - s_W[j]*grad[j,]
              newW = sign(tmp)*pmax(abs(tmp)-s_W[j]*lam[j,],0)
              Gt = (W[j,]-newW)/s_W[j] # newW = W[j,] - t*Gt
              Gt2 = crossprod(Gt,Gt)
              if ( Gt2 == 0 )
                break
              gradGt = crossprod(grad[j,],Gt)
              gGtSigmaZ = g[j]*matrix(t(Gt)%*%matrix(SigmaZ,L),L)
              g2WSigmaZGt = g[j]*t(Gt)%*%gWSigmaZ[j,,]
              g2Gt2SigmaZ = g[j]*t(Gt)%*%gGtSigmaZ
              newmu = mu[j,] - s_W[j]*g[j]*t(Gt)%*%Z
              newg2WWSigmaZ = g2WWSigmaZ[j,] - 2*s_W[j]*g2WSigmaZGt + s_W[j]^2*g2Gt2SigmaZ
              newrho = exp(newmu+newg2WWSigmaZ/2)
              newQ = sum(newrho - kappa[j,]*newmu)
              if ( newQ <= Q[j] - s_W[j]*gradGt + s_W[j]*Gt2/2 )
              {
                W[j,] = newW
                rho[j,] = newrho
                s_W[j] = s_W[j] * 1.2
                break
              }
              s_W[j] = s_W[j] / 2
            }
          }
        }
      }
      

      if ( sstype == "Gaussian" )
        df = sum(theta>0.5)
      else
        df = sum(W!=0)
      
      GW = W*g
      GWZ = GW%*%Z
      mu = m + GWZ
      
      gWSigmaZ = array(GW%*%matrix(SigmaZ,L),c(p,L,n))
      g2WWSigmaZ = apply(gWSigmaZ*rep(GW,n),c(1,3),sum)
      
      if ( max(abs(pW-W)) < eps )
        break
      
      
      # M-step for A and Rotation
      if ( PXL )
      {
        A = apply(SigmaZ,c(1,2),mean) + Z%*%t(Z)/n
        A_L = t(chol(A))
        W = W %*% A_L
        GW = GW %*% A_L
        Z = forwardsolve(A_L,Z)
        SigmaZ = array(forwardsolve(A_L,matrix(aperm(array(forwardsolve(A_L,matrix(SigmaZ,L)),c(L,L,n)),c(2,1,3)),L)),c(L,L,n))
        gWSigmaZ = array(GW%*%matrix(SigmaZ,L),c(p,L,n))
      }
    }
    
    Ws[,,etai] = W
    Zs[,,etai] = Z
    rhos[,,etai] = rho
    thetas[,,etai] = theta
    LLKs[etai] = llk(X,mu,type,param)
    BICs[etai] = -2*LLKs[etai] + log(n)*df
    iters[etai] = iter
  }
  
  print(s_Z)
  print(s_W)
  list(p=p,n=n,L=L,e=e,v0=v0,v1=v1,delta=delta,eta=eta,param=param,m=m,g=g,W=Ws,Z=Zs,rho=rhos,theta=thetas,LLK=LLKs,BIC=BICs,iter=iters)
}


init <- function(X,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  Y = matrix(0,p,n)
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
    {
      Y[j,] = X[j,]
    }
    if ( type[j] == 1 )
    {
      pbar = (X[j,]+1)/(param[j]+2)
      Y[j,] = log(pbar/(1-pbar))
    }
    if ( type[j] == 2 )
    {
      pbar = (X[j,]+1)/(param[j]+X[j,]+2)
      Y[j,] = log(pbar/(1-pbar))
    }
    if ( type[j] == 3 )
    {
      Y[j,] = log(pmax(X[j,],1/2))
    }
  }
  Y
}


# function deviance: calculate the deviance
llk <- function(X,mu,type,param)
{
  p = nrow(X)
  n = ncol(X)
  
  l = 0
  for ( j in 1:p )
  {
    if ( type[j] == 0 )
      l = l - n*log(mean((X[j,]-mu[j,])^2))/2
    else if ( type[j] == 1 )
      l = l + sum(dbinom(X[j,],param[j],1/(1+exp(-mu[j,])),TRUE))
    else if ( type[j] == 2 )
      l = l + sum(dnbinom(X[j,],param[j],1/(1+exp(mu[j,])),log=TRUE))
    else if ( type[j] == 3 )
      l = l + sum(dpois(X[j,],exp(mu[j,]),TRUE))
  }
  l
}
