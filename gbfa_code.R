
#source("../GBFA.R")
data_type = 1
data_size = 1
data_type_list  = list()
data_type_list[[1]]  = "gaussian" 
data_type_list[[2]]  = "binary" 
data_type_list[[3]]  = "binomial" 
data_type_list[[4]]  = "mix" 
data_size_list = list()
data_size_list[[1]] = "low"
data_size_list[[2]] = "large"

suffix_graph = list()
suffix_graph[[1]] = "g_1"
suffix_graph[[2]] = "g_2"
suffix_graph[[3]] = "g_3"

ddtt_list = matrix(rep(list(),8),nrow = 4, ncol =2)
ddtt_list[1,1] = list(rep(0,90))
ddtt_list[1,2] = list(rep(0,270))
ddtt_list[2,1] = list(rep(1,90))
ddtt_list[2,2] = list(rep(1,270))
ddtt_list[3,1] = list(rep(1,90))
ddtt_list[3,2] = list(rep(1,270))
ddtt_list[4,1] = list(c(rep(0,30),rep(1,60)))
ddtt_list[4,2] = list(c(rep(0,90),rep(1,180)))


  
    #data_1 =readRDS(paste0("../save_data.rds"))
    trials = data_1$trials
    W = data_1$V
    Z = t(data_1$U)
    X = data_1$X
    mu = W %*% Z
    n = 200
    H = 3
    L = dim(W)[2]
    p = dim(W)[1]
    #grid_L  = c(L-1,L,L+1)
    data_tbs = list('true_W'= W,'true_Z'=Z,'true_X'=X)
    
    p1 = p/3
    p2 = p/3
    p3 = p/3
    data_dim = c(p1,p2,p3)
    ind_s = c(1,p1+1,p2+p1+1,p1+p2+p3)
    ind_e = c(p1,p1+p2,p1+p2+p3)
    pathway_list=list()
    pathway_list[[1]] = rep(p1/3,3)
    pathway_list[[2]] = rep(p2/3,3)
    pathway_list[[3]] = rep(p3/3,3)
    
    ###
    #graph = readRDS(paste0("../graph.rds"))
    ############# argument for GBFA_EM #####################
    
    edge_index = which(graph==1,arr.ind = TRUE)
    edge_index= edge_index[order(edge_index[,1],edge_index[,2],decreasing=FALSE),]
    
    
    ########################################################
    grid_L  = c(2,3,4)
    w_ini_l = list()
    z_ini_l = list()
    L_seed = 1113
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      w_ini_l[[l]]=matrix(rnorm(LL*p,0,1),nrow = p,ncol = LL)
    }
    
    for (l in 1:length(grid_L)) {
      set.seed(l%%3+L_seed)
      LL= grid_L[l]
      z_ini_l[[l]] = matrix(rnorm(LL*n,0,1),nrow=LL,ncol=n)
    }
    
    # tuning EM
    param = trials[,1]
    ddtt = ddtt_list[data_type,data_size][[1]]
    
    grid_v0 =c(0.1,0.15,0.2)
    grid_v1 = c(1,0.5)
    bic_em = c()
    for (l in 1:length(grid_L)) {
      LL= grid_L[l]
      for (i in 1:length(grid_v0)) {
        v0 = grid_v0[i]
        for (j in 1:length(grid_v1)) {
          v1 = grid_v1[j]
          data_GBFAEM = GBFA_EM(X,ddtt,param=param,edge_index,L=LL,v0=v0,v1=v1,fit.m=TRUE,m.init=0,W.init=w_ini_l[[l]],delta=1,eta=0.5)
          names = paste0('GBFAEM_',l,'_',i,'_',j)
          bic_em =c(bic_em, data_GBFAEM$BIC[2])
          data_tbs[[names]] = data_GBFAEM
        }
      }
    }
    




