cond <- grepl('/hpc', getwd())
if (cond){
  source('/hpc/group/adrc/dcg27/african_american_multiome/scripts/config.R')
} else {
  source('C:/Users/danie/Desktop/african_american_multiome/scripts/config.R')
} ; rm(cond)
set.seed(1)

library(glmnet)
library(MatrixGenerics)
library(numDeriv)

slurm_array_task_id = Sys.getenv('SLURM_ARRAY_TASK_ID') %>% as.numeric()

setwd(apoe_ccre_objects)
gene_pairs <- list.files(pattern = '3.prep_lasso___aa_modules_.*.rds')[slurm_array_task_id]
chr <- gsub('.*_|.rds', '', gene_pairs)
gene_pairs <- readRDS(gene_pairs)
atac <- paste0('3.prep_lasso___aa_atac_metacells_', chr, '.rds')
atac <- readRDS(atac)
rna <- paste0('3.prep_lasso___aa_rna_metacells_', chr, '.rds')
rna <- readRDS(rna)

# make metacells into dataframe
atac <- t(atac)
ix <- colSums(atac) > 0
atac <- atac[, ix]
atac <- as.data.frame(atac)

rna <- t(rna)
ix <- colSums(rna) > 0
rna <- rna[, ix]
rna <- as.data.frame(rna)

# soft threshold operator
S = function(x, lambda, rho){
  thresh = lambda / rho
  sign(x) * pmax(abs(x) - thresh, 0)
}

# augmented lagrangian
lossFunction <- function(X, Y, B, D, z, d, u, w, L1, L2, r1, r2, DB){
  N = nrow(X)
  term1 = 1 / (2 * N) * sum((Y - X %*% B)^2)
  term2 = L1 * sum(abs(z)) + L2 * sum(abs(d)) 
  term3 = r1/2 * sum((B - z + u)^2) 
  term4 = r2/2 * sum((DB - d + w)^2)
  return(term1 + term2 + term3 + term4)
}

# fused lasso ADMM algorithm
fusedLASSO = function(X, Y, B, D, z, d, u, w, L1, L2, r1, r2,
                      XtY, XtX, DtD, Dt, verbose = FALSE){
  best_loss = Inf
  

  for (iter in 1:1000){
    
    # message(iter)
    if (iter == 1){
      counter = 0
    }
    
    # update B
    I = diag(1, nrow = length(B))
    if (iter == 1){
      A = 1/N * XtX + r1 * I + r2 * DtD
      R = chol(A)
    }
    b = 1/N * XtY + r1 * (z - u) + r2 * Dt %*% (d - w)
    # B = solve(A, b)
    B = backsolve(R, forwardsolve(t(R), b))
    DB = D %*% B
    
    # update z and d
    z = S(B + u, lambda = L1, rho = r1)
    d = S(DB + w, lambda = L2, rho = r2)
    
    # update u and w
    u = u + (B - z)
    w = w + (DB - d)
    
    loss = lossFunction(X, Y, B, D, z, d, u, w, L1, L2, r1, r2, DB)
    
    if (iter %% 10 == 0 & verbose){
      message('LOSS: ', round(loss, 5))
    }
    # stop if converge
    if (best_loss - loss < 1e-5){
      counter = counter + 1
    } else {
      counter = 0
      best_loss = loss
    }
    if (counter > 10){
      if (verbose){
        message('convergence reached')
      }
      break
    }
  }
  results = list(B = B, z = z, d = d, u = u, w = w)
  return(results)
}

genes <- unique(gene_pairs$id)
# genes <- sample(genes, 100)
master <- NULL
i = 1
j = 1
k = 1
for (i in 1:length(genes)){
  try(expr = {
    message('running model ', i, ' of ', length(genes))
    gene_pairs_i <- gene_pairs[gene_pairs$id == genes[i], ]
    df <- atac[, c('dgx', 'geno', gene_pairs_i$peak)]
    message(genes[i], ': running random lasso')
    
    # group specific terms
    X <- df[, -c(1,2)] %>% as.matrix()
    X <- scale(X)
    dgx = df$dgx == 1
    Z <- X * dgx
    X <- X * !dgx
    colnames(Z) <- paste0('dgx_', colnames(X))
    
    # intercept term
    W = df[, c('dgx', 'geno')]
    W = cbind(1, W)
    colnames(W)[1] = 'intercept'
    
    # data matrix 
    X = cbind(W, X, Z)
    
    # response variable
    target_gene <- unique(gene_pairs_i$gene)
    Y <- rna[, target_gene]
    Y <- scale(Y) %>% as.numeric()
    
    # initial lambda estimates
    penalty_factor = rep(1, ncol(X))
    penalty_factor[1:3] <- 0
    adaptive_lasso <- glmnet(x = X, y = Y, alpha = 1, nlambda = 20,
                             standardize = FALSE, penalty.factor = penalty_factor, 
                             intercept = FALSE) 
    lambda = adaptive_lasso$lambda
    
    load_cells = which(X[, 'dgx'] == 1)
    normal_cells = which(X[, 'dgx'] == 0)
    
    message('running random fused lasso stability selection')
    Bstore = Dstore = 0
    Bstability = Dstability = matrix(0, nrow = ncol(X)-3, ncol = 5)
    rownames(Bstability) = rownames(Dstability) = paste0(genes[i], '_', colnames(X)[4:ncol(X)])
    colnames(Bstability) = colnames(Dstability) = paste0('lambda', 1:ncol(Bstability))
    for (j in 1:20){
      message(j, '/20')
      ix1 = sample(normal_cells, length(normal_cells)/2)
      ix2 = sample(load_cells, length(load_cells)/2)
      ix = c(ix1, ix2)
      # group specific terms
      X <- df[ix, -c(1,2)] %>% as.matrix()
      X <- scale(X)
      dgx = df$dgx[ix] == 1
      Z <- X * dgx
      X <- X * !dgx
      colnames(Z) <- paste0('dgx_', colnames(X))
      
      # data matrix 
      X = cbind(X, Z)
      
      # response variable
      target_gene <- unique(gene_pairs_i$gene)
      Y <- rna[ix, target_gene]
      Y <- scale(Y) %>% as.numeric()
      
      # center Y to remove need for intercept terms
      init = glm(Y ~ ., family = 'gaussian', data = df[ix, c('dgx', 'geno')])
      Y = resid(init, type = 'response')

      N = nrow(X)
      p = ncol(X)
      
      X = as.matrix(X)
      XtY = t(X) %*% Y
      XtX = t(X) %*% X
      D = matrix(0, nrow = ncol(X), ncol = ncol(X))
      ix1 = 1:(p/2)
      ix2 = ix1 + length(ix1)
      D[cbind(ix1, ix1)] = 1
      D[cbind(ix1, ix2)] = -1
      DtD = t(D) %*% D
      Dt = t(D)
      
      if (j == 1){
        Bcache = zcache = dcache = ucache = wcache = list()
      }
      
      for (k in 1:5){
        lambda_k = lambda[k]
        
        if (j == 1){
          B = rep(0, ncol(X))
          z = B
          d = D %*% B
          u = rep(0, length(B))
          w = rep(0, nrow(D))
        } else {
          B = Bcache[[k]]
          z = zcache[[k]]
          d = dcache[[k]]
          u = ucache[[k]]
          w = wcache[[k]]
        }
        
        L1 = lambda_k
        L2 = 1 * L1
        r1 = r2 = L1

        results = fusedLASSO(X, Y, B, D, z, d, u, w, L1, L2, r1, r2,
                       XtY, XtX, DtD, Dt, verbose = FALSE)
        B = results$z
        Bint = results$d
        # plot(B)
        
        # cache results
        if (j == 1){
          Bcache[[k]] = results$B
          zcache[[k]] = results$z
          dcache[[k]] = results$d
          ucache[[k]] = results$u
          wcache[[k]] = results$w
        }
        Btally = B != 0
        Dtally = Bint != 0
        
        Bstability[, k] = Bstability[, k] + Btally
        Dstability[, k] = Dstability[, k] + Dtally
        
        if (k == 5){
          Bstore = Bstore + B
          Dstore = Dstore + Bint
        }
      }
    }
    
    Bstore = Bstore / j
    Dstore = Dstore / j 
    eff_sizes = cbind(Bstore, Dstore)
    rownames(eff_sizes) = rownames(Bstability)
    colnames(eff_sizes) = c('group_level_eff', 'interaction_eff')
    
    Bstability = Bstability / j
    Dstability = Dstability / j
    Bprobs = apply(Bstability, 1, max)
    Dprobs = apply(Dstability, 1, max)
    
    selection_probs = cbind(Bprobs, Dprobs)
    colnames(selection_probs) <- c('group_specific_prob', 'interaction_prob')

    if (i == 1){
      EFF_SIZES = eff_sizes
      SELECTION_PROBS = selection_probs
    } else {
      EFF_SIZES = rbind(EFF_SIZES, eff_sizes)
      SELECTION_PROBS = rbind(SELECTION_PROBS, selection_probs)
    }
    
    RESULTS = cbind(SELECTION_PROBS, EFF_SIZES)
    
    if (i %% 100 == 0){
      setwd(apoe_ccre_objects)
      saveRDS(RESULTS, paste0('4.fused_lasso___', chr, '.rds'))
    }
  })
}

setwd(apoe_ccre_objects)
saveRDS(RESULTS, paste0('4.fused_lasso___', chr, '.rds'))

