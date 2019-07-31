divideData <- function(dt_merged, desfecho){
  browser()
  # o desfecho deve ser numerico e dicotomico 0 ou 1
  require(purrr)
  
  y <- dt_merged[, desfecho]
  dt_merged <- dt_merged[, -which(colnames(dt_merged) == desfecho)]
  dt_merged <- dt_merged[!is.na(y), ]
  y <- y[!is.na(y)]
  
  w <- y == 0
  data1 <- dt_merged[w, ]
  data2 <- dt_merged[!w, ]
  
  nnzv1 <- map_lgl(data1, function(x){
    cond1 <- is.factor(x)
    cond2 <- min(table(x)) > 8
    
    if(cond1 & cond2) y <- TRUE else y <- FALSE
    y
  }
  )
  
  nnzv2 <- map_lgl(data2, function(x){
    cond1 <- is.factor(x)
    cond2 <- min(table(x)) > 8
    
    if(cond1 & cond2) y <- TRUE else y <- FALSE
    y
  }
  )
  
  nnzv1 <- names(nnzv1[nnzv1])
  nnzv2 <- names(nnzv2[nnzv2])
  
  return(list("data1" = data1, "data2" = data2, 
              "valid_vars_data1" = nnzv1, "valid_vars_data2" = nnzv2))
}

getVarsType <- function(dt){
  require(purrr)
  vars_type <- map_chr(dt, class)
  vars_type <- as.factor(vars_type)
  levels(vars_type) <- c("c", "g")
  vars_type <- as.character(vars_type)
  vars_type
}

getVarsLevels <- function(dt){
  require(purrr)
  vars_levels <- map_int(dt, function(x){length(table(x))})
  vars_levels
}

dealMissing <- function(dt){
  train_matrix <- dt
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  for (i in names(train_matrix)){
    if (is.factor(train_matrix[, i])) {
      print(i)
      mode_value = getmode(train_matrix[,i])
      #recorded_values[1,i] = mode_value
      train_matrix[is.na(train_matrix[,i]),i] = mode_value
      
    } else {
      mean_value = mean(train_matrix[,i], na.rm=TRUE)
      #recorded_values[1,i] = mean_value
      train_matrix[is.na(train_matrix[,i]),i] = mean_value
      
    }
    
  }
  return(train_matrix)
}

changeValue <- function(var, current_value, new_value){
  w <- which(var == current_value)
  var[w] <- new_value
  var
}

calculateNetworkMeasures <- function(data1, data2, n_iterations, vars_type, vars_levels){
  browser()
  require(mgm)
  require(qgraph)
  ##### Invariance measures #####
  x1 <- data1
  x2 <- data2
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  data.list <- list(x1,x2)
  b <- 1:(nobs1+nobs2)
  nvars <- ncol(x1)
  nedges <- nvars*(nvars-1)/2
  
  it <- n_iterations
  
  glstrinv.perm <- glstrinv.real <- nwinv.real <- nwinv.perm <- c()
  diffedges.perm <- matrix(0,it,nedges) 
  einv.perm.all <- array(NA,dim=c(nvars, nvars, it))
  corrpvals.all <- matrix(NA,nvars,nvars)
  edges.pvalmattemp <- matrix(0,nvars,nvars)
  
  
  fit_mgm1 <- mgm(data = x1, type = vars_type, 
                  levels = vars_levels, k = 2, lambdaSel = "CV", 
                  lambdaFolds = 10, ruleReg = "AND")
  
  fit_mgm2 <- mgm(data = x2, type = vars_type, 
                  levels = vars_levels, k = 2, lambdaSel = "CV", 
                  lambdaFolds = 10, ruleReg = "AND")
  
  Q1 <- qgraph(fit_mgm1$pairwise$wadj, edge.color = fit_mgm1$pairwise$edgecolor,
               nodeNames = colnames(x1), legend = TRUE, legend.cex = 0.1, label.cex = 3, vsize = 3, layout = "spring")
  
  Q2 <- qgraph(fit_mgm2$pairwise$wadj, edge.color = fit_mgm2$pairwise$edgecolor,
               nodeNames = colnames(x2), legend = TRUE, legend.cex = 0.1, label.cex = 3, vsize = 3, layout = "spring")
  
  nw1 <- fit_mgm1$pairwise$wadj
  nw2 <- fit_mgm2$pairwise$wadj
  
  ## Global strength invariance
  glstrinv.real <- abs(sum(abs(nw1[upper.tri(nw1)]))-sum(abs(nw2[upper.tri(nw2)])))
  glstrinv.real
  
  # global strength of individual networks
  glstrinv.sep <- c(sum(abs(nw1[upper.tri(nw1)])), sum(abs(nw2[upper.tri(nw2)])))
  glstrinv.sep
  
  ## Individual edge invariance
  diffedges.real <- abs(nw1-nw2)[upper.tri(abs(nw1-nw2))] 
  diffedges.realmat <- matrix(diffedges.real,it,nedges,byrow=TRUE)
  diffedges.realoutput <- abs(nw1-nw2)
  
  ## Network structure invariance
  nwinv.real <- max(diffedges.real)
  
  #####################################
  #####     Start permutations    #####
  #####################################
  
  for (i in 1:it)
  {
    diffedges.permtemp <- matrix(0, nvars, nvars)
    
    s <- sample(1:(nobs1+nobs2), nobs1, replace=FALSE)
    x1perm <- dataall[s,]
    x2perm <- dataall[b[-s],]
    
    
    # If not paired data ----
    fit_perm1 <- mgm(data = x1perm, type = vars_type, 
                     levels = vars_levels, k = 2, lambdaSel = "CV", 
                     lambdaFolds = 10, ruleReg = "AND")
    
    fit_perm2 <- mgm(data = x2perm, type = vars_type, 
                     levels = vars_levels, k = 2, lambdaSel = "CV", 
                     lambdaFolds = 10, ruleReg = "AND")
    
    r1perm <- fit_perm1$pairwise$wadj
    r2perm <- fit_perm2$pairwise$wadj
    
    
    
    ## Invariance measures for permuted data ----
    glstrinv.perm[i] <- abs(sum(abs(r1perm[upper.tri(r1perm)]))-sum(abs(r2perm[upper.tri(r2perm)])))
    diffedges.perm[i,] <- abs(r1perm-r2perm)[upper.tri(abs(r1perm-r2perm))]
    diffedges.permtemp[upper.tri(diffedges.permtemp, diag=FALSE)] <- diffedges.perm[i,]
    diffedges.permtemp <- diffedges.permtemp + t(diffedges.permtemp)
    einv.perm.all[,,i] <- diffedges.permtemp
    nwinv.perm[i] <- max(diffedges.perm[i,])
    
    #if (progressbar==TRUE) setTxtProgressBar(pb, i)
  }
  
  p.adjust.methods <- "BH"
  
  edges.pvaltemp <- colSums(diffedges.perm >= diffedges.realmat)/it 
  
  
  ## If all edges should be tested
  # corrected p-values (or not if p.adjust.methods='none')
  corrpvals.all.temp <- round(p.adjust(edges.pvaltemp, method=p.adjust.methods),3)
  # matrix with corrected p values
  corrpvals.all
  corrpvals.all[upper.tri(corrpvals.all,diag=FALSE)] <- corrpvals.all.temp 
  rownames(corrpvals.all) <- colnames(corrpvals.all) <- colnames(x1)
  
  library(reshape2)
  einv.pvals <- melt(corrpvals.all, na.rm=TRUE, value.name = 'p-value')
  einv.perm <- einv.perm.all
  einv.real <- diffedges.realoutput
  
  edges.tested <- colnames(einv.perm)
  
  res <- list(glstrinv.real = glstrinv.real,
              glstrinv.sep = glstrinv.sep,
              glstrinv.pval = sum(glstrinv.perm >= glstrinv.real)/it, 
              glstrinv.perm = glstrinv.perm,
              nwinv.real = nwinv.real,
              nwinv.pval = sum(nwinv.perm >= nwinv.real)/it, 
              nwinv.perm = nwinv.perm,
              edges.tested = edges.tested,
              einv.real = einv.real,
              einv.pvals = einv.pvals,
              einv.perm = einv.perm, 
              nw1 = nw1,
              nw2 = nw2)
  
  return(list("res" = res, "fit_mgm1" = fit_mgm1, "fit_mgm2" = fit_mgm2, "net1" = Q1, "net2" = Q2))
  
}

simulatePermutations <- function(data1, data2, n_iterations, vars_type){
  x1 <- data1
  x2 <- data2
  nobs1 <- nrow(x1)
  nobs2 <- nrow(x2)
  dataall <- rbind(x1,x2)
  b <- 1:(nobs1+nobs2)
  
  it <- n_iterations
  cat_vars <- which(vars_type == "c")
  perm_matrix <- matrix(nrow = it, ncol = length(cat_vars))
  perm_matrix2 <- perm_matrix
  
  for (i in 1:it){
    s <- sample(1:(nobs1+nobs2), nobs1, replace=FALSE)
    x1perm <- dataall[s,]
    x2perm <- dataall[b[-s],]
    x1perm <- x1perm[, cat_vars]
    x2perm <- x2perm[, cat_vars]
    
    nnzv1 <- map_lgl(x1perm, function(x){cond1 <- min(table(x)) >= 8; cond1})
    nnzv2 <- map_lgl(x2perm, function(x){cond1 <- min(table(x)) >= 8; cond1})
    
    perm_matrix[i, ] <- nnzv1
    perm_matrix2[i, ] <- nnzv2
    
  }
  
  near_zero <- map_lgl(as.data.frame(perm_matrix), function(x){if (length(which(x == FALSE)) >= 1) nz <- TRUE else nz <- FALSE; nz})
  near_zero2 <- map_lgl(as.data.frame(perm_matrix2), function(x){if (length(which(x == FALSE)) >= 1) nz <- TRUE else nz <- FALSE; nz})
  
  names(near_zero) <- colnames(data1)[cat_vars]
  near_zero <- names(near_zero)[near_zero == TRUE]
  
  names(near_zero2) <- colnames(data1)[cat_vars]
  near_zero2 <- names(near_zero2)[near_zero2 == TRUE]
  
  near_zero_all <- union(near_zero, near_zero2)
  
  return(list("perm_matrix" = perm_matrix, "perm_matrix2" = perm_matrix2, "near_zero" = near_zero_all))
}


