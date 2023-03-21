#-------------------------------------------------#
#            generate simulation data             #
#-------------------------------------------------#
SimuData <- function(NumSeed, NumObs, NumGene, NumDEGene, MinCount, MaxCount, herit, sigma2, PVE, MiuCons, RefK, Setting){
  # PVE dose not affect null settings
  
  set.seed(NumSeed)
  
  # K simulation
  DimK <- NumObs
  DimRefK <- dim(RefK)[1]
  SimuK <- diag(0, nrow=DimK, ncol=DimK)
  # fill its off-diagonal elements with the random drawn off-diagonal elements of RefK
  for(row_ind in 1 : (DimK-1)){
    for(col_ind in (row_ind+1) : DimK){
      rand_i <- sample(1:DimRefK, 1, replace=T)
      rand_j <- sample(setdiff(1:DimRefK, rand_i), 1, replace=T)
      SimuK[row_ind, col_ind] <- SimuK[col_ind, row_ind] <- RefK[rand_i, rand_j]
    }
  }
  diag(SimuK) <- apply(SimuK, 1, sum) * -1
  #scale_c <- NumObs / sum(diag(SimuK))
  #SimuK <- scale_c * SimuK
  # make K a positive definite matrix if it is not
  SimuK_eign <- eigen(SimuK)$values
  if (min(SimuK_eign) < 0){
    ApproK <- Matrix::nearPD(SimuK, ensureSymmetry=T)$mat@x
    SimuK <- matrix(ApproK, nrow=DimK, byrow=F)
  }
  SimuK <- cov2cor(SimuK)
  
  # predictor x simulation
  pred <- rnorm(NumObs, mean=0, sd=1)
  pred <- (pred - mean(pred)) / sd(pred)
  
  # total count simulation
  tcount <- sample(MinCount:MaxCount, NumObs, replace=T)
  tcount <- as.data.frame(t(tcount))
  colnames(tcount) <- paste("idv", 1:NumObs, sep="")
  
  # gene count simulation
  DE_ind <- sample(1:NumGene, NumDEGene, replace=F)
  gcount_null <- gcount_alt <- matrix(NA, nrow=NumGene, ncol=NumObs)
  isDE_null <- isDE_alt <- rep("unknown", NumGene)
  beta <- rep(0, NumGene)
  
  for (j in 1:NumGene){
    RanEffect <- MASS::mvrnorm(1, mu=rep(0, NumObs), Sigma=SimuK)
    EnvEffect <- rep(0, NumObs)
    for (i in 1:NumObs){
      EnvEffect[i] <- rnorm(1, mean=0, sd=1)
    }
    # scale random effect and environment effect
    var_g <- herit * sigma2
    var_e <- sigma2 - var_g
    ScaFactor_g <- sqrt(var_g / var(RanEffect))
    RanEffect <- RanEffect * ScaFactor_g
    ScaFactor_e <- sqrt(var_e / var(EnvEffect))
    EnvEffect <- EnvEffect * ScaFactor_e
    #var(RanEffect) / (var(RanEffect) + var(EnvEffect))
    
    # compute miu
    miu <- log(MiuCons / mean(as.numeric(tcount)))
    
    # under the null setting
    if (Setting == 'null'){
      lambda_null <- exp(miu + RanEffect + EnvEffect)
      for (i in 1:NumObs){
        N_lambda_i <- lambda_null[i] * tcount[1, i]
        gcount_null[j, i] <- rpois(1, N_lambda_i)
      }
      isDE_null[j] <- "non-DE"
    }
    
    # under the alternative setting
    if (Setting == 'alt' & j %in% DE_ind){
      sd_beta_j <- sqrt((PVE*(var(RanEffect) + var(EnvEffect))) / (var(pred)*(1-PVE)))
      beta_j <- rnorm(n=1, mean=0, sd=sd_beta_j)
      lambda_alt <- exp(miu + pred*beta_j + RanEffect + EnvEffect)
      isDE_alt[j] <- "DE"
      beta[j] <- beta_j
      for (i in 1:NumObs){
        N_lambda_i <- lambda_alt[i] * tcount[1, i]
        gcount_alt[j, i] <- rpois(1, N_lambda_i)
      }
    }
    
    if (Setting == 'alt' & (!(j %in% DE_ind))){
      isDE_alt[j] <- "non-DE"
      beta[j] <- 0
      lambda_null <- exp(miu + RanEffect + EnvEffect)
      for (i in 1:NumObs){
        N_lambda_i <- lambda_null[i] * tcount[1, i]
        gcount_alt[j, i] <- rpois(1, N_lambda_i)
      }
    }
  }
  
  
  if (Setting == 'null'){
    gcount_null <- data.frame(paste("gene", 1:NumGene, sep=""), gcount_null, isDE_null, stringsAsFactors=FALSE)
    colnames(gcount_null) <- c("Gene", paste("idv", 1:NumObs, sep=""), "isDE")
    return(list(GeneCountNull=gcount_null, TotalCount=tcount, predictor=pred, KMat=SimuK))
  }
    
  if (Setting == 'alt'){
    gcount_alt <- data.frame(paste("gene", 1:NumGene, sep=""), gcount_alt, beta, isDE_alt, stringsAsFactors=FALSE)
    colnames(gcount_alt) <- c("Gene", paste("idv", 1:NumObs, sep=""), "beta", "isDE")
    return(list(GeneCountAlt=gcount_alt, TotalCount=tcount, predictor=pred, KMat=SimuK))
  }
}

#RefK <- read.table('c:/users/83944/Desktop/MALAX-Poi-data/BaboonData/Baboon_kernel.txt', sep='\t', header=F)
#DataTest <- SimuData(NumSeed = 10086, NumObs = 20, NumGene = 100, NumDEGene = 10, MinCount = 1770083, MaxCount = 9675989,
#                  herit = 0.3, sigma2 = 0.5, PVE = 0.35, MiuCons = 5, RefK = RefK, Setting='alt')


#-------------------------------------------------#
#         generate simulation datasets            # 
#    based on the specified parameters grid       #
#-------------------------------------------------#
RunSimuData <- function(NumObs, herit, sigma2, PVE, MiuCons, repli, NumGene, NumDEGene, MinCount, MaxCount, RefK, OutDire, NumCore, Setting){
  
  SimuData <- function(NumSeed, NumObs, NumGene, NumDEGene, MinCount, MaxCount, herit, sigma2, PVE, MiuCons, RefK, Setting){
    # PVE dose not affect null settings
    
    set.seed(NumSeed)
    
    # K simulation
    DimK <- NumObs
    DimRefK <- dim(RefK)[1]
    SimuK <- diag(0, nrow=DimK, ncol=DimK)
    # fill its off-diagonal elements with the random drawn off-diagonal elements of RefK
    for(row_ind in 1 : (DimK-1)){
      for(col_ind in (row_ind+1) : DimK){
        rand_i <- sample(1:DimRefK, 1, replace=T)
        rand_j <- sample(setdiff(1:DimRefK, rand_i), 1, replace=T)
        SimuK[row_ind, col_ind] <- SimuK[col_ind, row_ind] <- RefK[rand_i, rand_j]
      }
    }
    diag(SimuK) <- apply(SimuK, 1, sum) * -1
    #scale_c <- NumObs / sum(diag(SimuK))
    #SimuK <- scale_c * SimuK
    # make K a positive definite matrix if it is not
    SimuK_eign <- eigen(SimuK)$values
    if (min(SimuK_eign) < 0){
      ApproK <- Matrix::nearPD(SimuK, ensureSymmetry=T)$mat@x
      SimuK <- matrix(ApproK, nrow=DimK, byrow=F)
    }
    SimuK <- cov2cor(SimuK)
    
    # predictor x simulation
    pred <- rnorm(NumObs, mean=0, sd=1)
    pred <- (pred - mean(pred)) / sd(pred)
    
    # total count simulation
    tcount <- sample(MinCount:MaxCount, NumObs, replace=T)
    tcount <- as.data.frame(t(tcount))
    colnames(tcount) <- paste("idv", 1:NumObs, sep="")
    
    # gene count simulation
    DE_ind <- sample(1:NumGene, NumDEGene, replace=F)
    gcount_null <- gcount_alt <- matrix(NA, nrow=NumGene, ncol=NumObs)
    isDE_null <- isDE_alt <- rep("unknown", NumGene)
    beta <- rep(0, NumGene)
    
    for (j in 1:NumGene){
      RanEffect <- MASS::mvrnorm(1, mu=rep(0, NumObs), Sigma=SimuK)
      EnvEffect <- rep(0, NumObs)
      for (i in 1:NumObs){
        EnvEffect[i] <- rnorm(1, mean=0, sd=1)
      }
      # scale random effect and environment effect
      var_g <- herit * sigma2
      var_e <- sigma2 - var_g
      ScaFactor_g <- sqrt(var_g / var(RanEffect))
      RanEffect <- RanEffect * ScaFactor_g
      ScaFactor_e <- sqrt(var_e / var(EnvEffect))
      EnvEffect <- EnvEffect * ScaFactor_e
      #var(RanEffect) / (var(RanEffect) + var(EnvEffect))
      
      # compute miu
      miu <- log(MiuCons / mean(as.numeric(tcount)))
      
      # under the null setting
      if (Setting == 'null'){
        lambda_null <- exp(miu + RanEffect + EnvEffect)
        for (i in 1:NumObs){
          N_lambda_i <- lambda_null[i] * tcount[1, i]
          gcount_null[j, i] <- rpois(1, N_lambda_i)
        }
        isDE_null[j] <- "non-DE"
      }
      
      # under the alternative setting
      if (Setting == 'alt' & j %in% DE_ind){
        sd_beta_j <- sqrt((PVE*(var(RanEffect) + var(EnvEffect))) / (var(pred)*(1-PVE)))
        beta_j <- rnorm(n=1, mean=0, sd=sd_beta_j)
        lambda_alt <- exp(miu + pred*beta_j + RanEffect + EnvEffect)
        isDE_alt[j] <- "DE"
        beta[j] <- beta_j
        for (i in 1:NumObs){
          N_lambda_i <- lambda_alt[i] * tcount[1, i]
          gcount_alt[j, i] <- rpois(1, N_lambda_i)
        }
      }
      
      if (Setting == 'alt' & (!(j %in% DE_ind))){
        isDE_alt[j] <- "non-DE"
        beta[j] <- 0
        lambda_null <- exp(miu + RanEffect + EnvEffect)
        for (i in 1:NumObs){
          N_lambda_i <- lambda_null[i] * tcount[1, i]
          gcount_alt[j, i] <- rpois(1, N_lambda_i)
        }
      }
    }
    
    
    if (Setting == 'null'){
      gcount_null <- data.frame(paste("gene", 1:NumGene, sep=""), gcount_null, isDE_null, stringsAsFactors=FALSE)
      colnames(gcount_null) <- c("Gene", paste("idv", 1:NumObs, sep=""), "isDE")
      return(list(GeneCountNull=gcount_null, TotalCount=tcount, predictor=pred, KMat=SimuK))
    }
    
    if (Setting == 'alt'){
      gcount_alt <- data.frame(paste("gene", 1:NumGene, sep=""), gcount_alt, beta, isDE_alt, stringsAsFactors=FALSE)
      colnames(gcount_alt) <- c("Gene", paste("idv", 1:NumObs, sep=""), "beta", "isDE")
      return(list(GeneCountAlt=gcount_alt, TotalCount=tcount, predictor=pred, KMat=SimuK))
    }
  }
  library(doParallel)
  
  params_sapce <- expand.grid(seed=repli[1]:repli[2], n=NumObs, h2=herit, sigma2=sigma2, miu=MiuCons, PVE=PVE)
  #write.table(params_sapce, paste(OutDire, '/AAAParamsGrid.txt', sep=''), row.names=F, col.names=T, sep='\t')
  
  cat('## Total scenarios:', nrow(params_sapce)/(repli[2] - repli[1] + 1), '\n')
  cat('## Replicates:', (repli[2] - repli[1] + 1), '\n')
  cat('## Generate simulation datasets:', '\n')
  
  cl <- parallel::makeCluster(NumCore)
  doSNOW::registerDoSNOW(cl)
  
  # do parallel to generate datasets
  pro_bar <- txtProgressBar(max=nrow(params_sapce), style=3, char='#')
  progress <- function(n) setTxtProgressBar(pro_bar, n)
  opts <- list(progress=progress)
  
  simu_data <- foreach::foreach(i=1:nrow(params_sapce), .options.snow=opts) %dopar% {
    SimuData(NumSeed=params_sapce$seed[i], NumObs=params_sapce$n[i], NumGene=NumGene, NumDEGene=NumDEGene,
             MinCount=MinCount, MaxCount=MaxCount, herit=params_sapce$h2[i], sigma2=params_sapce$sigma2[i],
             PVE=params_sapce$PVE[i], MiuCons=params_sapce$miu[i], RefK=RefK, Setting=Setting)
  }
  
  close(pro_bar)
  
  # do parallel to output datasets
  cat('## Save simulation datasets:\n')
  pro_bar <- txtProgressBar(max=nrow(params_sapce), style=3, char='#')
  progress <- function(n) setTxtProgressBar(pro_bar, n)
  opts <- list(progress=progress)
  
  file_index <- expand.grid(seed=repli[1]:repli[2], scenarios=1:(nrow(params_sapce) / (repli[2] - repli[1] + 1)))
  file_name <- paste('scenas', file_index$scenarios, file_index$seed, sep='_')
  write.table(cbind(file_name, params_sapce), paste(OutDire, '/AAAParamsGrid.txt', sep=''), row.names=F, col.names=T, sep='\t')
  
  out <- foreach::foreach(i=1:nrow(params_sapce), .options.snow=opts) %dopar% {
    dataset_i <- simu_data[[i]]
    if (Setting == 'null'){
      write.table(dataset_i$GeneCountNull[, -ncol(dataset_i$GeneCountNull)], paste(OutDire, '/gcount_null_', file_name[i], '.txt', sep=''), sep='\t',
                  row.names=F, col.names=T)
    }
    
    if (Setting == 'alt'){
      write.table(dataset_i$GeneCountAlt[, -c(ncol(dataset_i$GeneCountAlt)-1, ncol(dataset_i$GeneCountAlt))], paste(OutDire, '/gcount_alt_', file_name[i], '.txt', sep=''), sep='\t',
                  row.names=F, col.names=T)
      write.table(dataset_i$GeneCountAlt[, c('Gene', 'beta', 'isDE')], paste(OutDire, '/ginfo_alt_', file_name[i], '.txt', sep=''), sep='\t',
                  row.names=F, col.names=T)
    }

    write.table(dataset_i$TotalCount, paste(OutDire, '/tcount_', file_name[i], '.txt', sep=''), sep='\t',
                row.names=F, col.names=T)
    write.table(cbind('tcount', dataset_i$TotalCount), paste(OutDire, '/tcountmacau_', file_name[i], '.txt', sep=''), sep='\t',
                row.names=F, col.names=T)
    write.table(dataset_i$predictor, paste(OutDire, '/pred_', file_name[i], '.txt', sep=''), sep='\t',
                row.names=F, col.names=F)
    write.table(dataset_i$KMat, paste(OutDire, '/kernel_', file_name[i], '.txt', sep=''), sep='\t',
                row.names=F, col.names=F)
  }
  
  close(pro_bar)
  parallel::stopCluster(cl)
  
}

# RefK <- read.table('c:/users/83944/Desktop/MALAX-Poi-data/BaboonData/Baboon_kernel.txt', sep='\t', header=F)
# RunSimuData(NumObs = c(50, 100, 150), herit = c(0.1, 0.3, 0.6), sigma2 = c(0.1, 0.3, 0.5), PVE = c(0.15, 0.25, 0.35),
#             MiuCons = c(5, 50), repli = c(1, 1), NumGene = 100, NumDEGene = 10, MinCount = 1770083, MaxCount = 9675989,
#             RefK = RefK, OutDire = 'simudata', NumCore = 8, Setting='null')


#-------------------------------------------------#
#       run PQLseq for all simulation datasets    #
#-------------------------------------------------#
run_PQL <- function(DireData, DireOut, SimuRange, setting, NumCore, NumCorePQL=1){
  
  library(doParallel)
  DiffTime <- function(t0, unit='mins'){
    return(as.numeric(difftime(Sys.time(), t0, units=unit)))
  }
  
  NumSimu <- SimuRange[1]:SimuRange[2]
  
  cl <- parallel::makeCluster(NumCore)
  doSNOW::registerDoSNOW(cl)
  
  # do parallel to run PQLseq for each dataset
  pro_bar <- txtProgressBar(max=length(NumSimu), style=3, char='#')
  progress <- function(n) setTxtProgressBar(pro_bar, n)
  opts <- list(progress=progress)
  
  file_name <- read.table(paste(DireData, '/AAAParamsGrid.txt', sep=''), header=T)[, 1]
  res_pql <- foreach::foreach(i=NumSimu, .options.snow=opts) %dopar% {
    # read gene expression count
    if (setting == 'null'){
      gcount <- read.table(paste(DireData, '/gcount_null_', file_name[i], '.txt', sep=''), header=T, row.names=1, sep='\t', check.names=F)
    }
    if (setting == 'alt'){
      gcount <- read.table(paste(DireData, '/gcount_alt_', file_name[i], '.txt', sep=''), header=T, row.names=1, sep='\t', check.names=F)
    }
    
    # read total count
    tcount <- read.table(paste(DireData, '/tcount_', file_name[i], '.txt', sep=''), header=T, sep='\t', check.names=F)
    
    # read predictor
    pred <- read.table(paste(DireData, '/pred_', file_name[i], '.txt', sep=''), header=F, sep='\t', check.names=F)
    
    # read kernel
    kernel <- read.table(paste(DireData, '/kernel_', file_name[i], '.txt', sep=''), header=F, sep='\t')
    
    t0 <- Sys.time()
    res <- list('PQL'=(PQLseq::pqlseq(RawCountDataSet=gcount, LibSize=tcount, Phenotypes=pred, RelatednessMatrix=kernel, numCore=NumCorePQL, fit.model='PMM')),
                'time'=DiffTime(t0))
  }
  
  close(pro_bar)
  parallel::stopCluster(cl)
  
  ind <- 1
  simu_time <- simu_info <- rep(0, length(NumSimu))
  for (i in NumSimu){
    result_i <- res_pql[[ind]]$PQL
    simu_time[ind] <- res_pql[[ind]]$time
    simu_info[ind] <- file_name[i]
    file_name_i <- paste(DireOut, '/ResPQL_', file_name[i], '.txt', sep='')
    write.table(result_i, file_name_i, sep='\t', row.names=T, col.names=T)
    ind <- ind + 1
  }
  
  PQL_time <- data.frame('scenarios'=simu_info, 'time-mins'=simu_time)
  write.table(PQL_time, paste(DireOut, '/AAASimuTime_', SimuRange[1], 'to', SimuRange[2], '.txt', sep=''),
              sep='\t', row.names=F, col.names=T)
}

#run_PQL(DireData='simudata', DireOut='ResPQLNull', SimuRange=c(35, 55), setting='null', NumCore=8)


#-------------------------------------------------#
#        generate command file for MACAU          #
#-------------------------------------------------#
GenMacauCom <- function(ComMacau, ComData, ComOutDir, DireData, Setting){
  file_name <- read.table(paste(DireData, '/AAAParamsGrid.txt', sep=''), header=T)[, 1]
  main_com <- paste(ComMacau, '/macau', sep='')
  
  # generate sub command for each simulation dataset
  macau_com <- rep(0, length(file_name))
  for (i in 1:length(file_name)){
    if (Setting == 'null'){
      gcount_file <- paste(ComData, '/gcount_null_', file_name[i], '.txt', sep='')
    }
    if (Setting == 'alt'){
      gcount_file <- paste(ComData, '/gcount_alt_', file_name[i], '.txt', sep='')
    }
    
    tcount_file <- paste(ComData, '/tcountmacau_', file_name[i], '.txt', sep='')
    pred_file <- paste(ComData, '/pred_', file_name[i], '.txt', sep='')
    kernel_file <- paste(ComData, '/kernel_', file_name[i], '.txt', sep='')
    out_file <- paste('ResMACAU_', file_name[i], sep='')
    
    macau_com[i] <- paste(main_com, '  -g ', gcount_file, '  -t ', tcount_file, '  -p ', pred_file,
                          '  -k ', kernel_file, '  -pmm ', '  -o ', out_file, ' -outdir ', ComOutDir, sep='')
  }
  
  write.table(macau_com, paste(DireData, '/AAAComMacau_', file_name[1], '_to_', file_name[length(file_name)], '.txt', sep=''),
              sep='\t', row.names=F, col.names=F, quote=F)
  
}

# GenMacauCom(ComMacau='/public/home/SU_yinfei/code/MACAU130', ComData='/public/home/SU_yinfei/simudata/NullSeed1',
#             DireData='C:/Users/83944/Desktop', ComOutDir='/public/home/SU_yinfei/MacauNullSeed1', Setting='null')







