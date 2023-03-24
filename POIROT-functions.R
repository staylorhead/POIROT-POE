
###############################################################################
# function to center phenotypes by genotype
###############################################################################
center_pheno <- function(dat, by){
  npheno <- length(grep("^X",colnames(dat)))
  
  if(by=="mean"){
    l <- lapply(1:npheno, function(x) {group.center(dat[,paste0("X",x)], 
                                                    dat$off_geno_grp)}) # center at mean 
    dat <- cbind(dat,data.frame(t(matrix(unlist(l), nrow=length(l), byrow=T))))
    colnames(dat)[(ncol(dat)-npheno+1):ncol(dat)] <- paste0("pheno",1:npheno)
  }else if(by=="median"){
    l <- lapply(1:npheno, function(x) {
      tapply(dat[,paste0("X",x)],dat$off_geno_grp,median) # center at median
    })
    
    dat <- cbind(dat,dat[,grep("^X",colnames(dat))])
    colnames(dat)[(ncol(dat)-npheno+1):ncol(dat)] <- paste0("pheno",1:npheno)
    
    cent_phenos <- sapply(1:npheno,function(i){
      tmp <- {}
      tmp[dat$off_geno_grp=="AA"] <- dat[dat$off_geno_grp=="AA",paste0("X",i)]-l[[i]][1]
      tmp[dat$off_geno_grp=="Aa"] <- dat[dat$off_geno_grp=="Aa",paste0("X",i)]-l[[i]][2]
      tmp[dat$off_geno_grp=="aa"] <- dat[dat$off_geno_grp=="aa",paste0("X",i)]-l[[i]][3]
      return(tmp)
    })
    dat[,paste0("pheno",1:npheno)] <- cent_phenos
  }else if(by=="none"){
    colnames(dat)[grep("^X",colnames(dat))] <- paste0("pheno",1:npheno)
  }
  return(dat)
}

###############################################################################
# function to adjust for effects of covariates
###############################################################################
extract_residuals <- function(phen_df, covar_df){
  npheno <- ncol(phen_df)
  resid_mat = matrix(NA, nrow=nrow(phen_df), ncol=npheno)
  for(l in 1:npheno){
    tmp <- data.frame(pheno=phen_df[,l])
    tmp <- cbind(tmp,covar_df)
    fit <- lm(pheno~.,data=tmp)
    resid_mat[,l] <- resid(fit)
  }
  resid_mat <- data.frame(resid_mat)
  colnames(resid_mat) <- colnames(phen_df)
  return(resid_mat)
}

################################################################################
# function to R-Omnibus test for equality of phenotypic covariance matrices
################################################################################
do_r_omnibus_test <- function(dat, varnames, groupname){
  dat <- dat[,c(groupname,varnames)]
  npheno <- length(varnames)
  grplevels <- levels(dat[,groupname])
  ngrp1 <- sum(dat[,groupname]==grplevels[1])
  ngrp2 <- sum(dat[,groupname]==grplevels[2])
  
  # grab medians of each variable in each group
  M <- sapply(1:npheno, function(x) tapply(dat[,varnames[x]], dat[,groupname], median))
  
  X <- vector(mode = "list", length = 2)
  X[[1]] <- as.matrix(dat[dat[,groupname]==grplevels[1],varnames])
  X[[2]] <- as.matrix(dat[dat[,groupname]==grplevels[2],varnames])
  
  x_M_1 <- sapply(1:npheno, function(x) X[[1]][,x]-M[1,x])
  x_M_2 <- sapply(1:npheno, function(x) X[[2]][,x]-M[2,x])
  
  Z <- vector(mode = "list", length = 2)
  Z[[1]] <- vector(mode="list",length=ngrp1)
  Z[[2]] <- vector(mode="list",length=ngrp2)
  Z[[1]] <- lapply(1:ngrp1, function(j) unlist(sapply(1:npheno, function(k) sapply(k:npheno, function(kk) x_M_1[j,k]*x_M_1[j,kk]) ) ))
  Z[[2]] <- lapply(1:ngrp2, function(j) unlist(sapply(1:npheno, function(k) sapply(k:npheno, function(kk) x_M_2[j,k]*x_M_2[j,kk]) ) )) 
  W <- vector(mode = "list", length = 2)
  W[[1]] <- vector(mode="list",length=ngrp1)
  W[[2]] <- vector(mode="list",length=ngrp2)
  
  W[[1]] <- lapply(1:ngrp1, function(j) {
    v <- unlist(Z[[1]][j])/(abs(unlist(Z[[1]][j]))^(1/2))
    v[is.nan(v)] <- 0
    return(v)
  })
  W[[2]] <- lapply(1:ngrp2, function(j) {
    v <- unlist(Z[[2]][j])/(abs(unlist(Z[[2]][j]))^(1/2))
    v[is.nan(v)] <- 0
    return(v)
  })
  
  W_df <- data.frame(grp=c(rep(grplevels[1],ngrp1),rep(grplevels[2],ngrp2)))
  W_df$grp <- factor(W_df$grp)
  ntest <- (npheno^2+npheno)/2
  
  mat_1 <- matrix(unlist(W[[1]]),nrow=ngrp1,ncol=ntest,byrow=T)
  mat_2 <- matrix(unlist(W[[2]]),nrow=ngrp2,ncol=ntest,byrow=T)
  mat <- rbind(mat_1,mat_2)
  W_df <- cbind(W_df,mat)
  
  fit <- manova(as.matrix(W_df[,-1])~grp, data=W_df)
  
  pval <- summary(fit)$stats[1,"Pr(>F)"]
  stat <- summary(fit)$stats[1,"approx F"]
  
  return(c(pval=pval,stat=stat))
}

################################################################################
# function to perform POIROT test for one variant
################################################################################
do_POIROT_by_snp <- function(i,phenodat,genodat){
  dat <- phenodat
  dat$geno <- genodat[,i]
  dat$geno_grp <- NA
  dat$geno_grp[dat$geno%in%c(0,2)] <- "Homozygote"
  dat$geno_grp[dat$geno==1] <- "Heterozygote"
  dat$geno_grp <- factor(dat$geno_grp,levels=c("Homozygote","Heterozygote"))
  npheno <- ncol(phenodat)
  
  # center by median among genotypes groups (0,1,2)
  medians <- sapply(1:npheno, function(x) tapply(dat[,x], dat$geno, median))
  tmp <- {}
  for(j in 1:npheno){
    tmp <- data.frame(cbind(tmp,dat[,j]-medians[dat$geno+1,j]))
  }
  colnames(tmp) <- paste0("pheno",1:npheno,"_cent")
  dat <- cbind(dat,tmp)
  tmp <- do_r_omnibus_test(dat=dat, varnames=paste0("pheno",1:npheno,"_cent"), groupname="geno_grp")
  
  return(tmp)
}
