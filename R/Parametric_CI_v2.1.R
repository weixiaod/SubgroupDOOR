#' Calculate the DOOR probability from summary table
#'
#' @param DOORtable A dataframe contains the frequency crosstable of DOOR variable and Treatment variable
#' @param DOORvar Variable name of DOOR
#' @param trtvar Variable name of treatment
#' @param DOORlevel number of DOOR level pre-specified
#'
#' @param trt Character label for treatment in trtvar
#' @param ctrl character label for control in trtvar
#' @param alpha significance level
#' @param method calculate variance under null hypothesis or alternative hypothesis
#'
#' @return list of CIs
#' @export
#'
HalperinCI_summary <- function(DOORtable,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0",alpha=0.05,method=c("alternative","null")){
  #DOOR <-  as.numeric(dplyr::pull(data,DOORvar))
  #TRT <-  as.numeric(dplyr::pull(data,trtvar))
  if(missing(method)) method="alternative"
  method <- match.arg(method)

  Ntable <- DOORtable %>% group_by(.data[[trtvar]]) %>% summarise(N = sum(n))

  Pt <- dplyr::pull(DOORtable%>% filter(.data[[trtvar]]==trt,),"Freq")
  Pc <- dplyr::pull(DOORtable%>% filter(.data[[trtvar]]==ctrl,),"Freq")
  Nt <- dplyr::pull(Ntable %>% filter(.data[[trtvar]]==trt,),"N")
  Nc <- dplyr::pull(Ntable %>% filter(.data[[trtvar]]==ctrl,),"N")
  N <- c(Nt,Nc)
  K = DOORlevel

  empty=matrix(0,nrow=K,ncol=K)
  empty[lower.tri(empty)] <- 1
  indicator=diag(K)*0.5+empty
  est=sum((Pc%*%t(Pt))*indicator)

  name = paste("DOOR Probability of Group",trt,"more desirable than Group",ctrl,sep=" ")


  if(method=="null"){
    if(est == 1 | 0){
      message("Variance of DOOR prob is 0")
      return(list("Probability " = name,
                  "DOOR prob" = est,
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (1ST ORDER)" = c(est,est),
                  "NORMAL APPROXIMATION WITH CORRECTION (1ST ORDER)" = c(est,est),
                  "HARPERIN ET EL.(1989) WITHOUT CORRECTION" = c(est,est),
                  "HARPERIN ET EL.(1989) WITH CORRECTION" = c(est,est),
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (2ND ORDER)" = c(est,est),
                  "NORMAL APPROXIMATION WITH CORRECTION (2ND ORDER)" = c(est,est)
      ))

    } else{

      Rt = Nt/sum(N)
      Rc = Nc/sum(N)

      Pt0 = Pc0 = Rt*Pt + Rc*Pc

      ##replace diagonal elements of Cov with Variance##
      variance_c = Pc0*(1-Pc0)/Nc
      cov_c = -Pc0%*%t(Pc0)/Nc
      diag(cov_c)=variance_c

      variance_t = Pt0*(1-Pt0)/Nt
      cov_t = -Pt0%*%t(Pt0)/Nt
      diag(cov_t)=variance_t

      empty=matrix(0,nrow=K,ncol=K)
      cov = rbind(cbind(cov_c,empty),cbind(empty,cov_t))
      empty[lower.tri(empty)] <- 1
      indicator=diag(K)*0.5+empty

      delta_c <- rowSums(cbind((rep(1,K)%*%t(Pt0))*indicator))
      delta_t <- colSums(cbind(Pc0%*%t(rep(1,K))*indicator))
      delta <- cbind(t(delta_c),t(delta_t))
      testcov=delta%*%cov%*%t(delta)
      se=sqrt(testcov)

      lb_normal_fo <- est + qnorm(alpha/2)*se
      ub_normal_fo <- est - qnorm(alpha/2)*se

      lb_normal_fo_cc <- est + qnorm(alpha/2)*se - 1/(2*prod(N))
      ub_normal_fo_cc <- est - qnorm(alpha/2)*se + 1/(2*prod(N))

      ######Halperin CI
      M <- Pt0%*%t(Pc0)
      #est <- sum(sapply(1:K, function(j) dplyr::lead(M[j,], n = j, default = 0))) + sum(diag(M))/2

      A1  <-  sapply(2:K, function(j) ((sum(Pc0[j:K]) + Pc0[j-1]/2)^2))%*%Pt0[1:(K-1)]
      A2 <- Pt0[K]*Pc0[K]^2/4
      A_hat = A1 + A2
      B1 <- sapply(1:(K-1), function(j) ((sum(Pt0[1:j]) + Pt0[j+1]/2)^2))%*%Pc0[2:K]
      B2 <- Pt0[1]^2*Pc0[1]/4
      B_hat <- B1+B2

      AA_hat_1 = sapply(1:(K-1), function(j) (1-Pc0[j])*sum(Pc0[(j+1):K])-sum(Pc0[(j+1):K])^2)%*%Pt0[1:(K-1)]
      AA_hat_2 = sum(Pt0*Pc0*(1-Pc0))
      AA_hat = A_hat - (AA_hat_1)/(Nc-1) - (AA_hat_2)/(4*Nc-4)
      BB_hat_1 =  sapply(2:K, function(j) (1-Pt0[j])*sum(Pt0[1:(j-1)])-sum(Pt0[1:(j-1)])^2)%*%Pc0[2:K]
      BB_hat_2 = sum(Pt0*Pc0*(1-Pt0))
      BB_hat = B_hat - (BB_hat_1)/(Nt-1) -(BB_hat_2)/(4*Nt-4)


      theta_hat_top = (sum(N)-2)*est - (Nc-1)*AA_hat - (Nt-1)*BB_hat
      theta_hat_bot = ((prod(N)-sum(N)+2)*est - prod(N)*est^2)/((Nt-1)*(Nc-1)) + AA_hat/(Nt-1) + BB_hat/(Nc-1)
      theta_hat = theta_hat_top/((sum(N)-2)*theta_hat_bot)
      gamma_hat = (sum(N)-1) - (sum(N)-2)*theta_hat
      cc = gamma_hat*qchisq(1-alpha,df=1)/prod(N)

      theta_hat = ifelse(theta_hat<0,0,theta_hat)
      theta_hat = ifelse(theta_hat>1,1,theta_hat)

      #var_est <- (est-(sum(N)-1)*est^2+(Nc-1)*A_hat + (Nt-1)*B_hat - sum(diag(M))/4)/prod(N)
      var_est = (((sum(N)-1) - (sum(N)-2)*theta_hat)*est*(1-est))/prod(N)

      ## Normal approximation CI
      z_alpha = qnorm(1-alpha/2)
      lb_normal_so = est - z_alpha*sqrt(var_est)
      ub_normal_so = est + z_alpha*sqrt(var_est)

      ## Normal - Continuity correction CI
      lb_normal_so_cc = est - z_alpha*sqrt(var_est) - 1/(2*prod(N))
      ub_normal_so_cc = est + z_alpha*sqrt(var_est) + 1/(2*prod(N))

      ## Chi-square CI
      lb_chisq = (cc + 2*est - (cc^2+4*cc*est*(1-est))^(1/2))/(2*cc+2)
      ub_chisq = (cc + 2*est + (cc^2+4*cc*est*(1-est))^(1/2))/(2*cc+2)

      ## Chi-square continuity correction CI
      est1 = est - 1/(2*prod(N))
      est2 = est + 1/(2*prod(N))

      lb_chisq_cc = (cc + 2*est1 - (cc^2+4*cc*est1*(1-est1))^(1/2))/(2*cc+2)
      ub_chisq_cc = (cc + 2*est2 + (cc^2+4*cc*est2*(1-est2))^(1/2))/(2*cc+2)


      return(list("Probability " = name,
                  "DOOR prob" = est,
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (1ST ORDER)" = c(lb_normal_fo,ub_normal_fo),
                  "NORMAL APPROXIMATION WITH CORRECTION (1ST ORDER)" = c(lb_normal_fo_cc,ub_normal_fo_cc),
                  "HARPERIN ET EL.(1989) WITHOUT CORRECTION" = c(lb_chisq,ub_chisq),
                  "HARPERIN ET EL.(1989) WITH CORRECTION" = c(lb_chisq_cc,ub_chisq_cc),
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (2ND ORDER)" = c(lb_normal_so,ub_normal_so),
                  "NORMAL APPROXIMATION WITH CORRECTION (2ND ORDER)" = c(lb_normal_so_cc,ub_normal_so_cc))
      )
    }
  }
  if(method=="alternative"){
    if(est == 1 | 0){
      message("Variance of DOOR prob is 0")
      return(list("Probability " = name,
                  "DOOR prob" = est,
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (1ST ORDER)" = c(est,est),
                  "NORMAL APPROXIMATION WITH CORRECTION (1ST ORDER)" = c(est,est),
                  "HARPERIN ET EL.(1989) WITHOUT CORRECTION" = c(est,est),
                  "HARPERIN ET EL.(1989) WITH CORRECTION" = c(est,est),
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (2ND ORDER)" = c(est,est),
                  "NORMAL APPROXIMATION WITH CORRECTION (2ND ORDER)" = c(est,est)
      ))

    } else{

      ##replace diagonal elements of Cov with Variance##
      variance_c = Pc*(1-Pc)/Nc
      cov_c = -Pc%*%t(Pc)/Nc
      diag(cov_c)=variance_c

      variance_t = Pt*(1-Pt)/Nt
      cov_t = -Pt%*%t(Pt)/Nt
      diag(cov_t)=variance_t

      empty=matrix(0,nrow=K,ncol=K)
      cov = rbind(cbind(cov_c,empty),cbind(empty,cov_t))
      empty[lower.tri(empty)] <- 1
      indicator=diag(K)*0.5+empty

      delta_c <- rowSums(cbind((rep(1,K)%*%t(Pt))*indicator))
      delta_t <- colSums(cbind(Pc%*%t(rep(1,K))*indicator))
      delta <- cbind(t(delta_c),t(delta_t))
      testcov=delta%*%cov%*%t(delta)
      se=sqrt(testcov)

      lb_normal_fo <- est + qnorm(alpha/2)*se
      ub_normal_fo <- est - qnorm(alpha/2)*se

      lb_normal_fo_cc <- est + qnorm(alpha/2)*se - 1/(2*prod(N))
      ub_normal_fo_cc <- est - qnorm(alpha/2)*se + 1/(2*prod(N))

      ######Halperin CI
      M <- Pt%*%t(Pc)
      #est <- sum(sapply(1:K, function(j) dplyr::lead(M[j,], n = j, default = 0))) + sum(diag(M))/2

      A1  <-  sapply(2:K, function(j) ((sum(Pc[j:K]) + Pc[j-1]/2)^2))%*%Pt[1:(K-1)]
      A2 <- Pt[K]*Pc[K]^2/4
      A_hat = A1 + A2
      B1 <- sapply(1:(K-1), function(j) ((sum(Pt[1:j]) + Pt[j+1]/2)^2))%*%Pc[2:K]
      B2 <- Pt[1]^2*Pc[1]/4
      B_hat <- B1+B2

      AA_hat_1 = sapply(1:(K-1), function(j) (1-Pc[j])*sum(Pc[(j+1):K])-sum(Pc[(j+1):K])^2)%*%Pt[1:(K-1)]
      AA_hat_2 = sum(Pt*Pc*(1-Pc))
      AA_hat = A_hat - (AA_hat_1)/(Nc-1) - (AA_hat_2)/(4*Nc-4)
      BB_hat_1 =  sapply(2:K, function(j) (1-Pt[j])*sum(Pt[1:(j-1)])-sum(Pt[1:(j-1)])^2)%*%Pc[2:K]
      BB_hat_2 = sum(Pt*Pc*(1-Pt))
      BB_hat = B_hat - (BB_hat_1)/(Nt-1) -(BB_hat_2)/(4*Nt-4)


      theta_hat_top = (sum(N)-2)*est - (Nc-1)*AA_hat - (Nt-1)*BB_hat
      theta_hat_bot = ((prod(N)-sum(N)+2)*est - prod(N)*est^2)/((Nt-1)*(Nc-1)) + AA_hat/(Nt-1) + BB_hat/(Nc-1)
      theta_hat = theta_hat_top/((sum(N)-2)*theta_hat_bot)
      gamma_hat = (sum(N)-1) - (sum(N)-2)*theta_hat
      cc = gamma_hat*qchisq(1-alpha,df=1)/prod(N)

      theta_hat = ifelse(theta_hat<0,0,theta_hat)
      theta_hat = ifelse(theta_hat>1,1,theta_hat)

      #var_est <- (est-(sum(N)-1)*est^2+(Nc-1)*A_hat + (Nt-1)*B_hat - sum(diag(M))/4)/prod(N)
      var_est = (((sum(N)-1) - (sum(N)-2)*theta_hat)*est*(1-est))/prod(N)

      ## Normal approximation CI
      z_alpha = qnorm(1-alpha/2)
      lb_normal_so = est - z_alpha*sqrt(var_est)
      ub_normal_so = est + z_alpha*sqrt(var_est)

      ## Normal - Continuity correction CI
      lb_normal_so_cc = est - z_alpha*sqrt(var_est) - 1/(2*prod(N))
      ub_normal_so_cc = est + z_alpha*sqrt(var_est) + 1/(2*prod(N))

      ## Chi-square CI
      lb_chisq = (cc + 2*est - (cc^2+4*cc*est*(1-est))^(1/2))/(2*cc+2)
      ub_chisq = (cc + 2*est + (cc^2+4*cc*est*(1-est))^(1/2))/(2*cc+2)

      ## Chi-square continuity correction CI
      est1 = est - 1/(2*prod(N))
      est2 = est + 1/(2*prod(N))

      lb_chisq_cc = (cc + 2*est1 - (cc^2+4*cc*est1*(1-est1))^(1/2))/(2*cc+2)
      ub_chisq_cc = (cc + 2*est2 + (cc^2+4*cc*est2*(1-est2))^(1/2))/(2*cc+2)


      return(list("Probability " = name,
                  "DOOR prob" = est,
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (1ST ORDER)" = c(lb_normal_fo,ub_normal_fo),
                  "NORMAL APPROXIMATION WITH CORRECTION (1ST ORDER)" = c(lb_normal_fo_cc,ub_normal_fo_cc),
                  "HARPERIN ET EL.(1989) WITHOUT CORRECTION" = c(lb_chisq,ub_chisq),
                  "HARPERIN ET EL.(1989) WITH CORRECTION" = c(lb_chisq_cc,ub_chisq_cc),
                  "NORMAL APPROXIMATION WITHOUT CORRECTION (2ND ORDER)" = c(lb_normal_so,ub_normal_so),
                  "NORMAL APPROXIMATION WITH CORRECTION (2ND ORDER)" = c(lb_normal_so_cc,ub_normal_so_cc))
      )
    }
  }
}

#' Calculate DOOR prob confidence intervals under null hypotheses
#'
#' Calculate the DOOR probability from summary table
#'
#' @param data A dataframe contains DOOR variable and Treatment variable
#' @param DOORvar Variable name of DOOR
#' @param trtvar Variable name of treatment
#' @param DOORlevel number of DOOR level pre-specified
#'
#' @param trt Character label for treatment in trtvar
#' @param ctrl character label for control in trtvar
#' @param alpha significance level
#' @param method calculate variance under null hypothesis or alternative hypothesis
#'
#' @return list of DOOR probability and CIs
#' @export
#'
Halperin_CI <- function(data,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0",alpha=0.05,method=c("null","alternative")){
  if(missing(method)) method="alternative"
  method <- match.arg(method)
  DOORtable  <-  summarize_DOOR(data=data,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)
  HalperinCI_summary(DOORtable=DOORtable,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl,alpha = alpha,method=method)
}
