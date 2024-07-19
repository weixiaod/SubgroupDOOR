#' Test the difference in DOOR probabilities by subgroups
#'
#' @param data the dataframe contains all variables
#' @param DOORlevel Total DOOR level K
#' @param trtvar the variable name for treatment indicator
#' @param DOORvar the variable name for DOOR
#' @param grpvar the variable name for subgroup indicator
#' @param grp_labels the subgroup labels
#' @param trt the treatment label
#' @param ctrl the control label
#'
#' @return p_value the p-value for the difference in DOOR probabilities
#' @export
#'
test_diff_doorprob <- function(data,DOORlevel,trtvar,DOORvar,grpvar,grp_labels=c(1,2),trt = "1",ctrl = "0"){
  # calculate the DOOR probability for each subgroup
  est = vector(length = length(grp_labels))
  var_est = vector(length = length(grp_labels))
  DOORtable_null <- summarize_DOOR(data=data,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)

  for(i in 1:length(grp_labels)){
    data.grp <- data %>% filter(.data[[grpvar]]==grp_labels[i])

    DOORtable <- summarize_DOOR(data=data.grp,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)
    Ntable <- DOORtable %>% group_by(.data[[trtvar]]) %>% summarise(N = sum(n))

    Pt0 <- dplyr::pull(DOORtable_null%>% filter(.data[[trtvar]]==1,),"Freq")
    Pc0 <- dplyr::pull(DOORtable_null%>% filter(.data[[trtvar]]==0,),"Freq")
    Nt <- dplyr::pull(Ntable %>% filter(.data[[trtvar]]==1,),"N")
    Nc <- dplyr::pull(Ntable %>% filter(.data[[trtvar]]==0,),"N")
    N <- c(Nt,Nc)
    K = DOORlevel

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
    var_est[i]=delta%*%cov%*%t(delta)

    est[i] = calculate_DOORprob_fromsummary(DOORtable=DOORtable,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)$`Prob[trt more desirable]`

  }


  # calculate the p-value
  z_stat = abs(est[1]-est[2])/sqrt(var_est[1]+var_est[2])
  p_value = 2*pnorm(-abs(z_stat))

  return(p_value)
}
