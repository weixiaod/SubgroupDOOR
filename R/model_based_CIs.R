
#' Calculate the DOOR probabilities for new dataset
#'
#' @param the dataset that contains the predicted data from fit.clm
#'
#' @return DOOR probability for each row
#' @noRd

get_DOORprob_cont <- function(predictdata) {
  # separate the data into two groups to ensure that the concatenated data has the desired order
  # which is for each row, is c(Pt, Pc)
  data_PT <- predictdata %>%
    filter(ARMN == 1) %>%
    dplyr::select(-ARMN)
  data_PC <- predictdata %>%
    filter(ARMN == 0) %>%
    dplyr::select(-ARMN)
  data_both <- merge(data_PT, data_PC, by = "covariate", prefix = c("trt_", "ctrl_"))

  DOORprob <- apply(data_both[, -1], 1, DOORprob_byrow, DOORlevel = DOORlevel)
  return(cbind(data_both, DOORprob))
}


#' get the derivative given covariate
#'
#' @param x the covariate value
#' @param fit.clm the model fit
#' @param DOORlevel Total DOOR level K
#'
#' @return the derivative of the DOOR prob given covariate
#' @noRd
getDeri_clm <- function(x, fit.clm, DOORlevel) {
  # get the derivative of the DOOR prob given beta
  # x is the covariate value
  # fit.clm is the model fit
  m.coef <- matrix(fit.clm$coefficients, ncol = 4)

  newdata.clm <- data.frame(
    intercept = 1,
    ARMN = c(1, 0),
    covariate = x,
    interaction = c(1, 0) * x
  )
  prob <- t(predict(fit.clm, newdata = newdata.clm)$fit)

  # h: (DOORlevel-1) by 2 matrix
  h <- exp(m.coef %*% t(newdata.clm))
  # hh: (DOORlevel-1) by 2 matrix for derivative of P(Y<=k)
  hh <- h / ((1 + h)^2)

  # coefficient order:
  # 1: intercept
  # 2: ARMN
  # 3: covariate
  # 4: interaction
  der_1 <- vector(length = DOORlevel - 1)
  der_2 <- vector(length = DOORlevel - 1)
  der_3 <- vector(length = DOORlevel - 1)
  der_4 <- vector(length = DOORlevel - 1)

  for (k in 1:(DOORlevel - 1)) {
    der_1[k] <- hh[k, 1] * (prob[k, 2] + prob[k + 1, 2]) / 2 - hh[k, 2] * (prob[k, 1] + prob[k + 1, 1]) / 2
    der_2[k] <- hh[k, 1] * (prob[k + 1, 2] + prob[k, 2]) / 2
  }
  der_3 <- der_1 * x
  der_4 <- der_2 * x
  return(c(der_1, der_2, der_3, der_4))
}


#' get the variance-covariance matrix for the parameters
#'
#' @param covariate the covariate value
#' @param fit.clm the model fit
#' @param DOORlevel Total DOOR level K
#'
#' @return the variance-covariance matrix for the parameters
#' @noRd
conditionalSE_clm <- function(covariate, fit.clm, DOORlevel) {
  # get the conditional CI for the DOOR prob given beta
  # x is the covariate value
  # fit.model1 is the model fit for ARMN = 1
  # fit.model0 is the model fit for ARMN = 0
  # vcov_both is the combined var-cov matrix
  vcov <- fit.clm$vcov
  der <- getDeri_clm(x = covariate, fit.clm, DOORlevel)
  # get the conditional CI for the DOOR prob given beta
  vcov_DOORprob <- t(der) %*% vcov %*% (der)
  # get the conditional CI for the DOOR prob given beta
  return(sqrt(vcov_DOORprob))
}


#' Calculate the pointwise confidence interval for DOOR probabilities for new dataset
#'
#' @param fit.clm the model fit
#' @param newdata the new dataset contains arm indicator and covariate
#' @param DOORlevel Total DOOR level K
#'
#' @return a dataframe contains the covariate, DOOR probability, and confidence interval
#' @export
#'

getCI_clm <- function(fit.clm, newdata, DOORlevel) {
  se_covariate <- sapply(unique(sort(newdata$covariate)), conditionalSE_clm, fit.clm = fit.clm, DOORlevel = DOORlevel)
  data.pred <- cbind(newdata, predict(fit.clm, newdata = newdata))
  out.pred <- get_DOORprob_cont(data.pred)[, "DOORprob"]
  out.pred_lb <- sapply(out.pred - 1.96 * se_covariate, max, 0)
  out.pred_ub <- sapply(out.pred + 1.96 * se_covariate, min, 1)

  # get confidence interval by logit transformation
  # get CI for log odds of DOOR prob
  out.logodds_pred_lb <- log(out.pred/(1-out.pred)) - 1.96 * se_covariate / out.pred/(1-out.pred)
  out.logodds_pred_ub <- log(out.pred/(1-out.pred)) + 1.96 * se_covariate / out.pred/(1-out.pred)

  # get simultaneous confidence band by Scheffe's method
  out.pred_simullb <- sapply(out.pred - qnorm(1-0.05/200) * se_covariate, max, 0)
  out.pred_simulub <- sapply(out.pred + qnorm(1-0.05/200) * se_covariate, min, 1)
  plotdata <- data.frame(
    covariate = 1:100,
    DOORprob = out.pred,
    DOORprob_lb = out.pred_lb,
    DOORprob_ub = out.pred_ub,
    DOORprob_lb_lo = exp(out.logodds_pred_lb)/(1+exp(out.logodds_pred_lb)),
    DOORprob_ub_lo = exp(out.logodds_pred_ub)/(1+exp(out.logodds_pred_ub)),
    DOORprob_simullb = out.pred_simullb,
    DOORprob_simulub = out.pred_simulub,
    warningCode = fit.clm$convergence$code

  )
  return(plotdata)
}

###############################################################################
#######  For Proportional odds assumption                             #########
###############################################################################
#' calculate the derivative of covariate
#'
#' @param x covariate value
#' @param fit.po the model fit
#' @param DOORlevel Total DOOR level K
#'
#' @return the derivative of the DOOR prob given covariate
#' @noRd
#'
getDeri_po <- function(x, fit.po, DOORlevel) {
  # get the derivative of the DOOR prob given beta
  # x is the covariate value
  # fit.po is the model fit
  coef_intercept <- t(fit.po$coefficients[1:(DOORlevel - 1)])
  m.coef <- matrix(
    cbind(
      fit.po$coefficients[1:(DOORlevel - 1)],
      -rep(fit.po$coefficients[DOORlevel], DOORlevel - 1),
      -rep(fit.po$coefficients[DOORlevel + 1], DOORlevel - 1),
      -rep(fit.po$coefficients[DOORlevel + 2], DOORlevel - 1)
    ),
    ncol = 4
  )

  newdata.po <- data.frame(
    intercept = 1,
    ARMN = c(1, 0),
    covariate = x,
    interaction = c(1, 0) * x
  )
  prob <- t(predict(fit.po, newdata = newdata.po)$fit)

  # h: (DOORlevel-1) by 2 matrix
  h <- exp((m.coef) %*% t(newdata.po))
  # hh: (DOORlevel-1) by 2 matrix for derivative of P(Y<=k)
  hh <- h / ((1 + h)^2)
  hhh <- hh - lag(hh, default = 0)
  hhhh <- hh + lag(hh, default = 0)
  # coefficient order:
  # 1: intercept
  # 2: ARMN
  # 3: covariate
  # 4: interaction
  der_1 <- vector(length = DOORlevel - 1)
  der_2 <- 0
  der_3 <- 0
  der_4 <- 0

  for (k in 1:(DOORlevel - 1)) {
    der_1[k] <- hh[k, 1] * (prob[k, 2] + prob[(k + 1), 2]) / 2 - hh[k, 2] * (prob[k, 1] + prob[(k + 1), 1]) / 2

    der_2 <- der_2 + hhh[k, 1] * sum(prob[(k + 1):DOORlevel, 2], prob[k, 2] / 2)

    der_3 <- der_3 + hhh[k, 1] * sum(prob[(k + 1):DOORlevel, 2], prob[k, 2] / 2) - prob[k, 1] * hhhh[k, 2] / 2
  }
  der_2 <- der_2 - hh[DOORlevel - 1, 1] * prob[DOORlevel, 2] / 2
  der_3 <- (der_3 - hh[DOORlevel - 1, 1] * prob[DOORlevel, 2] / 2 - hh[DOORlevel - 1, 2] * prob[DOORlevel, 1] / 2) * x
  der_4 <- der_2 * x
  return(c(der_1, der_2, der_3, der_4))
}

#' calcluate the conditional standard error for the DOOR probability
#'
#' @param covariate the covaraite value
#' @param fit.po the proportional odds model fit
#' @param DOORlevel Total DOOR level K
#'
#' @return the conditional standard error for the DOOR probability
#' @noRd
conditionalSE_po <- function(covariate, fit.po, DOORlevel) {
  # get the conditional CI for the DOOR prob given beta
  # x is the covariate value
  # fit.model1 is the model fit for ARMN = 1
  # fit.model0 is the model fit for ARMN = 0
  # vcov_both is the combined var-cov matrix
  vcov <- fit.po$vcov
  der <- getDeri_po(x = covariate, fit.po, DOORlevel)
  # get the conditional CI for the DOOR prob given beta
  vcov_DOORprob <- t(der) %*% vcov %*% (der)
  # get the conditional CI for the DOOR prob given beta
  return(sqrt(vcov_DOORprob))
}



#' Calculate the pointwise confidence interval for DOOR probabilities for new dataset
#'
#' @param fit.po the fitted proportional odds model
#' @param newdata the new dataset contains arm indicator and covariate
#' @param DOORlevel Total DOOR level K
#'
#' @return a dataframe contains the covariate, DOOR probability, and confidence interval
#' @export
#'
getCI_po <- function(fit.po, newdata, DOORlevel) {
  se_covariate <- sapply(unique(sort(newdata$covariate)), conditionalSE_po, fit.po = fit.po, DOORlevel = DOORlevel)
  data.pred <- cbind(newdata, predict(fit.po, newdata = newdata))
  out.pred <- get_DOORprob_cont(data.pred)[, "DOORprob"]
  out.pred_lb <- sapply(out.pred - 1.96 * se_covariate, max, 0)
  out.pred_ub <- sapply(out.pred + 1.96 * se_covariate, min, 1)
  # get confidence band by logit transformation
  out.logodds_pred_lb <- log(out.pred/(1-out.pred)) - 1.96 * se_covariate / out.pred/(1-out.pred)
  out.logodds_pred_ub <- log(out.pred/(1-out.pred)) + 1.96 * se_covariate / out.pred/(1-out.pred)

  # get simultaneous confidence band by Scheffe's method
  out.pred_simullb <- sapply(out.pred - qnorm(1-0.05/200) * se_covariate, max, 0)
  out.pred_simulub <- sapply(out.pred + qnorm(1-0.05/200) * se_covariate, min, 1)

  plotdata <- data.frame(
    covariate = 1:100,
    DOORprob = out.pred,
    DOORprob_lb = out.pred_lb,
    DOORprob_ub = out.pred_ub,
    DOORprob_lb_lo = exp(out.logodds_pred_lb) / (1 + exp(out.logodds_pred_lb)),
    DOORprob_ub_lo = exp(out.logodds_pred_ub) / (1 + exp(out.logodds_pred_ub)),
    DOORprob_simullb = out.pred_simullb,
    DOORprob_simulub = out.pred_simulub,
    warningCode = fit.po$convergence$code
  )
  return(plotdata)
}



