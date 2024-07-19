#' Subgroup identification function
#'
#' @param data the data frame
#' @param DOORvar specify the DOORvar name such as "DOOR"
#' @param trtvar specify the treatment variable name such as "ARMN"
#' @param covar specify the covariate name such as "covariate"
#' @param DOORlevel numeric value of the DOOR level
#' @param trt treatment indicator
#' @param ctrl control indicator
#' @param minsize minimum size of the split
#' @param incre increament of the split
#'
#' @return a vector of the split range
#' @export

cartdoor <- function(data,DOORvar,trtvar,covar,DOORlevel,trt = "1", ctrl = "0",minsize = 30,incre = 5,midgroup = FALSE,alpha = 0.05){
  # range of covariate
  covariate = dplyr::pull(data,covar)
  range <- range(covariate)
  range_step <- seq(range[1],range[2],by = incre)
  length <- length(range_step)
  sigpvalue = 1


  if(midgroup){
    for(i in 1:length){
      for(j in i:length){
        # split the data
        ind = (covariate > range_step[i] & covariate < range_step[j])
        if(sum(ind) > minsize & sum(!ind) > minsize) {
          data1 = data[ind,]
          data2 = data[!ind,]
          DOORprob1 = calculate_DOORprob(data = data1,DOORvar = DOORvar,trtvar = trtvar,
                                         DOORlevel = DOORlevel,trt = trt, ctrl = ctrl)$`Prob[trt more desirable]`
          DOORprob2 = calculate_DOORprob(data = data2,DOORvar = DOORvar,trtvar = trtvar,
                                         DOORlevel = DOORlevel,trt = trt, ctrl = ctrl)$`Prob[trt more desirable]`
          diff = abs(DOORprob1 - DOORprob2)

          # give a group label by ind
          data_split = data
          data_split$stratum = ifelse(ind,1,2)
          p_value = test_diff_doorprob(data = data_split,DOORlevel = DOORlevel, DOORvar = DOORvar,grpvar = "stratum", trtvar = trtvar,grp_labels = c(1,2))

          if (p_value < sigpvalue) {
            sigpvalue = p_value
            splitrange = c(range_step[i],range_step[j])
            DPrecord1 = DOORprob1
            DPrecord2 = DOORprob2
          }
        }
      }
    }
    } else {
      for(i in 2:length){

          ind = data[,covar] < range_step[i]
          if(sum(ind) > minsize & sum(!ind) > minsize) {
            data1 = data[ind,]
            data2 = data[!ind,]
            DOORprob1 = calculate_DOORprob(data = data1,DOORvar = DOORvar,trtvar = trtvar,
                                           DOORlevel = DOORlevel,trt = trt, ctrl = ctrl)$`Prob[trt more desirable]`
            DOORprob2 = calculate_DOORprob(data = data2,DOORvar = DOORvar,trtvar = trtvar,
                                           DOORlevel = DOORlevel,trt = trt, ctrl = ctrl)$`Prob[trt more desirable]`
            #diff = abs(DOORprob1 - DOORprob2)
            #if (maxdiff < diff) {
              #maxdiff = diff
            # give a group label by ind
            data_split = data
            data_split$stratum = ifelse(ind,1,2)
            p_value = test_diff_doorprob(data = data_split,DOORlevel = DOORlevel, DOORvar = DOORvar,grpvar = "stratum", trtvar = trtvar,grp_labels = c(1,2))

            if (p_value < sigpvalue) {
              sigpvalue = p_value
              splitrange = c(range_step[1],range_step[i])
              DPrecord1 = DOORprob1
              DPrecord2 = DOORprob2
            }
          }
        }
    }
  if(sigpvalue < alpha){
  return(list(splitrange = splitrange,
              p.value =  sigpvalue,
              grp1prob = DPrecord1,
              grp2prob = DPrecord2))
  } else {
    return("No significant split found")
  }
}



