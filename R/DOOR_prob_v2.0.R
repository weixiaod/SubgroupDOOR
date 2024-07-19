#' Summarize DOOR variable by treatment groups
#'
#' @param data A dataframe contains DOOR variable and Treatment variable
#' @param DOORvar Variable name of DOOR
#' @param trtvar Variable name of treatment
#' @param DOORlevel number of DOOR level pre-specified
#'
#' @param trt Char label for treatment indicator
#' @param ctrl char label for control indicator
#'
#' @export
#'
summarize_DOOR <- function(data,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0"){

  if(missing(DOORlevel)) {
    stop("DOOR level is a required argument")
  }

  ARM <-  pull(data,trtvar)
  DOOR <-  as.numeric(pull(data,DOORvar))

  na_id = is.na(DOOR) | is.na(ARM)


  if(any(na_id)) message(paste("Warning: Missing data in",DOORvar,"\nDOOR probability is still calculated by removing all missing"))

  if(DOORlevel >= length(unique(na.omit(DOOR)))){
      data.f <- data %>%
        filter(!na_id) %>%
        mutate(DOOR.f =  factor(.data[[DOORvar]],levels = c(1:DOORlevel)))
  } else {
    stop(paste("DOOR levels in data is greater than input DOOR level ",DOORlevel))

  }

  dist.table <- data.f %>%
    group_by(.data[[trtvar]], .drop = FALSE) %>%
    count(DOOR.f) %>%
    mutate(Freq = n/sum(n)) %>%
    rename(!!DOORvar:= DOOR.f)

  return(dist.table)
}

#' Calculate the DOOR probability from summary table
#'
#' @param DOORtable A dataframe contains the frequency crosstable of DOOR variable and Treatment variable
#' @param DOORvar Variable name of DOOR
#' @param trtvar Variable name of treatment
#' @param DOORlevel number of DOOR level pre-specified
#'
#' @param trt Character label for treatment in trtvar
#' @param ctrl character label for control in trtvar
#'
#' @export

calculate_DOORprob_fromsummary <- function(DOORtable,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0"){

  N <- DOORtable %>% group_by(.data[[trtvar]]) %>% summarise(N = sum(n))

  Pt <- pull(DOORtable%>% filter(.data[[trtvar]]==trt,),"Freq")
  Pc <- pull(DOORtable%>% filter(.data[[trtvar]]==ctrl,),"Freq")
  Nt <- pull(N %>% filter(.data[[trtvar]]==trt,),"N")
  Nc <- pull(N %>% filter(.data[[trtvar]]==ctrl,),"N")

  K = DOORlevel

  M <- sapply(1:K, function(i) lead(Pc, n = i, default = 0))

  prob_larger = sum(M%*%Pt)
  prob_equal = sum(Pt*Pc)

  doorprob = prob_larger + 0.5 * prob_equal

  return(list("Prob[trt > ctrl]"=prob_larger,
              "Prob[trt = ctrl]"=prob_equal,
              "Prob[trt < ctrl]"=1-prob_larger-prob_equal,
              "Prob[trt more desirable]"=doorprob,
              "Prob[ctrl more desirable]"=1-doorprob)
         )
}

#' Calculate the DOOR probability
#'
#' @param DOORvar Variable name of DOOR
#' @param trtvar Variable name of treatment
#' @param DOORlevel number of DOOR level pre-specified
#'
#' @param trt Character label for treatment in trtvar
#' @param ctrl character label for control in trtvar
#' @param data dataframe that contains required variables
#'
#' @return list of DOOR probability outputs
#' @export

calculate_DOORprob <- function(data,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0"){
  DOORtable  <-  summarize_DOOR(data=data,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)
  calculate_DOORprob_fromsummary(DOORtable=DOORtable,DOORvar=DOORvar,trtvar=trtvar,DOORlevel=DOORlevel,trt=trt,ctrl=ctrl)

}
