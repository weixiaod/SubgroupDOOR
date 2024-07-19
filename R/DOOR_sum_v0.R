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
#'
#'
get_DOOR_sum <- function(data,DOORvar,trtvar,DOORlevel,trt="1",ctrl="0"){

  if(missing(DOORlevel)) {
    stop("DOORlevel is a required argument")
  }

  ARM <-  dplyr::pull(data,trtvar)
  DOOR <-  as.numeric(dplyr::pull(data,DOORvar))

  na_id = is.na(DOOR) | is.na(ARM)


  if(any(na_id)) message(paste("Warning: Missing data in data structure","\nDOOR probability is still calculated by removing all missing"))

  if(DOORlevel >= length(unique(na.omit(DOOR)))){
    data_f <- data %>%
      filter(!na_id) %>%
      mutate(DOOR_f =  factor(.data[[DOORvar]],levels = c(1:DOORlevel)))
  } else {
    stop(paste("DOOR levels in data is greater than input DOOR level ",DOORlevel))

  }

  DOOR_dist <- data_f %>%
    group_by(.data[[trtvar]], .drop = FALSE) %>%
    count(DOOR_f) %>%
    mutate(Freq = n/sum(n)) %>%
    rename(!!DOORvar:= DOOR_f) %>% ungroup()

  #Create S3 object - DOORsum
  DOOR_sum <- structure(list(),class="DOORsum")

  DOOR_sum[["label"]] <-   c(trt,ctrl)
  names(DOOR_sum[["label"]]) <- c("treatment","control")

  DOOR_sum[["treatment"]] <- DOOR_dist %>% filter(.data[[trtvar]]==trt)
  DOOR_sum[["control"]] <- DOOR_dist %>% filter(.data[[trtvar]]==ctrl)
  DOOR_sum[["n"]] <-   c(sum(DOOR_sum[["treatment"]]$n),sum(DOOR_sum[["control"]]$n))
  names(DOOR_sum[["n"]]) <- c("treatment","ctrl")

  attr(DOOR_sum,"DOOR") <- DOORvar
  attr(DOOR_sum,"ARM") <- trtvar
  attr(DOOR_sum,"DOORlevel") <- DOORlevel
  attr(DOOR_sum,"treatment") <- trt
  attr(DOOR_sum,"control") <- ctrl


  return(DOOR_sum)
}


#'  Print method for S3 object DOORsum
#'
#' @param obj S3 object DOORsum
#'
#' @export
#'
print.DOORsum <- function(obj) {
  #cat("# Active Group",obj$label[1], "\t Sample size =", obj$n[1],"\n")
  #cat("# Control Group",obj$label[2], "\t Sample size =", obj$n[2],"\n\n")

  cat("# DOOR Distribution\n")
  cat("## Active Group",obj$label[1], "\t Sample size =", obj$n[1])

  print(kable(obj$treatment %>%
          mutate(proportion = scales::percent(Freq,accuracy = 0.1)) %>%
          select(DOOR,n,proportion) %>%
          rename("DOOR rank" = DOOR,
                 !!paste("Group",obj$label[1],"N"):=n)
        ))

  cat("\n## Control Group",obj$label[2], "\t Sample size =", obj$n[2])
  print(kable(obj$control %>%
          mutate(proportion = scales::percent(Freq,accuracy = 0.1)) %>%
          select(DOOR,n,proportion) %>%
          rename("DOOR rank" = DOOR,
                 !!paste("Group",obj$label[2],"N"):=n)
  ))
}

#print(DOOR_sum)
