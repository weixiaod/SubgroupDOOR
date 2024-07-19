#' convert float to percentage
#'
#' @param num number as double
#' @param digits digits to round, default is 1 digit percentage
#'
#' @return a character string containing %
#' @noRd
#'
num_to_percent <- function(num,digits=2){
  input = paste0("%.",digits,"f")
  return(paste0(sprintf(input,num*100),"%"))
}
