#' @title   Converging DOOR categorical variable to a 1/0 matrix
#' @param DOOR A vector contains DOOR outcome for all subjects
#' @param nlevel defined by user: how many DOOR levels is used in analysis

#'
#' @return  A binary matrix which contains the DOOR category information
#' @export


DOOR_to.dummy <- function(DOOR,nlevel){
  DOOR = factor(DOOR,levels=c(1:nlevel))
  DOOR_level = levels(DOOR)
  n = length(unique(DOOR))
  p = length(DOOR)
  DOOR_matrix = matrix(0,nrow = p,ncol = nlevel)
  var_name=vector(length = nlevel)
  for (i in 1:nlevel){
    DOOR_matrix[,i]=(DOOR==DOOR_level[i])
    var_name[i] = paste("DOOR_",i,sep = "")
  }
  colnames(DOOR_matrix) <- var_name
  return(DOOR_matrix)
}
