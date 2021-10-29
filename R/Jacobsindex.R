#' Calculate Jacob's index.
#'
#' @param usin The community for which to calculate Jacob's index
#' @return A data frame with columns for predator, prey, and the value of Jacob's index.
#' @details
#' Jacob's index is calculated for all predators with more than one prey item. A value of -1 indicates maximum avoidance, 0 indicates no preference, and 1 indicates maximum preference. This index is a more user friendly output than the values in imat, which can be difficult to interpret.
#' @examples
#' # Calculate the minimum values of s for the introduction community.
#' Jacobsindex(intro_comm)
#' @export
Jacobsindex = function(usin){
  usin = checkcomm(usin)

  imat = usin$imat
  prop = usin$prop
  Nnodes = dim(imat)[1]

  r = (imat *matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)/rowSums(imat * matrix(prop$B, nrow = Nnodes,ncol = Nnodes, byrow = T)))

  imat2 = imat
  imat2[imat2 > 0] = 1

  p = (imat2 *matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)/rowSums(imat2 * matrix(prop$B, nrow = Nnodes,ncol = Nnodes, byrow = T)))

  JI = (r-p)/(r+p-2*p*r)

  JI = data.frame(JI)
  JI$Predator = row.names(JI)
  rownames(JI) = NULL

  JI = stats::reshape(data=JI, idvar="Predator",
               varying = colnames(JI)[colnames(JI) != "Predator"],
               times = colnames(JI)[colnames(JI) != "Predator"],
               v.names = "Jacobs",
               timevar = "Prey",
               new.row.names = 1:1e6,
               direction="long")

  JI = subset(JI, !is.na(JI[,"Jacobs"]))

  return(JI)

}
