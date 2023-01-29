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

  # Check the community for errors
  usin = checkcomm(usin)

  imat = usin$imat # Save the feeding matrix
  prop = usin$prop # Same the property data
  Nnodes = dim(imat)[1] # Calculate the number of nodes

  # Calculate the components of Jacob's Index r and p.
  r = (imat *matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)/rowSums(imat * matrix(prop$B, nrow = Nnodes,ncol = Nnodes, byrow = T)))

  imat2 = imat
  imat2[imat2 > 0] = 1

  p = (imat2 *matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)/rowSums(imat2 * matrix(prop$B, nrow = Nnodes,ncol = Nnodes, byrow = T)))

  # Jacob's index formula
  JI = (r-p)/(r+p-2*p*r)

  # Format the output
  JI = data.frame(JI)
  JI$Predator = row.names(JI)
  rownames(JI) = NULL

  # Reshape the output into a cleaner table
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
