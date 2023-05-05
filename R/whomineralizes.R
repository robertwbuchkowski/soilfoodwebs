#' Direct and indirect contributions to mineralizations
#'
#' @param usin The community in which we want to calculate mineralization rates.
#' @param selected A vector of names for which you want to calculate the direct and indirect effects. Default NULL means all of them. Useful for excluding nodes whose removal breaks the community (i.e., basal nodes)
#' @return A table of node effects on mineralization rates.
#' @details
#' The results are labeled as follows with direct contributions calculated from the full food web and indirect contributions calculated from the food web without that node. Indirect contributions do not include the direct contribution (i.e., it is subtracted).
#'
#'\describe{
#'   \item{DirectC}{The direct contribution to carbon mineralization.}
#'   \item{DirectN}{The direct contribution to nitrogen mineralization.}
#'   \item{IndirectC}{The indirect contribution to carbon mineralization.}
#'   \item{IndirectN}{The indirect contribution to nitrogen mineralization.}
#' }
#' The indirect contributions are calculated as the total mineralization of the community with the trophic species minus the trophic species direct mineralization minus the total mineralization without the trophic species all divided by the total mineralizaiton with the trophic species.
#'
#' @examples
#' # Basic example for the introductory community:
#' whomineralizes(intro_comm)
#' @export
whomineralizes <- function(usin, selected = NULL){
  usin = checkcomm(usin, verbose = FALSE) # Check the community for errors
  Nnodes = dim(usin$imat)[1] # Get the number of nodes
  Nnames = usin$prop$ID # Get the names

  # Select only the chosen nodes:
  if(!is.null(selected)){
    if(!all(selected %in% Nnames)) stop("All selected nodes must be present in the community.")
    Nnames = selected
  }

  res1 <- comana(usin) # Calculate the C and N fluxes

  if(any(unname(usin$prop$ID) != names(res1$usin$prop$ID))){
    stop("Sorting of trophic levels not matching up")
  }

  # Calculate the direct and indirect mineralization rates
  output = data.frame(
    ID = unname(res1$usin$prop$ID), # Name of the node
    DirectC = res1$Cmin/sum(res1$Cmin), # Direct C mineralization--the amount that the species respires over the total
    DirectN = rowSums(res1$Nminmat)/sum(res1$Nminmat), # Direct N mineralization--the amount excreted divided by the total
    IndirectC = NA, # Save space for indirect
    IndirectN = NA) # Save space for indirect
  rownames(output) = NULL # remove rownames

  output = output[(output$ID %in% Nnames),]

  # Calculate the indirect effects by removing a single node from the community and calculating the difference
  for(rmnode in Nnames){
    usinmod = removenodes(usin, rmnode) # Remove a code
    res2 = comana(usinmod) # Calculate the new fluxes

    # Indirect carbon flux: accounts for change caused by removing the node less the node direct effect.
    output[output$ID == rmnode, "IndirectC"] =
      (sum(res1$Cmin) - res1$Cmin[rmnode] - sum(res2$Cmin))/sum(res1$Cmin)

    # Indirect nitrogen flux: accounts for change caused by removing the node less the node direct effect.
    output[output$ID == rmnode, "IndirectN"] =
      (sum(res1$Nminmat)- res1$Nmin[rmnode] - sum(res2$Nminmat))/sum(res1$Nminmat)
  }
  return(output)
}

