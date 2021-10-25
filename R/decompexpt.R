#' Decomposition rates and effect of individual organisms
#'
#' @param usin The community in which we want to calculate decomposition rates.
#' @return A table of decomposition constants for each detritus pool.
#' @details
#' The output list indicates the baseline decomposition rate (k; basedecomp) and then a table showing the effect on decomposition rates when the direct and indirect effects of individual trophic species are removed (decompeffects). Removing direct effects means setting consumption of that species to zero, while indirect effects occur when the community is recalculated without that species before calculating the decompostion rate (k). A positive number means the species increases decomposition. The column ID contains the node having the ffect, while DetritusID indicates the detritus pool whose rate is being reported.
#' Direct effects are calculated as the decomposition rate in the full community minus decomposition rate without the species consumption of detritus, divided by the decomposition rate in the full community.
#' Indirect effects are calculated as the difference in decomposition rate of the full community minus decomposition rate without the species and minus, divided by the decomposition rate with the species and minus the direct effect. A negative value here means that removing the species reduces decomposition more than any direct consumption effects it has on detritus.
#'
#' @examples
#' # Basic example for the introductory community:
#' decompexpt(intro_comm)
#' @export
decompexpt <- function(usin){
  usin = checkcomm(usin)
  Nnodes = dim(usin$imat)[1]
  IDdetritus = which(usin$prop$isDetritus > 0)
  Nnames = usin$prop$ID[which(usin$prop$isDetritus == 0)] # List of species that aren't detritus

  res1 <- comana(usin)

  basedecomp <- (colSums(res1$fmat)/usin$prop$B)[IDdetritus]

  # Create a list for the results:
  results = vector("list", length(Nnames))


  for(rmnode in Nnames){

    # Direct effects by removing consumption of each species without changing the matrix:
    FMAT = res1$fmat
    FMAT[rmnode,] = 0 # Set consumption of that species to zero

    # The decomposition rate without the direct effects
    Direct = (colSums(FMAT)/usin$prop$B)[IDdetritus]

    # The direct effect (removing species flips the order of the numerator) from what you expect
    Direct2 = (basedecomp - Direct)/basedecomp

    # Indirect effects of removing the individual species.
    usinmod = removenodes(usin, rmnode)
    res2 = comana(usinmod)

    # The decomposition rate without any effects
    Indirect = (colSums(res2$fmat)/usinmod$prop$B)[which(usinmod$prop$isDetritus > 0)]

    # The formula occurs because it really is (k_base - k_indirect - (k_base - k_direct))/k_base = (k_direct - k_indirect)/k_base
    Indirect = (Direct - Indirect)/basedecomp

    results[[which(Nnames == rmnode)]] = data.frame(ID = rmnode, DetritusID = usin$prop$ID[IDdetritus], Direct = Direct2, Indirect = Indirect)
  }

  # Get rid of rownames
  results = do.call("rbind",results)
  rownames(results) = NULL

  return(list(basedecomp = basedecomp, decompeffects = results))
}

