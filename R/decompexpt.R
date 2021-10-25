#' Decomposition rates and effect of individual organisms
#'
#' @param usin The community in which we want to calculate decomposition rates.
#' @param overtime Do you want to return a decomposition trajectory overtime with and without all species? If so, set this to the number of time steps you want to see.
#' @return A list. The first item is the basic decomposition constants. The second is a table of decomposition constants for each detritus pool. The third only runs if asked for by `overtime` and is the predicted trajectory.
#' @details
#' The output list indicates the baseline decomposition rate (k; basedecomp) and then a table showing the effect on decomposition rates when the direct and indirect effects of individual trophic species are removed (decompeffects). Removing direct effects means setting consumption of that species to zero, while indirect effects occur when the community is recalculated without that species before calculating the decompostion rate (k). A positive number means the species increases decomposition. The column ID contains the node having the ffect, while DetritusID indicates the detritus pool whose rate is being reported.
#' Direct effects are calculated as the decomposition rate in the full community minus decomposition rate without the species consumption of detritus, divided by the decomposition rate in the full community.
#' Indirect effects are calculated as the difference in decomposition rate of the full community minus decomposition rate without the species and minus, divided by the decomposition rate with the species and minus the direct effect. A negative value here means that removing the species reduces decomposition more than any direct consumption effects it has on detritus.
#' If `overtime` is a number, a table of the proportion of detritus over time is returned with and without all species. Column "Original" is the original community, while all successive columns show the removal of each trophic species.
#'
#' @examples
#' # Basic example for the introductory community:
#' decompexpt(intro_comm)
#' @export
decompexpt <- function(usin, overtime = NA){
  usin = checkcomm(usin)
  Nnodes = dim(usin$imat)[1]
  IDdetritus = which(usin$prop$isDetritus > 0)
  Nnames = usin$prop$ID[which(usin$prop$isDetritus == 0)] # List of species that aren't detritus

  res1 <- comana(usin)

  basedecomp <- (colSums(res1$fmat)/usin$prop$B)[IDdetritus]

  # Create a list for the results:
  results = vector("list", length(Nnames))
  overtimesave = vector("list", length(Nnames))

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

    overtimesave[[which(Nnames == rmnode)]] = Indirect

    # The formula occurs because it really is (k_base - k_indirect - (k_base - k_direct))/k_base = (k_direct - k_indirect)/k_base
    Indirect2 = (Direct - Indirect)/basedecomp

    results[[which(Nnames == rmnode)]] = data.frame(ID = rmnode, DetritusID = usin$prop$ID[IDdetritus], Direct = Direct2, Indirect = Indirect2)
  }

  # Get rid of rownames
  results = do.call("rbind",results)
  rownames(results) = NULL

  # Run the simulations using the overtime saved from above
  if(!is.na(overtime)){

    # Make sure it is an integer.
    if(is.integer(overtime) & overtime >0) stop("overtime must be either NA or a positive integer.")

    overtimeout = vector("list", length(overtimesave[[1]]))
    names(overtimeout) = names(overtimesave[[1]])

    for(i in 1:length(overtimeout)){

      # Create vector for the results
      overtimesave2 = vector("list", length(Nnames))

      for(j in 1:length(overtimesave)){
        overtimesave2[[j]] = exp(-overtimesave[[j]][i]*1:overtime)
      }
      overtimeout[[i]] = do.call("cbind",overtimesave2)
      colnames(overtimeout[[i]]) = Nnames
      overtimeout[[i]] = data.frame(cbind(
        Day = 1:overtime,
        Original = exp(-basedecomp[i]*1:overtime),
        overtimeout[[i]]))

    }

    return(list(basedecomp = basedecomp, decompeffects = results,overtime = overtimeout))
  }

  return(list(basedecomp = basedecomp, decompeffects = results))
}

