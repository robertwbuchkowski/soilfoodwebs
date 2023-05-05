#' Parameter uncertainty returns community with new parameters drawn from a distribution of choice
#'
#' @param usin The community in which we want to calculate mineralization rates.
#' @param parameters A vector of parameters you want to vary.
#' @param replacetiny A number. All parameter draws less than this value are replaced with it to avoid numerical errors in the calculations. Set to zero if you want all values to be left as drawn. Default is 0.000001.
#' @param distribution A single string or matrix for the distribution from which to draw the parameters. If it is a matrix is has rownames of web nodes matching usin and column names matching parameters. The acceptable options are gamma, normal, uniform.
#' @param errormeasure A single value or matrix following the format of distribution recording the error. Value depends on errortype.
#' @param errortype A single value or matrix following the format of distribution recording the error type. This can be "CV" for coefficient of variation, "Variance" for the variance, and "Min" for the minimum value. The latter can only be used when the distribution is uniform.
#' @param fcntorun The function you want to run on the resulting communities and the result you want to return. Current options are comana, whomineralizes, CNsim, decompexpt. You can also include any of the outputs of comana or decompexpt to automatically subset the results to the vector of interest. For example, Cmin only returns carbon mineralization.
#' @param returnprops Boolean. Do you want to return the communities with parameter values or just the results of the function? Only used if returnresults is TRUE.
#' @param returnresults Boolean. Do you want to return the results of the function? If this is FALSE, the fcntorun is ignored and a list of communities with parameter draws is returned.
#' @param replicates The number of replicate communities you want to create and analyze.
#' @param rejectnegconsump Boolean. Should the draws reject communities with negative consumption rates?
#' @param correctstoich Boolean. Do you want to correct the stoichiometry of the community before running the fcntorun? This does NOT correct the community stoichiometry returned in communitylist, so the user can see the original result without the correction applied.
#' @param verbose Boolean. Do you want warning messages about the functionality?
#' @return A list of the results. See details.
#' @details
#' The results are always in a list. If returnprops = T, then the top lay is a list of length 2 with resultslist and communitylist attributes, otherwise only resultslist is returned. The communitylist has the communities with parameter draws in order. The resultslist has the results of the function indicated in fcntorun.
#'
#' @examples
#' # Basic example for the introductory community:
#' parameter_uncertainty(intro_comm)
#' @export
parameter_uncertainty <- function(usin, parameters = c("B"), replacetiny = 0.000001, distribution = "gamma", errormeasure = 0.2, errortype = "CV", fcntorun = "comana", replicates = 100, returnprops = FALSE, returnresults = TRUE, rejectnegconsump = TRUE, correctstoich = TRUE, verbose = TRUE){

  # Confirm inputs are correct:
  if(!all(errortype %in% c("CV", "Variance", "Min"))) stop("errortype must be either CV or Variance or Min")

  if(!all(distribution %in% c("gamma", "normal", "uniform"))) stop("distribution must be either gamma or normal or uniform")

  if(!(fcntorun %in% c("comana", "whomineralizes", "CNsim", "decompexpt", "Cmin", "Nmin", "consumption", "Nminmat", "fmat", "Nfmat", "basedecomp", "decompeffects"))) stop("fcntorun must be a recognized one")

  # Check and finalize community:
  usin = checkcomm(usin)
  Nnodes = dim(usin$imat)[1]
  Nnames = usin$prop$ID

  # Turn single values into matrices that are compatible with the community:
  if(length(distribution) == 1){
    distribution = matrix(distribution, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  if(length(errormeasure) == 1){
    errormeasure = matrix(errormeasure, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  if(length(errortype) == 1){
    errortype = matrix(errortype, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  # Get the parameter mean values and set their names correctly
  meanvalues = as.matrix(usin$prop[,parameters])
  rownames(meanvalues) = usin$prop$ID
  colnames(meanvalues) = parameters

  # Check to make sure matrices are the right size:
  if(!all(
    c(dim(distribution) == c(Nnodes, length(parameters)), dim(errormeasure) == c(Nnodes, length(parameters)),dim(errortype) == c(Nnodes, length(parameters)))
  )) stop("Distribution, errormeasure, and errortype all need to be length 1 or have the dimensions number of nodes x length of parameters.")

  # Build two empty lists for the communities and the results of the analysis
  communitylist <- vector(mode = "list", length = replicates)

  resultslist <- vector(mode = "list", length = replicates)

  # Set warning flag
  if(verbose){
    onlyonce = c(0,0)
  }else{
    onlyonce = c(1,1)
  }

  # Loop over the replicates and draw parameters for each one:
  for(i in 1:replicates){
    # Create an internal copy of the community to modify:
    usintemp = usin

    # Start a loop to decide whether draw should be accepted
    accept = 0
    while(accept == 0){
      # Draw each parameter for each node:
      for(curparam in parameters){
        for(curID in Nnames){

          # First make sure we aren't modifying detritus a or p parameters.
          if(curparam %in% c("a", "p", "d") &
             usintemp$prop[usintemp$prop$ID == curID, "isDetritus"] == 1){
            if(onlyonce[1] ==0){
              warning("Detritus always has d = 0, a = 1, and p = 1, so it is skipped in the paramter draws.")
              onlyonce[1] = 1
            }
          }else{
            # Now draw the parameters for this case:
            curmean = meanvalues[curID, curparam]
            curdistribution = distribution[curID, curparam]
            curerrormeasure = errormeasure[curID, curparam]

            # Deal with uniform first:
            if(errortype[curID, curparam] == "Min"){
              if(curdistribution != "uniform") {stop("The errormeasure must be Min if using uniform distribution.")
              }else{
                usintemp$prop[usintemp$prop$ID == curID, curparam] = stats::runif(1, min = curerrormeasure, max = curmean)

                if(onlyonce[2] == 0){
                  warning("Uniform distribution uses the errormeasure as the minimum value and the parameter value in the community as the maximum.")
                  onlyonce[2] = 1
                }
              }
            }

            if(errortype[curID, curparam] == "CV"){
              curerrormeasure = (curmean*curerrormeasure)^2
            }

            if(curdistribution == "normal"){
              usintemp$prop[usintemp$prop$ID == curID, curparam] = stats::rnorm(1, mean = curmean, sd = sqrt(curerrormeasure))
            }

            if(curdistribution == "gamma"){
              usintemp$prop[usintemp$prop$ID == curID, curparam] = stats::rgamma(1, scale = curerrormeasure/curmean, shape = curmean^2/curerrormeasure)
            }
          }
        }
      } # Done drawing the parameter values

      if(replacetiny!=0){
        # Replace the tiny values drawn to avoid numerical errors.
        usintemp$prop[
          which(usintemp$prop[,parameters] > 0 & usintemp$prop[,parameters] < replacetiny),
          parameters] = replacetiny

      }

      # Now decide whether draw should be accepted
      if(rejectnegconsump){
        if(
          all(
            comana(usintemp)$consumption[TLcheddar(usintemp$imat) != 1] >= 0
            )
          ) accept = 1
      }else{
        accept = 1
      }
    }
    communitylist[[i]] = usintemp

    # Save the results of the analysis. This is a series of if statements to ask which of the functions should be run.
    if(returnresults){

      # Correct the stoichiometry if requested
      if(correctstoich){
        usintemp = corrstoich(usintemp)
      }

      if(fcntorun =="consumption"){
        resultslist[[i]] = comana(usintemp)$consumption
      }

      if(fcntorun =="consumption"){
        resultslist[[i]] = comana(usintemp)$consumption
      }

      if(fcntorun =="Cmin"){
        resultslist[[i]] = comana(usintemp)$Cmin
      }

      if(fcntorun =="Nmin"){
        resultslist[[i]] = comana(usintemp)$Nmin
      }

      if(fcntorun =="Nminmat"){
        resultslist[[i]] = comana(usintemp)$Nminmat
      }

      if(fcntorun =="fmat"){
        resultslist[[i]] = comana(usintemp)$fmat
      }

      if(fcntorun =="Nfmat"){
        resultslist[[i]] = comana(usintemp)$Nfmat
      }

      if(fcntorun == "comana"){
        resultslist[[i]] = comana(usintemp)
      }

      if(fcntorun == "whomineralizes"){
        resultslist[[i]] = whomineralizes(usintemp)
      }

      if(fcntorun == "CNsim"){
        resultslist[[i]] = CNsim(usintemp)
      }

      if(fcntorun == "decompexpt"){
        resultslist[[i]] = decompexpt(usintemp)
      }

      if(fcntorun == "basedecomp"){
        resultslist[[i]] = decompexpt(usintemp)$basedecomp
      }

      if(fcntorun == "decompeffects"){
        resultslist[[i]] = decompexpt(usintemp)$decompeffects
      }
    }
  } # Close replicates loop

  if(returnresults){
    if(!returnprops){
      return(resultslist)
    }else{
      return(list(resultslist=resultslist,communitylist = communitylist))
    }
  }else{
    return(communitylist)
  }


}
