#' Parameter uncertainty returns community with new parameters drawn from a distribution of choice
#'
#' @param usin The community in which we want to calculate mineralization rates.
#' @param parameters A vector of parameters you want to vary.
#' @param replacetiny A number. All parameter draws less than this value are replaced with it to avoid numerical errors in the calculations. Set to zero if you want all values to be left as drawn. Default is 0.000001.
#' @param distribution A single string or matrix for the distribution from which to draw the parameters. If it is a matrix is has rownames of web nodes matching usin and column names matching parameters. The acceptable options are gamma, normal, uniform.
#' @param errormeasure A single value or matrix following the format of distribution recording the error. Value depends on errortype.
#' @param errortype A single value or matrix following the format of distribution recording the error type. This can be "CV" for coefficient of variation, "Variance" for the variance, and "Min" for the minimum value. The latter can only be used when the distribution is uniform.
#' @param fcntorun The function you want to run on the resulting communities. Current options are comana, whomineralizes, and CNsim.
#' @param returnprops Boolean. Do you want to return the communities with parameter values or just the results of the function? Only used if returnresults is TRUE.
#' @param returnresults Boolean. Do you want to return the results of the function? If this is FALSE, the fcntorun is ignored and a list of communities with parameter draws is returned.
#' @param replicates The number of replicate communities you want to create and analyze.
#' @return A list of the results. See details.
#' @details
#' The results are always in a list. If returnprops = T, then the top lay is a list of length 2 with resultslist and communitylist attributes, otherwise only resultslist is returned. The communitylist has the communities with parameter draws in order. The resultslist has the results of the function indicated in fcntorun.
#'
#' @examples
#' # Basic example for the introductory community:
#' parameter_uncertainty(intro_comm)
#' @export
parameter_uncertainty <- function(usin, parameters = c("B"), replacetiny = 0.000001, distribution = "gamma", errormeasure = 0.2, errortype = "CV", fcntorun = "comana", replicates = 100, returnprops = F, returnresults = T){

  # Confirm inputs are correct:
  if(!all(errortype %in% c("CV", "Variance", "Min"))) stop("errortype must be either CV or Variance or Min")

  if(!all(distribution %in% c("gamma", "normal", "uniform"))) stop("distribution must be either gamma or normal or uniform")

  if(!(fcntorun %in% c("comana", "whomineralizes", "CNsim"))) stop("fcntorun must be either comana or whomineralizes or CNsim")

  # Check and finalize community:
  usin = checkcomm(usin)
  Nnodes = dim(usin$imat)[1]
  Nnames = usin$prop$ID

  # Turn single selections into matrices:
  if(length(distribution) == 1){
    distribution = matrix(distribution, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  if(length(errormeasure) == 1){
    errormeasure = matrix(errormeasure, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  if(length(errortype) == 1){
    errortype = matrix(errortype, nrow = Nnodes, ncol = length(parameters), dimnames = list(Nnames, parameters))
  }

  meanvalues = as.matrix(usin$prop[,parameters])
  rownames(meanvalues) = usin$prop$ID
  colnames(meanvalues) = parameters

  # Check to make sure matrices are the right size:
  if(!all(
    c(dim(distribution) == c(Nnodes, length(parameters)), dim(errormeasure) == c(Nnodes, length(parameters)),dim(errortype) == c(Nnodes, length(parameters)))
  )) stop("Distribution, errormeasure, and errortype all need to be length 1 or have the dimensions number of nodes x length of parameters.")

  communitylist <- vector(mode = "list", length = replicates)

  resultslist <- vector(mode = "list", length = replicates)

  # Set warning flag
  onlyonce = c(0,0)

  # Loop over the replicates and draw parameters for each one:
  for(i in 1:replicates){
    # Create an internal copy of the community to modify:
    usintemp = usin

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
      usintemp$prop[,parameters][which(usintemp$prop[,parameters] > 0 & usintemp$prop[,parameters] < replacetiny)] = replacetiny

    }

    communitylist[[i]] = usintemp

    if(returnresults){
      if(fcntorun == "comana"){
        resultslist[[i]] = comana(usintemp)
      }

      if(fcntorun == "whomineralizes"){
        resultslist[[i]] = whomineralizes(usintemp)
      }

      if(fcntorun == "CNsim"){
        resultslist[[i]] = CNsim(usintemp)
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
