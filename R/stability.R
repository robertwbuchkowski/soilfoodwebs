#' Calculates the stability of the food web
#'
#' @param usin The community to calculate stability.
#' @param correctstoich Boolean: Should stability be calculated after the stoichiometry is corrected with corrstoich?
#' @param forceProd Boolean: Should production efficiency be the only way to correct stoichiometry?
#' @param smin The value of smin in the Moorecobian.
#' @param method One of two options: Jacobian or Moorecobian. The later uses the value of smin and adds density-dependence to the calculation at the strength of smin.
#' @return Returns the Jacobian (or Moorecobian), eigenvalues, and the maximum eignvalue as a list.
#' @examples
#' # Basic stability calculation
#' stability(intro_comm)
#' @seealso \code{\link{calc_smin}} for the full use of the Moorecobian and \code{\link{stability2}} for the estimation of stability using the functions from \code{\link[rootSolve]{jacobian.full}}
#' @export
stability <- function(usin, correctstoich = TRUE, forceProd = FALSE, smin = 1, method = "Jacobian"){ # Function only requires the community inputs

  # Confirm that the type of stability test is correct
  if(!method %in% c("Jacobian", "Moorecobian")) stop("method must be either Jacobian or Moorecobian")

  # Correc the stoichiometry of the community if set
  if(correctstoich){
    usin = corrstoich(usin, forceProd = forceProd)
  }

  imat = usin$imat # Get the feeding matrix
  prop = usin$prop # Get the node properties

  # Warn if more than 1 plant is included that the function is not set up for this yet
  try(if(sum(prop$isPlant) > 1) warning("Function is only set up for 1 plant pool at the current time"))

  # Rescale DetritusRecycling to sum to 1
  if(sum(prop$DetritusRecycling) != 0 ){
    prop$DetritusRecycling = prop$DetritusRecycling/sum(prop$DetritusRecycling)
  }

  Nnodes = dim(imat)[1] # Number of nodes

  Bpred = matrix(prop$B, ncol = Nnodes, nrow = Nnodes) # A matrix of predators
  Bprey = matrix(prop$B, ncol = Nnodes, nrow = Nnodes, byrow = TRUE) # A matirx of prey
  fmat = comana(usin, shuffleTL = FALSE)$fmat # Get the consumption matrix (units = gC/ time)
  cij = fmat/(Bpred*Bprey) # Get the consumption rate matrix (units 1/ (gC * time))

  # List which pools are detritus
  detritusPOS = which(prop$isDetritus >0)

  D1 = rep(NA, dim(cij)[1]) # Create a vector that will be the diagonal of the Jacobian. This is the partial derivative of each trophic species equation w.r.t. the bioamss of that species
  for(ig in 1: dim(cij)[1]){
    if(method == "Jacobian"){ # | ig %in% detritusPOS
      D1[ig] = sum(prop$a[ig]*prop$p[ig]*cij[ig,]*Bprey[ig,]- # Deriative of consumption term of the focal species
                     cij[,ig]*Bpred[,ig]) - # Deritative of terms where the focal species is a prey item
        prop$d[ig] # Deriative of the death rate term
    }else{
      if(method == "Moorecobian"){
        D1[ig] = - prop$d[ig]*smin # Deriative of the death rate term, times smin if this is being adjusted. It is default set to 1, which has no effect on the Jacobian
      }
    }

  }

  # Need to ID the basal trophic species and fix the terms if they are plants or detritus

  # If plant, add the plant growth rate so the growth rate is enough to balance losses
  if(any(prop$isPlant>0)){
    D1[prop$isPlant>0] = D1[prop$isPlant>0] +
      (comana(usin)$consumption/prop$B)[prop$isPlant>0]
  }

  # If detritus, the detritus equation has terms for the biomass recycling of all trophic species that eat detritus and recycle the unassimilated biomass back to that pool.

  if(any(prop$isDetritus>0) &
     method != "Moorecobian"){

    for(DPOS in detritusPOS){
      D1[DPOS] = D1[DPOS]  +
        prop$DetritusRecycling[DPOS]*sum((1-prop$a)*cij[,DPOS]*prop$B)
    }


  }

  # Create a Jacobian of the same dimensions as the consumption matrix
  J = cij

  for(ii in 1:dim(J)[1]){ # Focus equation
    for(jj in 1:dim(J)[1]){ # derivative w.r.t.

      if(cij[ii,jj] >0 & cij[jj,ii] ==0) J[ii,jj] = prop$a[ii]*prop$p[ii]*cij[ii,jj]*prop$B[ii] # The derivative is based on consumption if the focal species is eaten by the species whose biomass is in the derivative

      if(cij[ii,jj] ==0 & cij[jj,ii] > 0) J[ii,jj] = -cij[jj,ii]*prop$B[ii] # The derivative is based on loss to predation is the focal equation is the prey and the derivative is take w.r.t. the predator biomass

      if(cij[ii,jj] > 0 & cij[jj,ii] > 0) J[ii,jj] = prop$a[ii]*prop$p[ii]*cij[ii,jj]*prop$B[ii] - cij[jj,ii]*prop$B[ii] # Both terms are present for species that feed on each other.
    }
  }

  diag(J) = D1 # Replace the diagonal with the one calculated above

  # Fix the detritus calculations for the derivatives associated with detritus
  if(any(prop$isDetritus >0)){

    for(DPOS in detritusPOS){

      all_other_pools = 1:Nnodes
      all_other_pools = all_other_pools[all_other_pools != DPOS]
      J[DPOS,all_other_pools] =
        J[DPOS,all_other_pools] +

        prop$DetritusRecycling[DPOS]*(
          prop$d[all_other_pools] + # Death rates

            rowSums(cij[all_other_pools,]*matrix(prop$B, nrow = length(all_other_pools), ncol = Nnodes, byrow=T)*matrix((1-prop$a[all_other_pools]), nrow = length(all_other_pools), ncol = Nnodes, byrow=FALSE)) + # Add poop from when focal pool is predator

            colSums(cij[,all_other_pools]*matrix(prop$B, nrow = Nnodes, ncol = length(all_other_pools), byrow=FALSE)*matrix((1-prop$a), nrow = Nnodes, ncol = length(all_other_pools), byrow=FALSE))# add poop from where it is prey
        )

    }
  }

  return(list(J = J, # Return the Jacobian
              eigen = base::eigen(J), # Return the eigenvalues and vectors
              rmax = max(Re(base::eigen(J)$values)))) # Return the maximum eigenvalue. The system is stable if this is negative
}
