#' A function to correct stoichiometry dynamically
#'
#' @param biomasses A vector of biomasses.
#' @param cij The consumption matrix from Cijfcn
#' @param CN A vector of C:N ratios.
#' @param p A vector of production efficiency.
#' @param a A vector of assimilation efficiency.
#' @param canIMM A Boolean vector of whether the nodes can immobilize nitrogen.
#' @param dietlimits The diet limits matrix for the stoichiometry correction (proportion of diet)?
#' @param diet_correct Boolean: Does the organism correct it's diet?
#' @param Conly Boolean: Is the model meant for carbon only?
#' @param Immobilizationlimit This is the limit of the amount of nitrogen the food web can immobilize nitrogen (NOT PLANTS). This will impact the calculations of inorganic nitrogen dynamics.
#' @return Returns the consumption rates (FMAT) and production efficiencies (p).
#' @details
#' This function takes inputs from the ODE and outputs corrected consumption rates.
#' The key difference from 'corrstoich' is that the prey DO NOT correct their feeding rates to compensate for higher consumption from the predators, so the system can leave equilibrium if a diet shift occurs.
correction_function <- function(biomasses, cij, CN, p, a, canIMM, dietlimits,diet_correct = T, Conly = F, Immobilizationlimit = Inf){

  Nnodes = length(biomasses)

  Bpred = matrix(biomasses, ncol = Nnodes, nrow = Nnodes) # A matrix of predators
  Bprey = matrix(biomasses, ncol = Nnodes, nrow = Nnodes, byrow = T) # A matirx of prey

  # **Non-equilibrium** consumption matrix
  FMAT = cij*Bpred*Bprey
  # AIJ function calculation
  Nnodes = dim(FMAT)[1]
  FMAT2 = FMAT
  FMAT2[FMAT2>0] = 1
  aij = (matrix(1/CN, nrow = Nnodes, ncol = Nnodes, byrow = T) - matrix(p/CN, nrow=Nnodes, ncol = Nnodes))*matrix(a, nrow=Nnodes, ncol = Nnodes)*FMAT2

  if(!Conly){
    if(diet_correct){
      species = which(rowSums(aij*FMAT) < -1e-10  & # Species must have negative Nmin rate that is beyond a tiny value caused by a rounding error in the calaculation.
                        apply(aij,1, max) > 0 & # Species must have a food item that gives them a positive nitrogen balance
                        apply(FMAT2 > 0, 1, sum) > 1 & # Species must have more than one food item
                        !canIMM == 1 # The species is not flagged as a species that can immobilize N
      )

      for(sp in species){

        food = FMAT[sp,] > 0
        ai = aij[sp,food]
        biomass = biomasses[food]
        FT = rowSums(FMAT)[sp]
        limits = dietlimits[sp,food]

        # Try to see whether there is a solution where the diet can be optimized
        sol1 <- NULL
        try(
          sol1 <- quadprog::solve.QP(Dmat = diag(length(ai)), # We need the squared terms for each diet, use the identity matrix to get them f^T I f
                                     dvec = biomass/sum(biomass), # The diet is as close to the relative abundance as possible
                                     Amat = t(rbind(rep(1, length(ai)),ai*FT,-diag(nrow = length(ai)),diag(nrow = length(ai)))), # proportions sum to 1, Nmin is zero, and none of the limits are exceeded
                                     bvec = c(1, 0, -limits, rep(0,length(ai))), # proportions sum to 1, Nmin is zero, and none of the limits are exceeded, and all values greater than zero
                                     meq = 1 # first position in Amat is an equality constraint)
          ),
          silent = T)

        # If there is no solution, then run a linear program to find the closest diet within the limits
        if(is.null(sol1)){
          sol1 <- lpSolve::lp(direction = "max",
                              objective.in = c(ai*FT),
                              const.mat = rbind(rep(1, length(ai)),-diag(nrow = length(ai)),diag(nrow = length(ai))),
                              const.dir = c("=", rep(">=", 2*length(ai))),
                              const.rhs = c(1, -limits, rep(0, length(ai)))
          )
        }

        # Confirm that solution is positive. Sometimes the solution produces a tiny negative number because of a rounding error (e.g. 1e-18). Set this to zero.
        solcheck = sol1$solution
        if(any(solcheck < -1e-10)) stop("Solution is returning a large negative value. There has been an error in the optimization.")
        if(any(solcheck < 0)){
          solcheck[solcheck < 0] = 0
        }

        FMAT[sp,food] = solcheck*FT
      }
    }

    # Modify the production efficiency:
    CNj = CN
    CNi = matrix(CN, ncol = Nnodes, nrow = Nnodes, byrow = T)

    # If Immobilization limit is infinity then you don't need to correct the production efficiency of immobilizing species
    if(is.infinite(Immobilizationlimit)){
      Nmin = 0 # Nmin is set to zero for this case to make term disappear
      temp1 = CNj*(rowSums(FMAT/CNi)-Nmin/a)/rowSums(FMAT)

      temp1[is.na(temp1)] = 1

      tochange = p > temp1 & canIMM != 1

      p[tochange] = temp1[tochange]
    }else{
      # Start with Nmin floor for each trophic level. Is zero if canIMM = 0 and is a portion of Immobilizationlimit otherwise.
      Nminj = - Immobilizationlimit*canIMM*biomasses/sum(canIMM*biomasses)

      temp1 = CNj*(rowSums(FMAT/CNi)-Nminj/a)/rowSums(FMAT)

      temp1[is.na(temp1)] = 1

      tochange = p > temp1

      p[tochange] = temp1[tochange]
    }
  }

  return(list(FMAT= FMAT, p = p))
}
