#' A function to fix production efficiency.
#'
#' @param usin The input community to fix the production efficiency.
#' @param immobilizationlimit Set a limit for the maximum rate of immobilization of inorganic nitrogen for the food web per time step.
#' @return The modified community with new production efficiency.
productionadj = function(usin,immobilizationlimit = Inf
){
  temp0 = comana(usin, shuffleTL = FALSE)
  FMAT = temp0$fmat
  NMIN = sum(temp0$Nminmat)
  Nnodes = dim(FMAT)[1]
  CNj = usin$prop$CN
  CNi = matrix(usin$prop$CN, ncol = Nnodes, nrow = Nnodes, byrow = TRUE)
  aj = usin$prop$a

  if(-NMIN > immobilizationlimit){

    # First change the cases where Nmin must be greater than zero:
    temp1 = CNj*rowSums(FMAT/CNi)/rowSums(FMAT)
    temp1[is.na(temp1)] = 1

    tochange = usin$prop$p > temp1 & usin$prop$canIMM != 1

    usin$prop$p[tochange] = temp1[tochange]


    # Then loop through the cases where nodes can immobilize (canIMM = 1) and fix their production efficiencies:

    imm_props = usin$prop$canIMM*usin$prop$B/sum(usin$prop$canIMM*usin$prop$B) # calculate the biomass proportion for distributing the immobilization quota
    imm_limits = imm_props*immobilizationlimit*-1 # Assign the immobilization limits to each pool

    # If Nmin is limited and not zero, we need to include consumption rate modifications into the calculation.
    # This means the calculations of the new production efficiency values need to iterate down the food web to account for the corrections to consumption rate (Mj).

    # Identify the species that can immobilize and sort by TL:
    species = which(
      usin$prop$canIMM == 1,
      arr.ind = TRUE)

    TLs = TLcheddar(usin$imat)
    species = species[order(TLs[species])]

    for(sp in species){
      FMAT = comana(usin, shuffleTL = FALSE)$fmat # calculate fmat
      Mj = sum(FMAT[,sp])
      Bj = usin$prop$B[sp]
      dj = usin$prop$d[sp]
      CNj = usin$prop$CN[sp]
      # Sum of prey N:C ratio
      E_inv_CNi = sum(1/(usin$prop$CN[FMAT[sp,] > 0]))
      Nminj = imm_limits[sp]

      # Calculate the new production efficiency using the formula derived from subsituting the expression for Fij into the Nmin equation and solving for p. This means production efficiency is corrected to match the new consumption rate that is necessary to keep the equilibrium
      usin$prop$p[sp] = E_inv_CNi*(dj*Bj - Mj)/(Nminj + (dj*Bj - Mj)/CNj)
    }
  }else{
    temp1 = CNj*rowSums(FMAT/CNi)/rowSums(FMAT)
    temp1[is.na(temp1)] = 1

    tochange = usin$prop$p > temp1 & usin$prop$canIMM != 1

    usin$prop$p[tochange] = temp1[tochange]

  }
  return(usin)
}

#' A function to correct the diet of trophic species.
#'
#' @param usin The input community in which to fix the diets.
#' @param dietlimits # A matrix the same size as imat that gives the diet limits as a proportion of the total diet. All values must be between 0 and 1. Leaving it as NA sets the limits of all diet items to 1.
#' @return The modified community with new diet preferences.
correct_diet <- function(usin,dietlimits = c(NA)){

  # Setting up and verifying diet limits
  if(any(is.na(dietlimits))){
    dietlimits = usin$imat
    dietlimits[dietlimits > 0] = 1
  }else{
    if(!all(dim(dietlimits) == dim(usin$imat))) stop("dietlimits must have the same dimensions as imat")
    if(any(dietlimits > 1) | any(dietlimits < 0)) stop("dietlimits must be a proportion of the diet between 0 and 1")
  }

  #Identify the species that need correction
  AIJ = Aijfcn(usin)
  species = which(
    rowSums(comana(usin, shuffleTL = FALSE)$Nminmat) < 0 & # Species must have negative Nmin rate
      apply(AIJ,1, max) > 0 & # Species must have a food item that gives them a positive nitrogen balance
      apply(usin$imat > 0, 1, sum) > 1 & # Species must have more than one food item
      !usin$prop$canIMM == 1 # The species is not flagged as a species that can immobilize N
  )

  for(sp in species){

    while(TRUE){
      food = usin$imat[sp,] > 0
      ai = AIJ[sp,food]
      biomass = usin$prop$B[food]
      FT = rowSums(comana(usin, shuffleTL = F)$fmat)[sp]
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

      if(any(solcheck < 1e-10)){
        warning(paste0("Diet correction removes an item from the diet of species ",sp,". May get strange model behavior. Check outputs for errors in diet proportions! This code just deletes the food item and distirbutes evenly across the other food items."))

        usin$imat[sp,usin$imat[sp,] > 0][which(solcheck < 1e-10)] = 0

      }else{
        break()
      }
    }

    FLIST = solcheck
    FLIST1 = FLIST/biomass
    FLIST1 = min(FLIST1[FLIST1 > 0])
    FLIST2 = (FLIST/biomass)/FLIST1

    usin$imat[sp,food] = FLIST2

    if(sum((comana(usin,shuffleTL = F)$fmat[sp,food]/FT - FLIST)^2) > 1e-10) stop("Check quadratic optimization. There is an issue.")

  }

  return(usin)

}

#' Correct stoichiometry
#'
#' @param usin # Community in which to correct stoichiometry.
#' @param forceProd Boolean: Should we force organisms to only change their production efficiency and not their diet?
#' @param dietlimits A matrix the same size as imat that gives the diet limits as a proportion of the total diet. All values must be between 0 and 1. Leaving it as NA sets the limits of all diet items to 1.
#' @param Immobilizationlimit Set a limit for the maximum rate of immobilization of inorganic nitrogen for the food web per time step.
#' @details
#' A function that corrects trophic species production efficiency and/or diet composition to balance carbon and nitrogen demands.
#' @return The modified community with new production efficiencies and diets.
#' @examples
#' corrstoich(intro_comm)
#' # To force the correction to modify production efficiency only.
#'
#' corrstoich(intro_comm, forceProd = TRUE)
#'
#' # Set limits on the composition of the animal diet
#' DL = intro_comm$imat # copy over the feeding matrix to use for the diet limits.
#' # Oribatid 1 can only have up to 20% of its diet from each fungal species
#' DL["Orib1",] = c(0,0,0,0.2,0.2,1,1)
#' # Oribatid 2 can only have up to 10% of its diet from fungi 1
#' DL["Orib2",] = c(0,0,0,0.1,1,1,1)
#'
#' # Run them with the limits:
#' corrstoich(intro_comm, dietlimits = DL)
#' @export
corrstoich <- function(usin,
                       forceProd = FALSE,
                       dietlimits = c(NA),
                       Immobilizationlimit = Inf
){
  usin = checkcomm(usin)

  if(forceProd){
    usin3 = productionadj(usin = usin, immobilizationlimit = Immobilizationlimit)
  }else{
    usin2 = correct_diet(usin = usin, dietlimits = dietlimits)

    usin3 = productionadj(usin = usin2, immobilizationlimit = Immobilizationlimit)
  }
  return(usin3)
}
