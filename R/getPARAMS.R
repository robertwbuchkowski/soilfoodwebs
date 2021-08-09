#' A function to get the parameters for a food web model.
#'
#' @param usin The community you want to simulate.
#' @param DIETLIMTS The diet limits matrix for the stoichiometry correction (proportion of diet)?
#' @param diet_correct Boolean: Does the organism correct it's diet?
#' @param Conly Boolean: Is the model meant for carbon only?
#' @param userdefinedinputs Do you want to input a user defined vector of input functions? If NA the input values that keep the system at equilibrium are calculated. If not, put in a vector of input rates for each node.
#' @param returnCNmin Boolean: Do you want to add Cmin and Nmin to the list of the results?
#' @param has_inorganic_nitrogen Boolean: Is there an inorganic nitrogen pool?
#' @param densitydependence Which nodes have density dependence? NA default means none of them do. Should be a vector of 0 (no DD) and 1 (DD) for each node.
#' @param inorganic_nitrogen_properties A list of state variables for the inorganic nitrogen pool (INN = inputs, q = per capita loss of N, eqmN = equilibrium N). Must include a value for two of the three variables and has the final one as NA.
#' @param verbose Boolean: Do you want extra warnings updates?
#' @param Immbolizationlimit This is the limit of the amount of nitrogen the food web can immobilize nitrogen (NOT PLANTS). This will impact the calculations of inorganic nitrogen dynamics.
#' @return A list with two elements: (1) a vector of parameters to run the model away from equilibrium and (2) a vector of starting biomasses.
#' @details
#' A function to get the parameters of a food web model for simulation purposes. It does not correct stoichiometry, so the user must do this beforehand if they want.
#' @export
getPARAMS <- function(usin,
                      DIETLIMTS = NA,
                      diet_correct = T,
                      Conly = F,
                      userdefinedinputs = NA,
                      returnCNmin = F,
                      has_inorganic_nitrogen = F,
                      densitydependence = NA,
                      inorganic_nitrogen_properties = list(INN = NA, q = NA, eqmN = NA),
                      verbose = T,
                      Immbolizationlimit = Inf){
  if(any(is.na(DIETLIMTS))){
    DIETLIMTS = usin$imat
    DIETLIMTS[DIETLIMTS > 0] = 1
  }


  if(Conly & has_inorganic_nitrogen) stop("Cannot have a carbon-only  model with inorganic nitrogen.")

  # If the model contains nitrogen, then correct the stoichiometry before drawing parameters:
  if(!Conly){
    usin = corrstoich(usin, forceProd = !diet_correct, dietlimits = DIETLIMTS, Immobilizationlimit = Immbolizationlimit)
  }

  Nnodes = dim(usin$imat)[1] # Number of nodes

  if(any(is.na(densitydependence))){
    densitydependence = rep(0, Nnodes)
  }else{
    stopifnot(length(densitydependence) == Nnodes)
  }

  cij = Cijfcn(usin) # Get the consumption rate matrix for the base parameters (units 1/ (gC * time))

  # This won't make corrections, since they were made above, but is a convenient way to get FMAT. Never correct the production efficiency for immobilization limit here, because we retain maximum production efficiency parameters in the full model.
  corred = correction_function(usin$prop$B, cij, usin$prop$CN, usin$prop$p, usin$prop$a, usin$prop$canIMM, dietlimits = DIETLIMTS, diet_correct = diet_correct, Conly = Conly, Immobilizationlimit = Inf)

  if(isFALSE(all.equal(corred$p, usin$prop$p))){
    warning("Some secondary correction happening. Check the code running correction_function")
  }

  FMAT = corred$FMAT
  p = corred$p
  a = usin$prop$a
  CN =  usin$prop$CN

  # Death rate for each node is set my densitydependence. 1 means density dependent, 0 means density independent.
  d = densitydependence*usin$prop$d/usin$prop$B + (1-densitydependence)*usin$prop$d

  B = usin$prop$B
  isDetritus = usin$prop$isDetritus
  DetritusRecycling = usin$prop$DetritusRecycling
  isPlant = usin$prop$isPlant
  canIMM = usin$prop$canIMM

  # Check to make sure userdefinedinputs is either 1 or the same length as the nodes
  stopifnot(length(userdefinedinputs) %in% c(1, length(B)))

  # Rescale DetritusRecycling to sum to 1
  if(any(isDetritus > 0)){
    DetritusRecycling = DetritusRecycling/sum(DetritusRecycling)
  }

  # Calculate IN parameters for equilibrium
  IN = rep(0, length(B))
  delta = a*p*rowSums(FMAT) - colSums(FMAT) - (densitydependence*d*B*B + (1-densitydependence)*d*B)
  delta[isDetritus>0] =
    delta[isDetritus>0] +
    DetritusRecycling[isDetritus>0]*
    sum(c((1-a)*rowSums(FMAT),
          densitydependence*d*B*B + (1-densitydependence)*d*B))
  IN[rowSums(FMAT)==0] = -1*delta[rowSums(FMAT)==0]

  # Modify inputs for plants to make them per biomass if there is no inorganic nitrogen:
  cpn = r_i = rep(0, length(B)) # Create a new parameter class r_P for growth of a pool based on it's own biomass gain = r_i*X_i. Also create a vector cpn for inorganic nitrogen uptake rates if necessary. Used for plants.
  if(any(isPlant>0 & isFALSE(has_inorganic_nitrogen))){
    r_i[isPlant>0] = IN[isPlant>0]/B[isPlant >0]
    IN[isPlant>0] = 0
  }

  # If inorganic nitrogen is included, then calculate the amount available and the parameters at equilibrium:
  if(has_inorganic_nitrogen){

    if(sum(is.na(inorganic_nitrogen_properties)) > 1) stop("You must define two of the three parameters q, INN, and eqmN to solve the system.")

    if(any(isPlant>0)){
      IN_plant = sum(IN[isPlant >0]/CN[isPlant>0]) # Plant N uptake
    }else{
      IN_plant = 0 # Zero if there are no plants
    }
    # Get the inorganic nitrogen parameters (two of these must be values):
    Nmin = sum(comana(usin, shuffleTL = F)$Nmin)
    INN = inorganic_nitrogen_properties$INN
    if(!is.na(INN)){
      if(INN < 0) warning("INN < 0 is a constant loss rate of inorganic nitrogen.")
    }
    q = inorganic_nitrogen_properties$q
    if(!is.na(q)){if(q < 0) stop("q should be positive")}
    eqmN = inorganic_nitrogen_properties$eqmN
    if(!is.na(eqmN)){if(eqmN < 0) stop("eqmN (equilibrium nitrogen) should be positive")}

    if(is.na(INN)){
      INN = IN_plant + q*eqmN - Nmin
    }

    if(is.na(q)){
      q = (INN - IN_plant + Nmin)/eqmN
    }

    if(is.na(eqmN)){
      eqmN = (INN - IN_plant + Nmin)/q
    }

    if(verbose){
      print(paste0("Nitrogen parameters are: INN= ", INN, ", q= ", q, ", eqmN= ", eqmN))
    }

    # Get the plant coefficients cpn (cpn is indexed in units of plant carbon):
    cpn[isPlant > 0] = (IN[isPlant >0]/CN[isPlant>0])/(B[isPlant >0]*eqmN)
    IN[isPlant >0] = 0

    plant_N_growth = cpn*B*eqmN
  }else{
    plant_N_growth = rep(0,length(B))
    INN = inorganic_nitrogen_properties$INN
    q = inorganic_nitrogen_properties$q
  }

  # Add back in the system inputs
  delta = delta + IN + r_i*B + plant_N_growth*CN

  if(max(abs(delta)) > 1e-10) warning("Error in initial equilibrium. Check the code")

  #Replace any user defined input values

  if(length(userdefinedinputs) == 1 & any(is.na(userdefinedinputs))){
    userdefinedinputs = rep(userdefinedinputs, length(B))
  } # Make the vector the right length

  IN[!is.na(userdefinedinputs)] = userdefinedinputs[!is.na(userdefinedinputs)] # Replace user defined inputs whenever the value of the vector is not NA

  if(Conly){ # Use the modified p value for a carbon only model.
    pin = p
    warning("Running a carbon only model with modified production efficiency and diet determined by equilibrium conditions")
  }else{ # Use the original or maximum production efficiency for C and N model
    pin = usin$prop$p
  }

  PARS = c(Nnodes,d,a,pin,B,CN,isDetritus, isPlant,DetritusRecycling, canIMM, densitydependence,cij,DIETLIMTS, IN, r_i, sum(diet_correct), sum(Conly), has_inorganic_nitrogen, cpn,INN,q, NA)

  names(PARS) = c("Nnodes",
                  paste0("d_", 1:Nnodes),
                  paste0("a_", 1:Nnodes),
                  paste0("p_", 1:Nnodes),
                  paste0("B_", 1:Nnodes),
                  paste0("CN_", 1:Nnodes),
                  paste0("isDetritus_", 1:Nnodes),
                  paste0("isPlant_", 1:Nnodes),
                  paste0("DetritusRecycling_", 1:Nnodes),
                  paste0("canIMM_", 1:Nnodes),
                  paste0("densitydependence_", 1:Nnodes),
                  paste0("c_", rep(1:Nnodes, each = Nnodes), rep(1:Nnodes, Nnodes)),
                  paste0("DIETLIMITS_", rep(1:Nnodes, each = Nnodes), rep(1:Nnodes, Nnodes)),
                  paste0("IN_", 1:Nnodes),
                  paste0("r_", 1:Nnodes),
                  "diet_correct",
                  "Conly",
                  "hasIORGN",
                  paste0("cpn_", 1:Nnodes),
                  "INN",
                  "q",
                  "DetExpt")

  # Add an flag to determine whether the output is for a stability analysis:
  PARS["forstability"] = 0

  # Add an flag to determine whether the output should produce all nitrogen values even if C:N ratio stays the same:
  PARS["keepallnitrogen"] = 1

  if(has_inorganic_nitrogen){
    yint = c(usin$prop$B,usin$prop$B/usin$prop$CN, eqmN)
    names(yint) = c(paste0(colnames(usin$imat),"C"),
                    paste0(colnames(usin$imat),"N"), "IORGN")
  }else{
    yint = c(usin$prop$B,usin$prop$B/usin$prop$CN)
    names(yint) = c(paste0(colnames(usin$imat),"C"),
                    paste0(colnames(usin$imat),"N"))
  }

  if(returnCNmin){
    Cmin = a*(1-p)*rowSums(FMAT)
    Cfaeces = (1-a)*rowSums(FMAT)
    Nmin = rowSums(
      (matrix(1/CN, nrow = Nnodes, ncol = Nnodes, byrow = T) -
         matrix(p/CN, nrow=Nnodes, ncol = Nnodes))*
        matrix(a, nrow=Nnodes, ncol = Nnodes)*FMAT)
    Nfaeces = (1-a)*rowSums(FMAT*matrix(1/CN, nrow = Nnodes, ncol = Nnodes, byrow = T))
    names(Cmin) = names(Nmin) = names(Cfaeces) = names(Nfaeces) = usin$prop$ID
    return(list(PARS = PARS, yint = yint, Cmin = Cmin, Cfaeces = Cfaeces, Nmin = Nmin, Nfaeces = Nfaeces))
  }else{
    return(list(PARS = PARS, yint = yint))
  }
}
