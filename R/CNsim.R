#' A function to simulate the dynamics over time wrapping getPARAMS and foodwebode
#'
#' @param usin The community you want to simulate.
#' @param DIETLIMTS The diet limits matrix for the stoichiometry correction (proportion of diet)?
#' @param diet_correct Boolean: Does the organism correct it's diet?
#' @param Conly Boolean: Is the model meant for carbon only?
#' @param userdefinedinputs Do you want to input a user defined vector of input functions? If NA the input values that keep the system at equilibrium are calculated. If not, put in a vector of input rates for each node.
#' @param start_mod A vector of modifications to the starting conditions, which default to the biomass vector in the community. This vector is multiplied by the biomass vector.
#' @param TIMES The vector of times that you want to run the model. Defaults to 1 to 100 by 1.
#' @param keepallnitrogen Boolean: Keep all the nitrogen pools in the model output? Will be set to FALSE if you are running a stability analysis using the function stability2().
#' @param has_inorganic_nitrogen Boolean: Is there an inorganic nitrogen pool?
#' @param densitydependence Which nodes have density dependence? NA default means none of them do. Should be a vector of 0 (no DD) and 1 (DD) for each node.
#' @param inorganic_nitrogen_properties A list of state variables for the inorganic nitrogen pool (INN = inputs, q = per capita loss of N, eqmN = equilibrium N). Must include a value for two of the three variables and has the final one as NA.
#' @param DETEXPT The pool that should be used for the detritus experiment by name or position in the vector usin$prop$ID.
#' @param DETEXPTSTART The start of the detritus experiment. Defaults to 100.
#' @return The output of the simulation.
#' @details
#' A function that simulates the food web over the user defined times. If you do not modify the starting state using start_mod or add a detritus experiment using DETEXPT, then the result will just be a flat line if the food web is stable.
#' @seealso
#' The function uses \code{\link{getPARAMS}} and \code{\link{foodwebode}} to simulate the community.
#' @examples
#' # Basic example with a 5% reduction in predator biomass:
#' CNsim(intro_comm, start_mod = c(0.95, 1,1,1,1,1,1))
#'
#' # Simulate a decomposition experiment:
#' CNsim(intro_comm, DETEXPT = c("Detritus1"), DETEXPTSTART = 100)
#'
#' @export
CNsim <- function(usin,
                  DIETLIMTS = NA,
                  diet_correct = T,
                  Conly = F,
                  userdefinedinputs = NA,
                  start_mod = NA,
                  TIMES = 1:100,
                  keepallnitrogen = T,
                  has_inorganic_nitrogen = F,
                  densitydependence = NA,
                  inorganic_nitrogen_properties = list(INN = NA, q = NA, eqmN = NA),
                  DETEXPT = NA,
                  DETEXPTSTART = NA){
  # Correct stoichiometry if necessary:
  if(Conly){
    usin = usin
  }else{
    usin = corrstoich(usin, forceProd = !diet_correct, dietlimits = DIETLIMTS)
  }

  parEin = getPARAMS(usin = usin,
                     DIETLIMTS = DIETLIMTS,
                     diet_correct = diet_correct,
                     Conly = Conly,
                     userdefinedinputs = userdefinedinputs,
                     has_inorganic_nitrogen = has_inorganic_nitrogen,
                     densitydependence = densitydependence,
                     inorganic_nitrogen_properties = inorganic_nitrogen_properties)

  if(is.character(DETEXPT)){
    DETEXPT = which(usin$prop$ID == DETEXPT)
  }

  if(!is.na(DETEXPT) & is.na(DETEXPTSTART)){
    DETEXPTSTART = 100
    warning("Making the starting detritus experiment 100 units of C.")
  }

  PARS = parEin$PARS
  y0 = parEin$yint

  if(!is.na(DETEXPT) &!is.na(DETEXPTSTART)){
    PARS["DetExpt"] = DETEXPT
    y0["DetExptState"] = DETEXPTSTART
  }

  # Fix the starting modifier vector
  if(length(start_mod) == 1){
    if(is.na(start_mod)){
      start_mod = rep(1, PARS["Nnodes"])
    }
  }

  # Double the start_mod vector for C and N
  start_mod = c(start_mod, start_mod)

  # Correct the length of the start_mod to add DETEXPT or IORGN:
  nta = length(y0) - length(start_mod)
  start_mod = c(start_mod, rep(1, nta))

  # Modify the starting vector:
  yint = y0*start_mod

  # Subset the starting vector so it only includes nitrogen for the detirtus pools wherein C:N ratio is flexible
  logdet = which(usin$prop$isDetritus == 1) # logical for the detritus pools
  Nnodes = PARS["Nnodes"] # number of nodes

  # Add in the extra pools
  extrapools = ifelse(has_inorganic_nitrogen, 1,0)
  extrapools = ifelse(!is.na(DETEXPT), extrapools+1,extrapools)
  posextra = (length(y0) - extrapools + 1):length(y0)

  # Get rid of the nitrogen parameters if keeping all nitrogen is true
  if(!keepallnitrogen){
    PARS["keepallnitrogen"] = 0
    if(extrapools > 0){
      yint = yint[c(1:Nnodes, # keep all the carbon state variables
                    c((Nnodes + 1):(2*Nnodes))[logdet], # keep only the detritus pools as nitrogen state variables
                    posextra # Add in any further state variables (detritus experiment or inorganic N)
      )]
    }else{
      yint = yint[c(1:Nnodes, # keep all the carbon state variables
                    c((Nnodes + 1):(2*Nnodes))[logdet] # keep only the detritus pools as nitrogen state variables
      )]
    }
  }

  # Simulate the model
  output = deSolve::ode(y = yint, times = TIMES, func = foodwebode, parms = PARS)

  output = data.frame(output)

  # Prepare for naming:

  Nnodes = PARS["Nnodes"]

  SVnames = sapply(names(y0)[1:Nnodes], function(XX) substr(XX,1,nchar(XX)-1))

  if(keepallnitrogen){
    SVnames2 = SVnames
  }else{
    SVnames2 = sapply(names(y0)[logdet], function(XX) substr(XX,1,nchar(XX)-1))
  }

  NdetritusPools = sum(PARS[grepl("isDetritus",names(PARS))])

  if(!is.na(DETEXPT)){
    DETINFO = c("DetExpt",paste0("k",seq(1:NdetritusPools),"_decompconst"))
  }else{
    DETINFO = paste0("k",seq(1:NdetritusPools),"_decompconst")
  }

  # Add inorgnaic nitrogen pool name if applicable:
  if(has_inorganic_nitrogen){
    DETINFO = c("IORGN", DETINFO)
  }

  colnames(output) = c("Day",
                       paste0(SVnames, "_Carbon"),
                       paste0(SVnames2, "_Nitrogen"),
                       DETINFO,
                       paste0(SVnames,"_Cmin"),
                       paste0(SVnames,"_Nmin"),
                       paste0(SVnames,"_prodeff"),
                       paste0(rep(SVnames, times = Nnodes), "-eats-",rep(SVnames, each = Nnodes), "_consrate"))

  # Get rid of all output columns that are just all zero (usually consumption rates between pools that don't eat each other.)
  output = output[,!colSums(abs(output)) == 0]

  return(output)
}
