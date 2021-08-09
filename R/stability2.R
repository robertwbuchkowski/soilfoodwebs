#' A function to run the stability analysis using the numerical simulation of the Jacobain matrix.
#'
#' @param usin The community for which to calculate stability.
#' @param DIETLIMTS The diet limits matrix for the stoichiometry correction (proportion of diet)?
#' @param diet_correct Boolean: Does the organism correct it's diet?
#' @param Conly Boolean: Is the model meant for carbon only?
#' @param userdefinedinputs Do you want to input a user defined vector of input functions? If NA the input values that keep the system at equilibrium are calculated. If not, put in a vector of input rates for each node.
#' @param has_inorganic_nitrogen Boolean: Is there an inorganic nitrogen pool?
#' @param inorganic_nitrogen_properties A list of state variables for the inorganic nitrogen pool (INN = inputs, q = per capita loss of N, eqmN = equilibrium N). Must include a value for two of the three variables and has the final one as NA.
#' @param forstabilityonly Boolean: Do you only want to check the stability of the carbon equations and inorganic nitrogen? If FALSE then changes in detritus nitrogen are also included.
#' @param densitydependence Which nodes have density dependence? NA default means none of them do. Should be a vector of 0 (no DD) and 1 (DD) for each node.
#' @return The stability calculated by the Jacobian as a list of the Jacobian, eigenvalues, and the maximum eigenvalue.
#' @details
#' The stability as defined by the Jacobian is estimated using the ordinary differential equations. Consequenlty, the parameters are retrieved using \code{\link{getPARAMS}} and then added to the function \code{\link[rootSolve]{jacobian.full}} using the ODE defined by \code{\link{foodwebode}} This stability function allows the user to define density-dependence by node as present or absent. If you want to apply a uniform level of density-dependence use the function 'stability' and choose the option Moorecobain.
#' @seealso
#' \code{\link{getPARAMS}} and \code{\link{foodwebode}} for functions which are called internally and \code{\link[rootSolve]{jacobian.full}} for the method used.
#'
#'@examples
#'# Basic stability calculation:
#'stability2(intro_comm)
#'
#'# Stability calculation with density-dependence for the animals:
#'stability2(intro_comm, densitydependence = c(1, 1, 1, 0, 0, 0, 0))
#' @export
stability2 <- function(usin,
                       DIETLIMTS = NA,
                       diet_correct = T,
                       Conly = F,
                       userdefinedinputs = NA,
                       has_inorganic_nitrogen = F,
                       inorganic_nitrogen_properties = list(INN = NA, q = NA, eqmN = NA),
                       forstabilityonly = T,
                       densitydependence = NA){

  if(!Conly) usin = corrstoich(usin, dietlimits = DIETLIMTS)

  PARSandY = getPARAMS(usin = usin,
                       DIETLIMTS = DIETLIMTS,
                       diet_correct = diet_correct,
                       Conly = Conly,
                       densitydependence = densitydependence,
                       userdefinedinputs = userdefinedinputs,
                       has_inorganic_nitrogen =has_inorganic_nitrogen,
                       inorganic_nitrogen_properties = inorganic_nitrogen_properties
  )

  parameters = PARSandY$PARS

  parameters["forstability"] = ifelse(forstabilityonly, 1,0)
  parameters["keepallnitrogen"] = 0

  Nnodes = parameters["Nnodes"]

  if(forstabilityonly){
    if(has_inorganic_nitrogen){
      YINT = PARSandY$yint[c(1:Nnodes,length(PARSandY$yint))]
    }else{
      YINT = PARSandY$yint[1:Nnodes]
    }
  }else{
    detinyint = c((Nnodes+1):(2*Nnodes))[usin$prop$isDetritus == 1]
    if(has_inorganic_nitrogen){
      YINT = PARSandY$yint[c(1:Nnodes,detinyint,length(PARSandY$yint))]
    }else{
      YINT = PARSandY$yint[c(1:Nnodes,detinyint)]
    }
  }

  J = rootSolve::jacobian.full(y=YINT,
                               func = foodwebode,
                               parms = parameters)

  return(list(J = J, # Return the Jacobian
              eigen = base::eigen(J), # Return the eigenvalues and vectors
              rmax = max(Re(base::eigen(J)$values)))) # Return the maximum eigenvalue. The system is stable if this is negative

}
