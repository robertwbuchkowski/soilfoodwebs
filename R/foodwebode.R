#' A function to simulation the food webs away from equilibrium.
#'
#' @param t The ODE time.
#' @param y The ODE simulation start.
#' @param pars The ODE parameters.
#' @details
#' The food web model ode for simulating the model over time: requires y inputs and parameters from getPARAM function: it adds a detritus decomposition experiment. The only pools with nitrogen values are the detritus, because all others have fixed C:N ratios.
#' @return The changes in each node biomass along with parameters and mineralization rates.
foodwebode <- function(t,y,pars){
  Nnodes = pars[1]
  d = pars[grepl("d_",names(pars))]
  a = pars[grepl("a_",names(pars))]
  p = pars[grepl("p_",names(pars))]
  IN_CN = pars[grepl("CN_",names(pars))]
  isDetritus = pars[grepl("isDetritus_",names(pars))]
  isPlant = pars[grepl("isPlant_",names(pars))]
  DetritusRecycling = pars[grepl("DetritusRecycling_",names(pars))]
  canIMM = pars[grepl("canIMM_",names(pars))]
  densitydependence = pars[grepl("densitydependence_",names(pars))]
  cij = matrix(pars[grepl("c_",names(pars))], nrow = Nnodes, ncol = Nnodes)
  DIETLIMITS = matrix(pars[grepl("DIETLIMITS_",names(pars))], nrow = Nnodes, ncol = Nnodes)
  IN = pars[grepl("IN_",names(pars))]
  r_i = pars[grepl("r_",names(pars))]
  diet_correct = pars[grepl("diet_correct",names(pars))] == 1
  Conly = pars[grepl("Conly",names(pars))] == 1

  # Load state variables:
  B = y[1:Nnodes]

  if(pars["keepallnitrogen"]){
    Nbiomass = y[(Nnodes+1):(2*Nnodes)]
  }else{
    Nbiomass = B/IN_CN # Nitrogen is fixed for most state variables
    # Fixed nitrogen for detritus pools:
    Numdetritus = sum(isDetritus)
    if(!pars["forstability"]){
      Nbiomass[isDetritus == 1] = y[(Nnodes+1):(Nnodes + Numdetritus)]
    }
  }

  # Load information about whether there is inorgnaic N or a detritus experiment:
  DetExpt = pars[grepl("DetExpt",names(pars))]
  hasIORGN = pars[grepl("hasIORGN", names(pars))]

  # Detritus experiment:
  if(!is.na(DetExpt)){
    DetExptState = y[length(y)] # Always the last state variable
  }

  # Inorganic nitrogen:
  cpn = pars[grepl("cpn_", names(pars))]
  INN = pars[grepl("INN", names(pars))]
  q = pars[grepl("q", names(pars))]
  if(hasIORGN){ # The last state variable without DetExpt, otherwise the second last
    if(!is.na(DetExpt)){
      IORGN = y[length(y)-1]
    }else{
      IORGN = y[length(y)]
    }
  }

  # Calculate new C:N ratios:
  CN = round(B/Nbiomass, 5) # This rounding occurs because of numeric instability in the calculation of C:N ratio. Small changes can lead to large differences and cause errors in our assessment of change. C:N ratio differences after 5 decimal places probably don't matter

  # Calculate the inorganic nitrogen that will be available for immobilization based on the current amount less losses and plant uptake, currently plants have precident for N uptake!!
  if(hasIORGN){
    N_avail = INN + IORGN - q*IORGN - sum(cpn[isPlant == 1]*B[isPlant == 1]*IORGN)
  }else{
    N_avail = Inf
  }

  # produce the correct vector of plant N and C uptake rates:
  #Set up rates regardless
  plant_growthC = plant_growthN = rep(0, Nnodes)

  if(hasIORGN){
    plant_growthN[isPlant == 1] = cpn[isPlant == 1]*B[isPlant == 1]*IORGN
    plant_growthC[isPlant == 1] = plant_growthN[isPlant == 1]*CN[isPlant == 1]
    if(!all(r_i == 0)) warning("Plants shouldn't grow via r_i and inorganic nitrogen uptake.")
  }
  # Correct the diet and production efficiency
  Corred = correction_function(B, cij, CN, p, a, canIMM, dietlimits = DIETLIMITS, diet_correct = diet_correct, Conly = Conly, Immobilizationlimit = N_avail)
  FMAT = Corred$FMAT
  p = Corred$p

  delta = a*p*rowSums(FMAT) - colSums(FMAT) - densitydependence*d*B^2 - (1-densitydependence)*d*B + IN + r_i*B + plant_growthC
  delta[isDetritus>0] = delta[isDetritus>0] + DetritusRecycling[isDetritus>0]*sum(c((1-a)*rowSums(FMAT), (densitydependence*d*B^2 + (1-densitydependence)*d*B))) # add back in dead biomass and rejected food to detritus

  # Nitrogen equations (C:N constant for all except detritus)
  deltaN = a*FMAT%*%(1/CN) - colSums(FMAT)/CN - densitydependence*d*(B^2)/CN - (1-densitydependence)*d*B/CN + (IN + r_i*B)/IN_CN + plant_growthN
  deltaN[isDetritus>0] = deltaN[isDetritus>0] + DetritusRecycling[isDetritus>0]*sum(c((1-a)*FMAT%*%(1/CN), (densitydependence*d*(B^2)/CN + (1-densitydependence)*d*B/CN))) # add back in dead biomass and rejected food to detritus

  # calculate k, which is the sum of all consumed detritus divided by the detritus pool size
  if(sum(isDetritus > 0) > 1){
    k = colSums(FMAT[,isDetritus > 0])/B[isDetritus > 0]
  }else{
    k = sum(FMAT[,isDetritus > 0])/B[isDetritus > 0]
  }

  # calculate consumption out of "Test detritus" pool
  if(!is.na(DetExpt)) dDetExpt = - sum(FMAT[,DetExpt])*DetExptState/B[DetExpt]

  # Calculate carbon and nitrogen mineralization
  Cmin = a*(1-p)*rowSums(FMAT)
  Nmin = rowSums(
    (matrix(1/CN, nrow = Nnodes, ncol = Nnodes, byrow = T) -
       matrix(p/CN, nrow=Nnodes, ncol = Nnodes))*
      matrix(a, nrow=Nnodes, ncol = Nnodes)*FMAT)

  # Take the mineralized N out of the pools
  deltaN = deltaN - Nmin

  # Confirm that fixed C:N ratio is staying fixed

  changenotdet = abs(delta) > 1e-04 & isDetritus == 0 # Identify the non-detritus pools that are actually changing in pool size.

  if(sum(changenotdet) >0 &
     !isTRUE(all.equal(unname(c(delta/deltaN)[changenotdet]), unname(IN_CN)[changenotdet], tolerance = 1e-4)) # Check if the ratio of change for non-Detritus pools that are changing size is still the C:N ratio.
  ){
    warning("Constant C:N ratios appear to be changing. Check model formulation to make sure this isn't an issue by simulating the model with and without keepallnitrogen = TRUE.")
    toprint = abs(c(delta/deltaN)[changenotdet] - IN_CN[changenotdet]) > 1e-4
    arrind = which(toprint, arr.ind = T)
    for(i in 1:sum(toprint)){
      print(paste0("At time step ",t, " the actual C:N ratio of change was  ", c(delta/deltaN)[changenotdet][arrind[i]], " but should have been ",IN_CN[changenotdet][arrind[i]]))
    }
  }

  # Add N to inorganic pool if necessary:
  if(hasIORGN) dIORGN = INN -q*IORGN - sum(cpn[isPlant == 1]*B[isPlant == 1]*IORGN) + sum(Nmin)

  # Check to make sure the user isn't trying to run stability AND keep all the nitrogen pools:
  if(pars["forstability"] & pars["keepallnitrogen"]) stop("You cannot check stability and keep all the nitrogen state variables. This will produce a Jacobian that is not full rank.")

  # Calculate what to return:
  if(pars["forstability"]){
    if(!is.na(DetExpt)) stop("Cannot run a detritus experiment and test stability at the same time.")

    if(hasIORGN){
      deltaF = c(delta,dIORGN)
    }else{
      deltaF = c(delta)
    }

    OUTPUT = c(list(deltaF))
  }else{
    if(pars["keepallnitrogen"]){
      if(hasIORGN){
        if(!is.na(DetExpt)){
          deltaF = c(delta, deltaN,dIORGN, dDetExpt)
        }else{
          deltaF = c(delta, deltaN,dIORGN)
        }
      }else{
        if(!is.na(DetExpt)){
          deltaF = c(delta, deltaN,dDetExpt)
        }else{
          deltaF = c(delta, deltaN)
        }
      }
    }else{
      if(hasIORGN){
        if(!is.na(DetExpt)){
          deltaF = c(delta, deltaN[isDetritus == 1],dIORGN, dDetExpt)
        }else{
          deltaF = c(delta, deltaN[isDetritus == 1],dIORGN)
        }
      }else{
        if(!is.na(DetExpt)){
          deltaF = c(delta, deltaN[isDetritus == 1],dDetExpt)
        }else{
          deltaF = c(delta, deltaN[isDetritus == 1])
        }
      }
    }

    OUTPUT = c(list(unname(deltaF)),k = k,Cmin = Cmin,Nmin = Nmin, p = p, cij = cij)
  }
  return(OUTPUT)
}
