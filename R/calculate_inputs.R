#' A function to calculate the inputs and outputs at equilibrium and print them for the user.
#'
#' @param usin The community you want to simulate.
#' @param DIETLIMTS The diet limits matrix for the stoichiometry correction (proportion of diet)?
#' @param diet_correct Boolean: Does the organism correct it's diet?
#' @param Conly Boolean: Is the model meant for carbon only?
#' @param verbose Boolean: Do you want to print the summary in the Console?
#' @param toround Boolean: Should the answer be rounded?
#' @return A list of inputs and mineralization rates.
#' @examples
#' # Calculate the inputs and outputs of a community. Prints a summary by default and saves a list.
#' calculate_inputs(intro_comm)
#' @export
calculate_inputs <- function(usin,
                             DIETLIMTS = NA,
                             diet_correct = T,
                             Conly = F,
                             verbose = T,
                             toround = T
){

  parEI = getPARAMS(usin = usin,
                    DIETLIMTS = DIETLIMTS,
                    diet_correct = diet_correct,
                    Conly = Conly,
                    returnCNmin = T)

  Cmin = parEI$Cmin
  Nmin = parEI$Nmin

  # Identify detritus using the detritus recycling variable
  DetritusRecycling = parEI$PARS[grepl("DetritusRecycling_",names(parEI$PARS))]

  # Get IN
  IN = parEI$PARS[grepl("IN_",names(parEI$PARS))]
  names(IN) = names(Cmin)

  # For detritus, set IN to zero because I calculate the flux below because of biased flows:
  INsave = IN
  IN[DetritusRecycling > 0] = 0

  # Get CN
  CN = usin$prop$CN

  # Get growth rates
  r_i = parEI$PARS[grepl("r_",names(parEI$PARS))]
  names(r_i) = names(Cmin)

  # Carbon rates
  IN2 = IN[IN > 0]
  IN3 = IN[IN < 0]
  r_i2 = r_i[r_i > 0]*usin$prop$B[r_i >0]

  # Nitrogen rates
  IN2N = IN[IN > 0]/CN[IN > 0]
  IN3N = IN[IN < 0]/CN[IN < 0]
  r_i2N = r_i[r_i > 0]*usin$prop$B[r_i >0]/CN[r_i > 0]

  # Calculate the death rate to report when there is not detritus:

  deathC = sum(parEI$PARS[grepl("d_",names(parEI$PARS))]*usin$prop$B)
  deathN = sum(parEI$PARS[grepl("d_",names(parEI$PARS))]*usin$prop$B/usin$prop$CN)

  # C and N for lost detritus (calculate the net flow):
  if(sum(DetritusRecycling>0)> 0){
    if(sum(DetritusRecycling>0) == 1){
      flux_det_C = sum((deathC + sum(parEI$Cfaeces) - sum(comana(corrstoich(usin))$fmat[,DetritusRecycling >0]))*-1)

      flux_det_N = sum((deathN + sum(parEI$Nfaeces) - sum(comana(corrstoich(usin))$fmat[,DetritusRecycling >0])/usin$prop$CN[DetritusRecycling >0])*-1)

      if(isFALSE(all.equal(flux_det_C, unname(sum(INsave[DetritusRecycling >0]))))){
        warning("Different rates for detritus calculations. Check the calculations")
      }
    }else{
      flux_det_C = (deathC + sum(parEI$Cfaeces) - sum(colSums(comana(corrstoich(usin))$fmat[,DetritusRecycling >0])))*-1

      flux_det_N = (deathN + sum(parEI$Nfaeces) - sum(colSums(comana(corrstoich(usin))$fmat[,DetritusRecycling >0])/usin$prop$CN[DetritusRecycling >0]))*-1

      if(isFALSE(all.equal(flux_det_C, unname(INsave[DetritusRecycling >0])))){
        warning("Different rates for detritus calculations. Check the calculations")
      }
    }
  }

  if(toround){
    if(verbose){
      print("------------------------------------------------")
      print("Equilibrium is maintained by these inputs:")
      if(any(IN >0)){
        for(i in 1:length(IN2)){
          print(paste0(round(IN2[i])," units of C to ", names(IN2)[i]))
          print(paste0(round(IN2N[i])," units of N to ", names(IN2)[i]))
        }
      }
      if(any(r_i >0)){
        for(i in 1:length(r_i2)){
          print(paste0(round(r_i2[i])," units of C to ", names(r_i2)[i]))
          print(paste0(round(r_i2N[i])," units of N to ", names(r_i2)[i]))
        }
      }
      print("------------------------------------------------")
      if(any(IN < 0)){
        print("The system looses these organic compounds:")
        for(i in 1:length(IN3)){
          print(paste0(round(IN3[i])," units of C from ", names(IN3)[i]))
          print(paste0(round(IN3N[i])," units of N from ", names(IN3)[i]))
        }
        print("------------------------------------------------")
      }
      if(all(DetritusRecycling == 0)){
        print("The system looses these organic compounds because there is no recycling to detritus:")

        print(paste0("Looses of C from death rate:", round(deathC)))
        print(paste0("Looses of N from death rate:", round(deathN)))
        print(paste0("Looses of C from faeces:", round(sum(parEI$Cfaeces))))
        print(paste0("Looses of N from faeces:", round(sum(parEI$Nfaeces))))
        print("------------------------------------------------")
      }else{
        if(flux_det_C < 0){
          print("The system looses organic compounds from the detritus pool:")

          print(paste0("Looses of C from detritus:", round(-1*flux_det_C)))
          print(paste0("Looses of N from detritus:", round(-1*flux_det_N)))
          print(paste0("Looses detritus at a C:N ratio of:", round(flux_det_C/flux_det_N)))
        }else{
          print("The system gains organic compounds from the detritus pool:")

          print(paste0("Gains of C from detritus:", round(flux_det_C)))
          print(paste0("Gains of N from detritus:", round(flux_det_N)))
          print(paste0("Gains detritus at a C:N ratio of:", round(flux_det_C/flux_det_N)))
        }
        print("------------------------------------------------")
      }
      print(paste0("The web mineralizes ", round(sum(Cmin)), " units of carbon"))
      if(sum(Nmin)>0){
        print(paste0("The web mineralizes ", round(sum(Nmin)), " units of nitrogen"))
      }else{
        print(paste0("The web must immobilize ", round(abs(sum(Nmin))), " units of nitrogen"))
      }
    }
  }else{
    if(verbose){
      print("------------------------------------------------")
      print("Equilibrium is maintained by these inputs:")
      if(any(IN >0)){
        for(i in 1:length(IN2)){
          print(paste0(IN2[i]," units of C to ", names(IN2)[i]))
          print(paste0(IN2N[i]," units of N to ", names(IN2)[i]))
        }
      }
      if(any(r_i >0)){
        for(i in 1:length(r_i2)){
          print(paste0(r_i2[i]," units of C to ", names(r_i2)[i]))
          print(paste0(r_i2N[i]," units of N to ", names(r_i2)[i]))
        }
      }
      print("------------------------------------------------")
      if(any(IN < 0)){
        print("The system looses these organic compounds:")
        for(i in 1:length(IN3)){
          print(paste0(IN3[i]," units of C from ", names(IN3)[i]))
          print(paste0(IN3N[i]," units of N from ", names(IN3)[i]))
        }
        print("------------------------------------------------")
      }
      if(all(DetritusRecycling == 0)){
        print("The system looses these organic compounds because there is no recycling to detritus:")

        print(paste0("Losses of C from death rate:", deathC))
        print(paste0("Losses of N from death rate:", deathN))
        print(paste0("Losses of C from faeces:", sum(parEI$Cfaeces)))
        print(paste0("Losses of N from faeces:", sum(parEI$Nfaeces)))
        print("------------------------------------------------")
      }else{
        if(flux_det_C < 0){
          print("The system looses organic compounds from the detritus pool:")

          print(paste0("Looses of C from detritus:", -1*flux_det_C))
          print(paste0("Looses of N from detritus:", -1*flux_det_N))
          print(paste0("Looses detritus at a C:N ratio of:", flux_det_C/flux_det_N))
        }else{
          print("The system gains organic compounds from the detritus pool:")

          print(paste0("Gains of C from detritus:", flux_det_C))
          print(paste0("Gains of N from detritus:", flux_det_N))
          print(paste0("Gains detritus at a C:N ratio of:", flux_det_C/flux_det_N))
        }
        print("------------------------------------------------")
      }
      print(paste0("The web mineralizes ", sum(Cmin), " units of carbon"))
      if(sum(Nmin)>0){
        print(paste0("The web mineralizes ", sum(Nmin), " units of nitrogen"))
      }else{
        print(paste0("The web must immobilize ", abs(sum(Nmin)), " units of nitrogen"))
      }
    }
  }


  return(list(IN = IN, r_i = r_i,Cmin = Cmin, Nmin = Nmin))

}
