#' Calculate the strength of stability as smin.
#'
#' @param usin The community for which to calculate smin
#' @return The values of smin.
#' @details
#' Calculates the value of smin in each web using the metric used by Moore and de Ruiter.
#' @examples
#' # Calculate the minimum values of s for the introduction community.
#' calc_smin(intro_comm)
#' @export
calc_smin <- function(usin){
  # Check the community for any errors before calculations
  usin = checkcomm(usin)

  # Start the calculations:
  Smin = 1 # Set the starting value for Smin at 1
  stabval = calculate_smin(SMIN = Smin, usin, isSQR = F)

  # Increase Smin until the community is stable
  while(stabval >= 0){
    Smin = Smin*10
    stabval = calculate_smin(SMIN = Smin, usin, isSQR = F)

  }

  # Decrease Smin until the community becomes unstable
  while(stabval < 0){
    Smin = Smin/10
    stabval = calculate_smin(SMIN = Smin, usin, isSQR = F)
  }

  # Increase Smin again until it is stable, this is the final step
  while(stabval >=0){
    Smin = Smin+Smin
    stabval = calculate_smin(SMIN = Smin, usin, isSQR = F)
  }
  return(Smin)
}

#' Function used inside the calc_smin function
#'
#' @param SMIN The value of smin.
#' @param usin The community used in the calculation
#' @param isSQR Boolean: Should the rmax value be squared?
#' @return The value of rmax for the community and given SMIN
calculate_smin <- function(SMIN = 1, usin, isSQR = T){
  if(isSQR){
    (stability(usin, method = "Moorecobian", smin = SMIN)$rmax*1e4)^2
  }else{
    stability(usin, method = "Moorecobian", smin = SMIN)$rmax
  }
}
