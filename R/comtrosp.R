#' A function to combine trophic species.
#'
#' @param usin # The community were you are combining trophic species
#' @param selected Select trophic species to combine, which ones do you want to combine (vector of names)? If left as NA, the most similar trophic species are combined.
#' @param  deleteCOMBOcannibal Boolean: Do you want to delete the cannibalism that may have been created by combining two trophic species (T) or leave it in the model (F)?
#' @param allFEEDING1 Boolean: Do you want to return all feeding preferences to 1 (T), or would you like to set the feeding preferences of the newly combined trophic species as the biomass-weighted average of the old ones (F)?
#' @return The new community with the seltf or most similar trophic species combined.
#' @details
#' The function combines trophic species by merging them in both the community matrix (imat) and the properties database (prop). The user can select the two or more trophic species to be combined using the selected option. If this is left as the default (NA), then the two most similar trophic species are combined by comparing their feeding relationships.
#' @examples
#' # Combine two trophic species that are the most similar.
#' comtrosp(intro_comm)
#'
#' # Combine two selected trophic species.
#' comtrosp(intro_comm, selected = c("Orib1", "Predator"))
#'
#' # Combine two selected trophic species
#' # Remove the cannibalism that is created
#' # Rescale all feeding preferences to 1.
#' comtrosp(intro_comm,
#' selected = c("Orib1", "Predator"),
#' deleteCOMBOcannibal = TRUE,
#' allFEEDING1 = TRUE)
#' @export

comtrosp <- function(usin,
                     selected = NA,
                     deleteCOMBOcannibal = F,
                     allFEEDING1 = F
){
  seltf = !any(is.na(selected))
  if(seltf == T & length(selected) > 2){
    # If the length is longer than 2, we need to run the function comtrosp_backend multiple times.
    NTC = length(selected) - 1 # Number of combinations
    usin_temp = usin # Save temporary file
    for(NTCsel in 1:NTC){ # Loop over the combinations
      # Find the new combination
      newTS = usin_temp$prop$ID[grepl(selected[1],usin_temp$prop$ID)]
      # Build a temporary selected vector
      selected_temp = c(newTS,selected[NTCsel + 1])
      # Run combination
      usin_temp = comtrosp_backend(usin = usin_temp, seltf = seltf, selected = selected_temp, deleteCOMBOcannibal = deleteCOMBOcannibal, allFEEDING1 = allFEEDING1)
    }
    return(checkcomm(usin_temp))
  }else{
    # Only need to run the function once
    usin_temp = comtrosp_backend(usin = usin, seltf = seltf, selected = selected, deleteCOMBOcannibal = deleteCOMBOcannibal, allFEEDING1 = allFEEDING1)
    return(checkcomm(usin_temp))
  }
}

#' A function to combine trophic species with only two seltf options available
#'
#' @param usin, # The community were you are combining trophic species
#' @param seltf Boolean: Do you want to select the trophic species to combine or do you want to combine two based on their similarity?
#' @param selected A vector of two trophic species to combine. The default is NA, which is not appropriate if seltf is T.
#' @param  deleteCOMBOcannibal Boolean: Do you want to delete the cannibalism that may have been created by combining two trophic species (T) or leave it in the model (F)?
#' @param allFEEDING1 Boolean: Do you want to return all feeding preferences to 1 (T), or would you like to set the feeding preferences of the newly combined trophic species as the biomass-weighted average of the old ones (F)?
#' @return The new community with the seltf or most similar trophic species combined.

comtrosp_backend <- function(usin, seltf = F, selected = NA, deleteCOMBOcannibal = F, allFEEDING1 = F){

  imat = usin$imat
  prop = usin$prop

  imat_mod_update = imat_mod = imat # Create some replicates of the feeding matrix to store intermediate steps

  # If seltf trophic species are combined, then just set them in the idmax, which is the matrix identifying their position in the web
  if(seltf){
    idmax = matrix(which(rownames(imat) %in% selected), ncol = 2, nrow = 1)
  }else{

    # If you don't select the species to be combined, then the most similar species are calcualted by comparing the number of matching trophic links

    # Matrix to save how similar each pair of species are
    simidx = matrix(NA, ncol = dim(imat_mod)[1], nrow = dim(imat_mod)[1])
    colnames(simidx) = rownames(simidx) = rownames(imat_mod)

    # Loop through and calculate that similarity
    for(ii in 1:dim(imat_mod)[1]){
      for(jj in 1:dim(imat_mod)[2]){
        simidx[ii,jj] = unname(table(c(imat_mod[ii,] == imat_mod[jj,],imat_mod[,ii] == imat_mod[,jj]))["TRUE"])
      }
    }

    diag(simidx) = -1 # Set the diagonal equal to -1, because each species is most similar to itself and we want to identify the most different ones

    simidx[lower.tri(simidx)] = -1 # Set the lower triangle of the matrix to -1, because these values repeat the upper triangle and create duplicates that are not appropriate

    # Find the most similar tropho-species
    idmax = which(simidx == max(simidx), arr.ind=T)

    # Randomly pick one combination if there are two trophic species that have the same trophic similarity.
    if(dim(idmax)[1] > 1){
      idmax = idmax[sample(seq(1,dim(idmax)[1]),1),]
    }

  }

  # The next section of code combines the two trophic species listed in idmax

  # Take the biomass average of their feeding
  imat_mod[idmax,] = imat[idmax,]*prop$B[idmax]/sum(prop$B[idmax])

  imat_mod[,idmax[1]] = imat[,idmax[1]]*prop$B[idmax[1]]/sum(prop$B[idmax])

  imat_mod[,idmax[2]] = imat[,idmax[2]]*prop$B[idmax[2]]/sum(prop$B[idmax])

  # Combine them
  imat_mod_update[idmax[1],] = imat_mod[idmax[1],] + imat_mod[idmax[2],]
  imat_mod_update[,idmax[1]] = imat_mod[,idmax[1]] + imat_mod[,idmax[2]]

  # Remove old species
  imat_mod_update = imat_mod_update[-idmax[2],]
  imat_mod_update = imat_mod_update[,-idmax[2]]

  #Update name: The new name is name1/name2
  rownames(imat_mod_update)[idmax[1]] = colnames(imat_mod_update)[idmax[1]] = paste0(colnames(imat_mod)[idmax], collapse = "/")

  # Rebalance the feeding relationships to 1 if asked
  if(allFEEDING1){
    imat_mod_update[imat_mod_update > 0] = 1
  }

  # Modify the parameters: Weighted average based on biomass

  # necessary function that is used to add a factor to the property list of trophic species names
  addlevel <- function(x, newlevel){
    if(is.factor(x)) return(factor(x, levels=c(levels(x), newlevel)))
    return(x)
  }

  # Find the properties of the combined trophic species
  prop_mod = prop[prop$ID %in% colnames(imat_mod)[idmax],]

  # Create an empty property vector to fill with the new data
  emptyprop = prop[1,]

  # The biomass is the sum of the two combined species biomass
  emptyprop$B = sum(prop_mod$B)

  # The other values are the biomass weighted average
  emptyprop[,c("d", "a", "p", "CN")] = colSums(prop_mod[,c("d", "a", "p", "CN")]*prop_mod[,"B"]/sum(prop_mod$B))

  emptyprop[,"canIMM"] = max(prop_mod$canIMM)
  if(length(unique(prop_mod$canIMM))>1) warning("Combined trophic species have different parameters for canIMM. Choosing canIMM == 1")

  emptyprop[,"isDetritus"] = max(prop_mod$isDetritus)
  if(length(prop_mod$isDetritus > 0)<2) warning("Combined trophic species have different parameters for isDetritus Choosing isDetritus == 1")

  emptyprop[,"DetritusRecycling"] = sum(prop_mod$DetritusRecycling) # Adding the contributions of detritus if combining two detrial pools

  emptyprop[,"isPlant"] = max(prop_mod$isPlant)
  if(length(unique(prop_mod$isPlant))>1) warning("Combined trophic species have different parameters for isPlant Choosing isPlant == 1")

  # Create the property updated data frame
  prop_update = prop

  # Add the level of the ID column
  prop_update[,1] = addlevel(prop_update[,1], newlevel = paste0(colnames(imat_mod)[idmax], collapse = "/"))

  # Add in the new row for the new species to replace the first species that was combined
  prop_update[prop_update$ID %in% colnames(imat_mod)[idmax[1]],-1] = emptyprop[,-1]
  prop_update[prop_update$ID %in% colnames(imat_mod)[idmax[1]],1] = paste0(colnames(imat_mod)[idmax], collapse = "/")

  # Remove the second species that was combined
  prop_update = prop_update[!(prop_update$ID %in% colnames(imat_mod)[idmax[2]]),]

  # Set cannibalism to zero
  if(deleteCOMBOcannibal){
    imat_mod_update = as.matrix(imat_mod_update)
    diag(imat_mod_update)[which(rownames(imat_mod_update) %in% paste0(colnames(imat_mod)[idmax], collapse = "/"))] = 0
  }

  # Reset DetritusRecycling so it sums to one:
  if(sum(prop_update$DetritusRecycling) != 0){
    prop_update$DetritusRecycling = prop_update$DetritusRecycling/sum(prop_update$DetritusRecycling)
  }

  # Double check that imat and prop have the same order
  stopifnot(all(colnames(imat_mod_update) == prop_update$ID))

  return(list(imat = imat_mod_update, prop = prop_update))
}
