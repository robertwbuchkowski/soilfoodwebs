#' A function to combine trophic species.
#'
#' @param usin The community were you are combining trophic species
#' @param selected Select trophic species to combine, which ones do you want to combine (vector of names)? If left as NA, the two most similar trophic species are combined. Similarity is determined by shared feeding relationships.
#' @param  deleteCOMBOcannibal Boolean: Do you want to delete the cannibalism that may have been created by combining two trophic species (T) or leave it in the model (F)?
#' @param allFEEDING1 Boolean: Do you want to return all feeding preferences to 1 (T), or would you like to set the feeding preferences of the newly combined trophic species as the biomass-weighted average of the old ones (F)?
#' @param newname The name you want to use for the combined trophic species. Default replaces combines the names of the original trophic species divided by "/".
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
                     allFEEDING1 = F,
                     newname = NA
){
  imat = usin$imat
  prop = usin$prop
  Nnodes = dim(imat)[1]

  imat_mod = imat # Create some replicates of the feeding matrix to store intermediate steps

  # If seltf trophic species are combined, then just set them in the idmax, which is the matrix identifying their position in the web
  if(!any(is.na(selected))){

    # Check to make sure all selected nodes are in the web
    if(!all(selected %in% prop$ID)){
      stop("All selected species must be in the food web.")
    }

    idmax = which(rownames(imat) %in% selected)
  }else{

    # If you don't select the species to be combined, then the most similar species are calculated by comparing the number of matching trophic links

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
    # Clean up the names. Product is the row/column numbers of the species to be combined.
    idmax = unname(idmax)
  }

  #Update name: The new name is name1/name2
  if(!is.na(newname)){
    newname1 = newname
  }else{
    newname1 = paste0(colnames(imat)[idmax], collapse = "/")
  }

  # Modify the parameters: Weighted average based on biomass

  # necessary function that is used to add a factor to the property list of trophic species names
  addlevel <- function(x, newlevel){
    if(is.factor(x)) return(factor(x, levels=c(levels(x), newlevel)))
    return(x)
  }

  # Find the properties of the combined trophic species
  prop_mod = prop[prop$ID %in% colnames(imat)[idmax],]

  # Create an empty property vector to fill with the new data
  emptyprop = prop[1,]

  # The biomass is the sum of the two combined species biomass
  emptyprop$B = sum(prop_mod$B)

  # The other values are the biomass weighted average
  emptyprop[,c("d", "a", "p", "CN")] = colSums(prop_mod[,c("d", "a", "p", "CN")]*prop_mod[,"B"]/sum(prop_mod$B))

  emptyprop[,"canIMM"] = max(prop_mod$canIMM)
  if(length(unique(prop_mod$canIMM))>1) warning("Combined trophic species have different parameters for canIMM. Choosing canIMM == 1")

  emptyprop[,"isDetritus"] = max(prop_mod$isDetritus)
  if(length(unique(prop_mod$isDetritus)) > 1) warning("Combined trophic species have different parameters for isDetritus. Choosing isDetritus == 1")

  emptyprop[,"DetritusRecycling"] = sum(prop_mod$DetritusRecycling) # Adding the contributions of detritus if combining detritus pools

  emptyprop[,"isPlant"] = max(prop_mod$isPlant)
  if(length(unique(prop_mod$isPlant))>1) warning("Combined trophic species have different parameters for isPlant Choosing isPlant == 1")

  # Create the property updated data frame
  prop_update = prop

  # Add the level of the ID column
  prop_update[,1] = addlevel(prop_update[,1], newlevel = newname1)

  # Add in the new row for the new species to replace the first species that was combined
  prop_update[idmax[1],] = emptyprop
  prop_update[idmax[1],"ID"] = newname1

  # Remove all other species
  prop_update = prop_update[-idmax[-1],]


  # Reset DetritusRecycling so it sums to one:
  if(sum(prop_update$DetritusRecycling) != 0){
    prop_update$DetritusRecycling = prop_update$DetritusRecycling/sum(prop_update$DetritusRecycling)
  }

  # The next section of code combines the trophic species listed in idmax in the imat

  # Calculate diet proportions in the original model:

  dietprop = (imat * matrix(prop$B, nrow = Nnodes, ncol = Nnodes,
                            byrow = T)/rowSums(imat * matrix(prop$B, nrow = Nnodes,
                                                             ncol = Nnodes, byrow = T)))

  # Set NA as zero
  dietprop[is.na(dietprop)] = 0

  dietprop2 = dietprop
  # Combine the trophic species: Predators of them have their diet portions summed, the combined trophic species have prey diet portions summed
  # Result is saved in the row and column of the first species to be combined
  dietprop2[,idmax[1]] = rowSums(dietprop[,idmax])
  dietprop2[idmax[1],] = colSums(dietprop[idmax,])/sum(dietprop[idmax,])

  # Take any feeding within the combined trophic species and reapply it
  dietprop2[idmax[1],idmax[1]] = sum(dietprop2[idmax[1],idmax])

  # Remove combined trophic species (leave the one that received the new parameters)
  dietprop2 = dietprop2[-idmax[-1],]
  dietprop2 = dietprop2[,-idmax[-1]]

  # Update the name:
  rownames(dietprop2)[idmax[1]] = colnames(dietprop2)[idmax[1]] = newname1


  # Make sure diet proportions are still summing to 1:
  # Select only non-zero row sums
  checkrowsums = rowSums(dietprop2)
  checkrowsums = checkrowsums[checkrowsums != 0]

  if(!all.equal(unname(checkrowsums), rep(1, length(checkrowsums)))){
    stop("Error calculating diet proportions. Check the math unless allFEEDING1 is true.")
  }

  # Convert diet proportions back into feeding preferences with new properties matrix
  imat_update = dietprop2/matrix(prop_update$B, nrow = dim(prop_update)[1], ncol = dim(prop_update)[1],byrow = T)


  # Confirm that the diet proportions turned out right
  checkdietprop = (imat_update * matrix(prop_update$B, nrow = dim(prop_update)[1], ncol = dim(prop_update)[1],byrow = T)/rowSums(imat_update * matrix(prop_update$B, nrow = dim(prop_update)[1],ncol = dim(prop_update)[1], byrow = T)))
  checkdietprop[is.na(checkdietprop)] = 0

  if(!all.equal(checkdietprop, dietprop2)){
    warning("Diet proportion NOT maintained by combining trophic species. Check the math unless allFEEDING1 is true.")
  }

  # Rebalance the feeding relationships to 1 if asked
  if(allFEEDING1){
    imat_update[imat_update > 0] = 1
  }

  # Set cannibalism to zero
  if(deleteCOMBOcannibal){
    imat_update = as.matrix(imat_update)
    diag(imat_update)[which(rownames(imat_update) %in% newname1)] = 0
  }

  # Double check that imat and prop have the same order
  if(!all(colnames(imat_update) == prop_update$ID)){
    stop("New species names not lining up. Check calculations.")
  }

  usin_temp = list(imat = imat_update,
                   prop = prop_update)

  return(checkcomm(usin_temp))
}
