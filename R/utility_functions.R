#' Calculates the trophic level for each tropospecies
#'
#' @param W A matrix of trophic species. Rows eat columns.
#' @return A vector of trophic level assignments. The base of the food chain is 0.
#' @details
#' This function is a subset of code provided in the package Cheddar written by Lawrence Hudson, Dan Reuman and Rob Emerson. It is licensed under a BSD_2_clause, the text of which is provided as a comment in the function code. The original package can be found on CRAN or at https://github.com/quicklizard99/cheddar/
#' @export
#' @examples
#' TLcheddar(intro_comm$imat)
TLcheddar <- function(W){

  # Original license for this code:
  # YEAR:2014-2016
  # COPYRIGHT HOLDER:  Lawrence Hudson, Dan Reuman and Rob Emerson
  #
  #   and specify as
  #
  # License: BSD_2_clause + file LICENSE
  #
  # Copyright (c) 2014-2016, Lawrence Hudson, Dan Reuman and Rob Emerson
  #
  #   Redistribution and use in source and binary forms, with or without
  # modification, are permitted provided that the following conditions are
  # met:
  #
  #   Redistributions of source code must retain the above copyright
  # notice, this list of conditions and the following disclaimer.
  #
  # Redistributions in binary form must reproduce the above copyright
  # notice, this list of conditions and the following disclaimer in
  # the documentation and/or other materials provided with the
  # distribution.
  #
  # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
  # "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
  # LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
  # A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
  # HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
  # SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
  # LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
  # DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
  # THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
  # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


  # Make each row sum to one.
  rs <- rowSums(W)
  W <- W / matrix(rs, ncol=ncol(W), nrow=nrow(W))
  W[0==rs,] <- 0      # Fix NA resulting from div by zero

  # Identity matrix
  I <- diag(ncol(W))

  # Invert the matrix. From Rich 2011-08-17: "If this matrix inversion
  # fails, there is an important problem with the network topology.
  # For a food web to be energetically feasible, every node must be
  # connected to a basal node. When the inversion fails it is because
  # there is at least one node that has no connection to a basal node."
  result <- tryCatch(solve(I-W), error = function(e) e)

  if('error' %in% class(result))
  {
    tl <- rep(NA, ncol(W))
    names(tl) <- colnames(W)
  }
  else
  {
    # Returned vector is named by node
    tl <- rowSums(result)   # Isolated nodes have tl of 1
  }

  return (tl)
}

#' Sorts the trophic levels from lowest to highest
#'
#' @param usin The community to be sorted.
#' @return The community returned after sorting
#' @examples
#' TLsort(intro_comm)
#' @export
TLsort <- function(usin){

  # Identify the cannibals and mutual predators if not already done
  if(!("MutualPred" %in% colnames(usin$prop))){
    usin = can_mutfeed(usin)
  }

  imat = usin$imat # row values of imat sets predator feeding preferences!
  prop = usin$prop

  Nnodes = dim(imat)[1] # Number of nodes in the food web

  # Calculate the trophic level
  TL = seq(Nnodes,1,-1)
  names(TL) = colnames(imat)

  TL_temp = TLcheddar(imat)

  if(!(any(is.na(TL_temp)))){
    TL = TL_temp
  }else{
    print("Error in TL sorting: usually comes from a poorly defined matrix where feeding relationships are not possible at equilibrium.")
  }

  rm(TL_temp)

  # Shuffle the trophic levels

  # First sort by trophic level to get the order close
  imat = imat[order(-TL),order(-TL)]

  # Then confirm that there are no cases where a low trophic level predator eats a higher trophic level prey
  positions <- 1:dim(imat)[1]
  rid = min(positions)
  stuckinloop = 0
  while(rid <= max(positions)){

    # Build a list of predators
    predatorlist = imat[,rid] > 0
    # Turn any predators that are on the mutual feeding list to false...these can't be sorted hierarchically
    predatorlist[names(predatorlist) %in%  stringr::str_split(prop$MutualPred[rid], "/")[[1]]] = FALSE
    if(sum(predatorlist) == 0){
      rid = rid + 1
    }else{
      max_predator = max(which(predatorlist))
      if(max_predator > rid){
        part1 = min(positions):max_predator
        if(max_predator == max(positions)){
          part2 = NULL
        }else{
          part2 = (max_predator+1):max(positions)
        }
        new_order = c(part1[part1 != rid], rid, part2[part2 != rid])
        imat = imat[new_order,new_order]
        stuckinloop = stuckinloop + 1
        if(stuckinloop > 100){
          break
          warning("Fine level sorting not converging, returning early. Likely caused by extensive mutual feeding.")
        }
      }else{
        rid = rid + 1
        stuckinloop = 0
      }
    }
  }

  # Order of the colnames:
  prop = prop[match(colnames(imat), prop$ID),]

  stopifnot(all(colnames(imat) == prop$ID))

  return(list(imat = imat, prop = prop))
}

#' A function that identifies cannibalism and mutual feeding
#'
#' @param usin The community in which to identify cannibalism and mutual feeding
#' @return The community with mutual feeding added to the properties database.
#' @examples
#' intro_comm_mod = intro_comm
#' intro_comm_mod$imat["Predator", "Predator"] = 1
#' can_mutfeed(intro_comm_mod)
#' @export
can_mutfeed <- function(usin){
  imat = usin$imat # row values of imat sets predator feeding preferences!
  prop = usin$prop

  Nnodes = dim(imat)[1] # Number of nodes in the food web

  MutualPred = rep(NA,dim(imat)[1])
  for(i in 1:dim(imat)[1]){
    MutualPred[i] = paste0(names(which(imat[imat[i,] > 0,i]>0)), collapse = "/")

  }

  MutualPred[MutualPred == ""] = NA

  prop[,"MutualPred"] = MutualPred

  return(list(imat = imat, prop = prop))
}


#' A function to calculate the nitrogen surplus or deficit each species gets from consuming another species
#'
#' @param usin The community on which nitrogen calculations are made.
#' @return A matrix of nitrogen surpluses or deficits.
#' @examples
#' Aijfcn(intro_comm)
#' @export
Aijfcn <- function(usin){
  tt = comana(usin, shuffleTL = FALSE)
  FMAT = tt$fmat
  Nnodes = dim(FMAT)[1]
  FMAT2 = FMAT
  FMAT2[FMAT2>0] = 1

  # Below matricies are multipled to calculate the net excess or deficit of nitrogen in the organisms diet using their C:N raito and the C:N ratio of the food that they consume.
  Aij = (matrix(1/usin$prop$CN, nrow = Nnodes, ncol = Nnodes, byrow = TRUE) - matrix(usin$prop$p/usin$prop$CN, nrow=Nnodes, ncol = Nnodes))*matrix(usin$prop$a, nrow=Nnodes, ncol = Nnodes)*FMAT2

  return(Aij)
}

#' A utility function to calculate the consumption rate of each species on all prey assuming a type I functional response.
#'
#' @param usin The community on which feeding rate calculations are made.
#' @param shuffleTL A Boolean stating whether the community should be sorted before consumption rates are calculated.
#' @param rmzeros A Boolean determining whether trophic species with zero biomass should be removed from the community.
#' @return A matrix of consumption rates with units set by the the biomass input units in biomass and time.
#' @examples
#' Cijfcn(intro_comm)
#' @export
Cijfcn <- function(usin, shuffleTL = FALSE, rmzeros = TRUE){ # Function only requires the community inputs

  # Check the community
  usin = checkcomm(usin, shuffleTL = shuffleTL, rmzeros = rmzeros)

  imat = usin$imat # Feeding matrix
  prop = usin$prop # Properties of the species

  Nnodes = dim(imat)[1] # Number of nodes

  Bpred = matrix(prop$B, ncol = Nnodes, nrow = Nnodes) # A matrix of predators
  Bprey = matrix(prop$B, ncol = Nnodes, nrow = Nnodes, byrow = T) # A matrix of prey
  fmat = comana(usin, shuffleTL = FALSE)$fmat # Get the consumption matrix (units = gC/ time)
  cij = fmat/(Bpred*Bprey) # Get the consumption rate matrix (units 1/ (gC * time))
  return(cij) # Return the consumption rates (gC^-1 time^-1)
}

#' A function check the community for errors before it is used in calculations.
#'
#' @param usin The community to check.
#' @param shuffleTL A Boolean stating whether the community should be sorted.
#' @param rmzeros A Boolean determining whether trophic species with zero biomass should be removed from the community.
#' @param verbose A Boolean. Do you want to see all the warnings?
#' @return The checked community.
#' @examples
#' checkcomm(intro_comm)
#' @export
checkcomm <- function(usin, shuffleTL = FALSE, rmzeros = TRUE, verbose = TRUE){

  # Check that all properties are included:
  if(!all(c("ID", "d", "a", "p", "B", "CN", "isDetritus", "isPlant", "canIMM","DetritusRecycling") %in% colnames(usin$prop))) stop("The community needs all the following properties in the database: ID, d, a, p, B, CN, isDetritus, isPlant, canIMM, and DetritusRecycling. MutualPred column will be created if it is not included.")

  # Remove zeros from the community
  if(rmzeros & any(usin$prop$B == 0)){

    usin$imat = usin$imat[usin$prop$B != 0,usin$prop$B != 0]

    usin$prop = subset(usin$prop, usin$prop$B !=0)

    if(verbose) warning("Removing zero biomass nodes.")
  }

  # Sort by trophic level
  if(shuffleTL) usin = TLsort(usin)

  # Calculate trophic level
  TL = TLcheddar(usin$imat)

  # Check for several errors in the community and return to user if any is true:

  if(!is.matrix(usin$imat)) stop("Feeding matrix must be of class matrix.")

  if(any(usin$prop$isDetritus >0)){
    # Make sure TL of detritus is 1
    if(!all(TL[usin$prop$isDetritus >0] ==1)) stop("Detirtus must have a trophic level position of 1.")


    # Detritus a and p are 1 and death is zero.
    if(!all(
      subset(usin$prop, usin$prop$isDetritus >0)$a == 1 & subset(usin$prop, usin$prop$isDetritus >0)$p == 1 & subset(usin$prop, usin$prop$isDetritus >0)$d == 0
    )){
      stop("Detritus must have a = 1, p = 1, and d = 0.")
    }

    # Make sure detritus recycling proportion sums to 1
    if(sum(usin$prop$DetritusRecycling) != 1 & verbose){
      warning("Rescaling Detritus Recycling to sum to 1.")

      usin$prop$DetritusRecycling = usin$prop$DetritusRecycling/sum(usin$prop$DetritusRecycling)
    }
  }

  # Properties and matrix have the same individuals
  if(!all(c(colnames(usin$imat) == rownames(usin$imat),colnames(usin$imat) == usin$prop$ID))) stop("Column names, row names, and property data frame IDs must all be the same and in the same order.")

  # Identify the cannibals and mutual predators after a reset in case it has changed since community was created.
  usin$prop$MutualPred = NULL
  if(!("MutualPred" %in% colnames(usin$prop))){
    usin = can_mutfeed(usin)
  }
  return(usin)
}

#' A function to rescale a vector.
#'
#' @param invec The vector to re-scale.
#' @param a The lower limit of the new scale.
#' @param b The upper limit of the new scale.
#' @return The scaled vector.
RESCALE <- function(invec, a = 0, b = 1){
  (b-a)*(invec - min(invec))/(max(invec)-min(invec)) + a
}

#' Remove nodes from community.
#'
#' @param COMM The community from which to remove nodes.
#' @param toremove A vector of nodes to remove using their names.
#' @return The community without the removed nodes.
#' @examples
#' removenodes(intro_comm, c("Predator"))
#' @export
removenodes <- function(COMM, toremove){
  whichtorm = !(COMM$prop$ID %in% toremove)

  COMM$imat = COMM$imat[whichtorm,whichtorm]

  COMM$prop = subset(COMM$prop, !(COMM$prop$ID %in% toremove))

  stopifnot(all(COMM$prop$ID == rownames(COMM$imat)))
  stopifnot(all(COMM$prop$ID == colnames(COMM$imat)))

  return(COMM)
}

#' Add node to the community
#'
#' @param COMM The community to which to add nodes.
#' @param newname The new node ID.
#' @param prey A vector of prey preferences with names.
#' @param predator A vector of predators and their preferences with name.
#' @param newprops A vector of the new properties with the appropriate names.
#' @return The community with the new node.
#' @seealso \code{\link{removenodes}}
#' @examples
#' # Add a node to the introductory community:
#' newnode(intro_comm, "NewNode",
#' prey = c(Detritus1 = 1),
#' predator = c(Predator = 2, Orib1 = 0.1),
#' newprops = c(d = 1, a = 0.1, p = 0.1,
#' B = 10, CN = 10, DetritusRecycling = 0,
#' isDetritus = 0, isPlant = 0, canIMM = 0))
#' @export
newnode <- function(COMM, newname, prey = NA, predator = NA, newprops){

  Nnodes = dim(COMM$imat)[1]

  # Add the level of the ID column
  COMM$prop$ID = factor(COMM$prop$ID, levels=c(COMM$prop$ID, newname))

  # Get rid of MutualPred column for now. It will be re-added to the properties when the community is checked below:
  if("MutualPred" %in% colnames(COMM$prop)) COMM$prop$MutualPred = NULL

  # Add the properties to the properties database
  COMM$prop = rbind(COMM$prop,
                    data.frame(ID = newname,
                               d = newprops["d"],
                               a = newprops["a"],
                               p = newprops["p"],
                               B = newprops["B"],
                               CN = newprops["CN"],
                               DetritusRecycling = newprops["DetritusRecycling"],
                               isDetritus = newprops["isDetritus"],
                               isPlant = newprops["isPlant"],
                               canIMM = newprops["canIMM"]))

  # Add the new species to the matrix

  COMM$imat = rbind(cbind(COMM$imat,rep(0, Nnodes)), rep(0, Nnodes + 1))
  row.names(COMM$imat)[(Nnodes +1)] = newname
  colnames(COMM$imat)[(Nnodes +1)] = newname

  # Identify prey
  if(!all(is.na(prey))){
    for(ii in 1:length(prey)){
      COMM$imat[newname, names(prey)[ii]] = unname(prey[ii])
    }
  }

  # Identify predators
  if(!all(is.na(predator))){
    for(ii in 1:length(predator)){
      COMM$imat[names(predator)[ii],newname] = unname(predator[ii])
    }
  }
  COMM = checkcomm(COMM)
  return(COMM)
}

#' Rename a node in a  community.
#'
#' @param COMM The community from which to remove nodes.
#' @param oldname The node's old name
#' @param newname The node's new name
#' @return The community with the new name.
#' @examples
#' renamenode(intro_comm, oldname = "Predator", newname = "NewPredator")
#' @export
renamenode <- function(COMM, oldname,newname){
  whichtorm = COMM$prop$ID %in% oldname

  colnames(COMM$imat)[whichtorm] = rownames(COMM$imat)[whichtorm] = newname

  COMM$prop$ID[whichtorm] = newname

  stopifnot(all(COMM$prop$ID == rownames(COMM$imat)))
  stopifnot(all(COMM$prop$ID == colnames(COMM$imat)))

  COMM = checkcomm(COMM)

  return(COMM)
}

#' Check the carbon flux equilibrium output by comana.
#'
#' @param cares The output from comana.
#' @param eqmtolerance A value used to set the equilibrium tolerance for the food web verification. If NA, the default value used by the function all.equal is used, which is approximately 1.5e-8.
#' @return Boolean: Is the community at equilibrium with inputs to the first trophic level?
#' @examples
#' checkeqm(comana(intro_comm))
#' @export
checkeqm <- function(cares, eqmtolerance = NA){

  netchange = # Consumption rate
    cares$usin$prop$a*cares$usin$prop$p*rowSums(cares$fmat) -

    # Natural death
    cares$usin$prop$d*cares$usin$prop$B -

    # Predation
    colSums(cares$fmat) +

    # Detritus recycling
    cares$usin$prop$DetritusRecycling*sum(cares$usin$prop$d*cares$usin$prop$B + (1 - cares$usin$prop$a)*cares$usin$prop$p*rowSums(cares$fmat))

  # Get rid of basal trophic levels where inputs are OK
  netchange = netchange[TLcheddar(cares$usin$imat) != 1]

  if(is.na(eqmtolerance)){
    return(all.equal(unname(netchange),rep(0, length(netchange))))
  }else{
    return(all.equal(unname(netchange),rep(0, length(netchange)), tolerance = eqmtolerance))
  }
}
