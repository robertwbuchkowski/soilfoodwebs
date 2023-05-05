#' A function to calculate carbon and nitrogen fluxes in the food web.
#'
#' @param  usin The community that you are analyzing: contains a matrix of interactions and a data frame of properties in a list.
#' @param mkplot Boolean: Should the plots be output?
#' @param whattoplot A vector of what to plot. Food web typology (web), nitrogen mineralization (Nmin), and/or carbon mineralization (Cmin).
#' @param showCN Boolean: # Should the food web show the C:N ratio of each trophic species next to it's name in the food web plot?
#' @param BOX.SIZE Size of boxes in the food web plot
#' @param BOX.PROP Proportion of box length and width in the food web plot
#' @param BOX.CEX Size of box text in the food web plot
#' @param PLOT.CEX Size of plot
#' @param edgepos Where to put the far left and far right boxes for trophic species in the food web plot (range = 0, 1).
#' @param TCK The size of the ticks on the C.min and N.min plots.
#' @param shuffleTL A Boolean stating whether the community should be sorted.
#' @param prettynames Alternative names in order for the food web plot trophic species. Cannot be used if showCN = T.
#' @param fwdlwdcust # A matrix of arrow line widths, same dimensions at the food web for plot customization.
#' @param arrowlog Boolean: Should relative arrow widths in the food web plot be on a log scale?
#' @param arrowsizerange The range of arrow sizes in the food web plot.
#' @param rmzeros A Boolean determining whether trophic species with zero biomass should be removed from the community before analysis.
#' @param eqmtolerance A value used to set the equilibrium tolerance for the food web verification. If NA, the default value used by the function all.equal is used, which is approximately 1.5e-8.
#' @return A list of consumption rates, carbon mineralization, nitrogen mineralization, carbon and nitrogen consumption rates, and the modified community if zeros where removed or sorting occurred.
#' @examples
#' comana(intro_comm)
#' @export

comana <- function(usin,
                   mkplot = FALSE,
                   whattoplot = c("web","Nmin","Cmin"),
                   showCN = FALSE,
                   BOX.SIZE = 0.1,
                   BOX.PROP = 0.3,
                   BOX.CEX = 1,
                   PLOT.CEX = 1,
                   edgepos = c(0.1, 0.9),
                   TCK = 0.05,
                   shuffleTL = FALSE,
                   prettynames = NA,
                   fwdlwdcust = NULL,
                   arrowlog = FALSE,
                   arrowsizerange = c(0.1,30),
                   rmzeros = TRUE,
                   eqmtolerance = NA
){

  # Check the community:
  usin = checkcomm(usin, shuffleTL = shuffleTL, rmzeros =rmzeros, verbose = FALSE)

  imat = usin$imat # row values of imat sets predator feeding preferences!
  prop = usin$prop # properties of each trophic species

  # Save the trophic levels for later use
  TL = TLcheddar(imat)
  TL2 = TL[order(-TL)]

  Nnodes = dim(imat)[1] # Number of nodes in the food web

  # Create a vector for the consumption rates
  temp_mat =
    -1*t(imat)*matrix(prop$B, nrow = Nnodes, ncol = Nnodes)/matrix(rowSums(imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow = T)), nrow = Nnodes, ncol = Nnodes, byrow = T)

  temp_mat[!is.finite(temp_mat)] = 0 # Replace non-finite values with 0 because total consumption was zero in this case

  diag(temp_mat) = prop$a*prop$p + diag(temp_mat) # Add in a*p term

  consumption = base::solve(temp_mat,(prop$d*prop$B))

  # Confirm that this solution is unique by showing Ax = 0 produces x = 0
  if(any(solve(temp_mat,rep(0, Nnodes)) != 0)){
    warning("Solution to the web is not unique!")
  }

  names(consumption) = colnames(imat) # Names match the trophic species names

  # Create a new matrix for feeding rates with the same dimensions as the food web matrix
  fmat = (imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)/rowSums(imat*matrix(prop$B, nrow = Nnodes, ncol = Nnodes, byrow=TRUE)))*consumption

  fmat[!is.finite(fmat)] = 0 # Replace NaN with 0 for no feeding

  # Fix the detritus calculations: detritus receives dead material from all other trophic levels and so consumption is the losses minus the inputs it already gets from inefficient eating and dead biomass. This value can be negative indicating that inputs from outside the ecosystem are necessary to meet internal demand for C.
  if(any(prop$isDetritus >0)){
    detritusPOS = which(prop$isDetritus >0)

    for(i in detritusPOS){
      consumption[i] = sum(fmat[,i]) - prop$DetritusRecycling[i]*(sum((1-prop$a)*consumption) + sum(prop$d*prop$B))
    }
  }

  # Calculate carbon mineralization using the production efficiency
  Cmin = prop$a*(1-prop$p)*consumption

  # Calculate nitrogen mineralization as a matrix so that the nitrogen mineralized by each consumption interaction is recorded. Summing the rows of this matrix gives the nitrogen that is mineralized overall.
  Nmin = (matrix(1/prop$CN, nrow = Nnodes, ncol = Nnodes, byrow = T) - matrix(prop$p/prop$CN, nrow=Nnodes, ncol = Nnodes))*matrix(prop$a, nrow=Nnodes, ncol = Nnodes)*fmat

  # Calculate the nitrogen flux throughout the matrix as the amount leaving a node including the N that will be mineralized.
  Nfmat = matrix(1/prop$CN, nrow = Nnodes, ncol = Nnodes, byrow = TRUE)*fmat

  # Make the plot if necessary
  if(mkplot){

    TL3 = round(TL2)
    # Position matrix
    posmat = matrix(NA, ncol = 2, nrow = Nnodes)
    posmat[,2] = RESCALE(TL2, a = 0.1, b=0.9)
    for(tl in sort(unique(TL3))){
      posmat[TL3 == tl,1] = seq(edgepos[1], edgepos[2], length = table(TL3)[as.character(tl)])
    }

    if(is.null(fwdlwdcust)){
      fmatlwd = fmat
      if(arrowlog) fmatlwd[fmatlwd!=0]  = log(fmatlwd[fmatlwd!=0])
      fmatlwd = RESCALE(fmatlwd, a = arrowsizerange[1], b = arrowsizerange[2])
    }else{
      fmatlwd = fwdlwdcust
    }

    if("web" %in% whattoplot){

      if(showCN){
        namevec = paste0(colnames(imat), ": ", prop$CN)
      }else{
        if(is.na(prettynames[1])){
          namevec = colnames(imat)
        }else{
          namevec = prettynames
        }

      }

      diagram::plotmat(imat, pos = posmat,name = namevec, box.size = BOX.SIZE, shadow.size = 0, box.type = "rect", box.prop = BOX.PROP, arr.lwd = fmatlwd, cex.txt = 0, arr.length = 0.1, box.cex = BOX.CEX)
    }


    if("Nmin" %in% whattoplot){
      graphics::barplot(rowSums(Nmin), ylab = "N min.", names.arg = prop$ID)
    }

    if("Cmin" %in% whattoplot){
      graphics::barplot(Cmin, ylab = "C min.", names.arg = prop$ID)
    }


  }

  # Create a series of carbon and nitrogen outputs with consistent names and formats by changing Nmin to Nminmat and creating vector of Nmin:
  Nminmat = Nmin
  Nmin = rowSums(Nmin)

  whattoreturn = list(consumption = consumption, # consumption vector
                      Cmin = Cmin, # carbon mineralization vector
                      Nmin = Nmin, # nitrogen mineralization vector
                      Nminmat = Nminmat, # nitrogen mineralization matrix to show contributions of each trophic species to mineralization
                      fmat = fmat, # Feeding matrix
                      Nfmat = Nfmat, # Nitrogen feeding matrix negative values are turned to zero and indicate too much cannibalism for equilibrium
                      usin =list(imat = imat, prop = prop)) # revised community with sorted trophic levels

  if(!checkeqm(whattoreturn,eqmtolerance = eqmtolerance)){
    stop("Equilibrium not successfully reached at tolerance provided. Check the calculations or increase the tolerance.")
  }

  if(any(consumption[TLcheddar(imat) != 1] < 0)){
    warning("Negative consumption is being predictd. This web probably does not have a stable solution.")
  }

  return(whattoreturn)
}
