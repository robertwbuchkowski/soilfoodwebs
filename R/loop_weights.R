#' A function to identify a loop in a food web vector
#'
#' @param VEC The vector to check.
#' @return Boolean
is_loop <- function(VEC){
  VEC[1] == VEC[length(VEC)]
}

#' A function to calculate the geometric mean
#'
#' @param VEC The vector whose mean you want to calculate.
#' @return The geometric mean.
geoMEAN <- function(VEC){
  exp(mean(log(VEC)))
}

#' A function to calculate the loop weights for one node.
#'
#' @param mat The matrix for which to calculate loop weights.
#' @param node The node for which to calculate loop weights.
#' @param max_length The maximum length of paths to check.
#' @return The loop weights for the chosen node.
path_n_weights <- function(mat, node=1, max_length = 100){

  Nnodes = dim(mat)[1]
  final_vec <- list()

  cur_node = node
  links <- which(mat[cur_node,]>0)

  # Add in the first connection
  vec <- vector(mode = "list", length = length(links))
  for(i in 1:length(links)){ vec[[i]] = c(cur_node, links[i])}

  numreps = 0
  while(T){
    # Add in the next connection
    for(i in 1:length(vec)){
      cur_path = vec[[i]]
      cur_node = utils::tail(cur_path, 1)
      links <- which(mat[cur_node,]>0) # unique(c(which(mat[cur_node,]>0), which(mat[,cur_node]>0)))
      vec[[i]] <- vector(mode = "list", length = length(links))
      if(length(links) > 0) for(j in 1:length(links)){ vec[[i]][[j]] = c(cur_path, links[j])}
    }

    # Flatten list into a single vector
    vec = rlang::flatten(vec)

    # Add the loops to the final vector
    final_vec = rlang::flatten(list(final_vec, vec[unlist(lapply(vec, is_loop))]))

    # Keep the non-loops
    vec_loops = unlist(lapply(vec, is_loop))

    if(length(vec_loops)>0){
      vec = vec[!vec_loops]

      # Get rid of the repeated steps
      vec = vec[unlist(lapply(vec, function(X) !any(table(X)>1)))]
    }else{
      vec = c()
    }

    numreps = numreps + 1

    if(length(vec) == 0 | numreps == max_length) break
  }

  # Get the path weights
  weights <- vector(mode = "list", length = length(final_vec))
  for(k in 1:length(weights)){
    cur_path2 <- final_vec[[k]]
    weights[[k]] <- rep(NA, length(cur_path2)-1)
    for(h in 1:length(weights[[k]])){
      weights[[k]][h] = mat[cur_path2[h],cur_path2[h+1]]
    }
  }
  return(list(path = final_vec, weights = weights))
}

#' A function to calculate loop weights for a matrix
#'
#' @param mat The matrix for which to calculate loop weights.
#' @param max_length The maximum length of paths to check.
#' @param max_return The number of loop weights per node you want to return.
#' @return The geometric mean.
#' @details
#' This function calculates the loop weights for each of possible loop up to the length max_length. This calculation can take a long time for large matricies, so a counter is provided. The function returns on the top 10 loop weights.
#' @export
loop_weights <- function(mat, max_length = 100, max_return = 10){
  results <- vector(mode = "list", length = dim(mat)[1])
  for(ij in 1:length(results)){
    output <- path_n_weights(mat,ij, max_length)
    output2 = data.frame(path = unlist(lapply(output$path, paste0, collapse = "-")),
                         loopweight = unlist(lapply(output$weights, geoMEAN)),
                         looplength = unlist(lapply(output$path, length)))

    output2 = output2[order((output2$loopweight*-1)),]

    results[[ij]] = output2[(1:max_return),]

    print(paste("Done", ij, "of", length(results)))

  }
  results <- do.call("rbind", results)
  results <- subset(results, results$looplength > 3)
  return(results)
}
