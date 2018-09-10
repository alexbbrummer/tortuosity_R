subsample <- function(vessel_coords, density = 2){
  ## This function will subsabmple vessel coordinates to a desired point density (default is 2 points per micrometer).  Currently this function will be designed specifically for the mouse brain datasets (units of micrometers).  Future iterations will need to be adapted for user defined densities with variable units (for example, centimeters).  Higher precision will be attained by starting with high density interpolation data.  This method functions by cumulatively summing nearest neightbor euclidean distances until reaching the desired density.  Indeces of desired densities are then stored and returned to the user.
  
  # Initialize pointwise_distance
  pointwise_distance <- mat.or.vec(length(vessel_coords[,1]), 1)
  # pointwise_distance[] <- NaN
  
  # Initialize subsample_indeces
  subsample_indeces <- c()
  
  # define arclength distance tolerance between points as inverse of density.
  arc_tol <- 1/density
  
  ## For testing/examing method, uncomment arc_length calls.
  # arc_length <- mat.or.vec(length(vessel_coords[,1]),1)
  # arc_length[] <- NaN
  
  # Loop through vessel coords
  for(i in 1:(length(pointwise_distance)-1)){
    pointwise_distance[i] <- euclidean_distance(vessel_coords[i,], vessel_coords[i+1,])
    # arc_length[i] <- sum(pointwise_distance)
  }
  
  arcsum <- 0
  for(i in 1:length(pointwise_distance)){
    if(arcsum < arc_tol){
      arcsum <- arcsum + pointwise_distance[i]
    }else if(arcsum >= arc_tol){
      arcsum <- 0
      subsample_indeces <- c(subsample_indeces, (i-1))
    }
  }
  return(subsample_indeces)
}