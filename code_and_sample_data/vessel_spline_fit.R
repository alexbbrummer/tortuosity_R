vessel_spline_fit <- function(vessel, number_samples = 10000, spline = "pspline", m = c(2,2), smoothing = NULL, subsample_density = 2, plot = TRUE){
  # This function performs a spline smoothing fit of the provided voxels for a given vessel.  The vessel parameter needs to be an N X 3 matrix of the (x, y, z) coordinates of the vessel voxels.  The number_samples parameter determines how many new points will be formed for defining the output vessel.  Currently we just define it as a single value, but in the future this should be modified to be representative of samples per unit arclength to maintain a constant number of samples through curves.  The smoothing parameter is a place holder for now that will be used to toggle the extent of smoothing.  This can be more complicated than just a single number.  See ?smooth.spline for more information on the use of the smoothing parameter, which in the documentation is referred to as "spar".  The subsample_density is specific to a call to the function subsample, where interpolated vessel coordinates are subsampled to a desired point density. The default is 2 points per micrometer, where micrometers are chosen to adhere to the current analysis of the mouse brain datasets.
  
  library("pracma")
  library("mgcv")
  library("tictoc")
  
  # source("./single_vessel_plotter.R")
  
  tic("total")
  
  vessel_steps <- nrow(vessel)
  
  if(plot){
    tic("plot points")
    nopen3d()
    single_vessel_plotter(vessel_coords = vessel)
    toc()
  }
  
  # Here is the Pspline method for performing the smoothing.  Following the literature of Auiguiler et al., "Comparative study of different B-spline approaches for functional data", Mathematical and Computer Modelling, 58 (2013) 1568-1579, we were motivated to experiment with this method.  Current code implementation in R for this method makes it easily scaleable for AIC testing to automate selection of "the most ideal fit".
  
  if(spline == "pspline"){
    tic("psplinefit")
    ## Add parametric variable to vessel data.frame
    parametric <- 1:vessel_steps
    vessel <- data.frame(vessel, "parametric" = parametric)
    
    ## Initialize gam fitting lists/vectors for scanning AIC scores as function of number of knots.
    fitxgamlist <- list()
    fitxaicscore <- mat.or.vec((vessel_steps-4), 1)
    
    fitygamlist <- list()
    fityaicscore <- mat.or.vec((vessel_steps-4), 1)
    
    fitzgamlist <- list()
    fitzaicscore <- mat.or.vec((vessel_steps-4), 1)
    
    ## Find ideal number of knots based on sign of slope of AIC values using centered difference method.  This approach is based on that taken by David Ruppert (2002) Selecting the Number of Knots for Penalized Splines, Journal of Computational and Graphical Statistics, 11:4, 735-757, DOI: 10.1198/106186002853. However, we base our method on an AIC slope comparison method instead of that taken by Ruppert.  Additionally, we use equally spaced knot locations justified by the fact that our original data are nearly equally spaced to begin with, a decision justitied by the work of E. T. Y. Lee (1989) Choosing nodes in parametric curve interpolation, Computer-aided design, 21:6 363-370.
    
    ## Loop through range of knots. Variable k essentially defines number of equally spaced knots (it is also counting number of model variables (1) + degree - 1 of splines). 
    for(i in 1:2){
      fitxgamlist[[i]] <- gam(formula = x ~ s(parametric, bs = "ps", k = (i+m[1]+1), fx = FALSE, m = m), data = vessel)
      fitxaicscore[i] <- AIC(fitxgamlist[[i]])
      
      fitygamlist[[i]] <- gam(formula = y ~ s(parametric, bs = "ps", k = (i+m[1]+1), fx = FALSE, m = m), data = vessel)
      fityaicscore[i] <- AIC(fitygamlist[[i]])
      
      fitzgamlist[[i]] <- gam(formula = z ~ s(parametric, bs = "ps", k = (i+m[1]+1), fx = FALSE, m = m), data = vessel)
      fitzaicscore[i] <- AIC(fitzgamlist[[i]])
    }
    
    xknots <- 1
    for(i in 2:(vessel_steps-4)){
      fitxgamlist[[i+1]] <- gam(formula = x ~ s(parametric, bs = "ps", k = (i+1+m[1]+1), fx = FALSE, m = m), data = vessel)
      fitxaicscore[i+1] <- AIC(fitxgamlist[[i+1]])
      local_fitxaic_slope <- (fitxaicscore[i+1] - fitxaicscore[i-1])/2
      # print(local_fitxaic_slope)
      if(local_fitxaic_slope >= 0){
        # print(i-1)
        xknots <- i-1
        break
      }
    }
    
    yknots <- 1
    for(i in 2:(vessel_steps-4)){ 
      fitygamlist[[i+1]] <- gam(formula = y ~ s(parametric, bs = "ps", k = (i+1+m[1]+1), fx = FALSE, m = m), data = vessel)
      fityaicscore[i+1] <- AIC(fitygamlist[[i+1]])
      local_fityaic_slope <- (fityaicscore[i+1] - fityaicscore[i-1])/2
      # print(local_fityaic_slope)
      if(local_fityaic_slope >= 0){
        # print(i-1)
        yknots <- i-1
        break
      }
    }
    
    zknots <- 1
    for(i in 2:(vessel_steps-4)){ 
      fitzgamlist[[i+1]] <- gam(formula = z ~ s(parametric, bs = "ps", k = (i+1+m[1]+1), fx = FALSE, m = m), data = vessel)
      fitzaicscore[i+1] <- AIC(fitzgamlist[[i+1]])
      local_fitzaic_slope <- (fitzaicscore[i+1] - fitzaicscore[i-1])/2
      # print(local_fitzaic_slope)
      if(local_fitzaic_slope >= 0){
        # print(i-1)
        zknots <- i-1
        break
      }
    }
    toc()
    ## Build 3D interpolated vessel for plotting.
    tic("interpolating")
    fit_length <- number_samples
    newdata <- data.frame(parametric = seq(from = 0, to = vessel_steps, length.out = fit_length))
    x <- predict.gam(fitxgamlist[[xknots]], newdata)
    y <- predict.gam(fitygamlist[[yknots]], newdata)
    z <- predict.gam(fitzgamlist[[zknots]], newdata)
    
    smth_vessel <- cbind(x, y, z)
    toc()
    
    ## Call to subsample function to subsample points to desired density.
    if(!is.null(subsample_density)){
      tic("subsampling")
      subsample_indeces <- subsample(vessel_coords = smth_vessel, density = subsample_density)
      smth_vessel <- smth_vessel[subsample_indeces,]
      toc()
    }
    
    ## Trim endpoints
    
    ## Call to single_vessel_ploter to plot smth_vessel.
    if(plot){
      tic("plot pspline")
      single_vessel_plotter(vessel_coords = smth_vessel, new = TRUE, col = "purple")
      toc()
    }
  }
  
  ## Call to alternate smoothing methods.
  else if(spline == "smoothing"){
    if(!is.null(smoothing)){
      fitx <- smooth.spline(x = c(1:vessel_steps), y = vessel[,1], spar = smoothing)
      fity <- smooth.spline(x = c(1:vessel_steps), y = vessel[,2], spar = smoothing)
      fitz <- smooth.spline(x = c(1:vessel_steps), y = vessel[,3], spar = smoothing)
      
      # Here is where we interpolate additional points to yield (ideally) more contiunuous measurements of vessel coordinates.
      fit_length <- number_samples
      
      x <- predict(fitx,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      y <- predict(fity,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      z <- predict(fitz,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      
      smth_vessel <- cbind(x, y, z)
      
      if(plot){
        single_vessel_plotter(vessel_coords = smth_vessel, new = TRUE, col = "purple")
      }
      
    }else{
      fitx <- smooth.spline(x = c(1:vessel_steps), y = vessel[,1], cv = TRUE)
      fity <- smooth.spline(x = c(1:vessel_steps), y = vessel[,2], cv = TRUE)
      fitz <- smooth.spline(x = c(1:vessel_steps), y = vessel[,3], cv = TRUE)
      
      # Here is where we interpolate additional points to yield (ideally) more contiunuous measurements of vessel coordinates.
      fit_length <- number_samples
      
      x <- predict(fitx,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      y <- predict(fity,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      z <- predict(fitz,seq(from = 1, to = vessel_steps, length.out = fit_length))$y
      
      smth_vessel <- cbind(x, y, z)
      
      if(plot){
        single_vessel_plotter(vessel_coords = smth_vessel, new = TRUE, col = "purple")
      }
    }
  }
   
  # Here we initialize the Frenet-Serrat frame vectors.
  tic("frenet")
  tangent_array <- mat.or.vec(length(smth_vessel[,1]), 3)
  normal_array <- mat.or.vec(length(smth_vessel[,1]), 3)
  binormal_array <- mat.or.vec(length(smth_vessel[,1]), 3)

  tangent_array[,] <- NaN
  normal_array[,] <- NaN
  binormal_array[,] <- NaN

  for(i in 3:(length(smth_vessel[,1]) - 2)){
    tangent_array[i,] <- tangent_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    normal_array[i,] <- normal_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
    binormal_array[i,] <- binormal_vector(smth_vessel[i-2,], smth_vessel[i-1,], smth_vessel[i,], smth_vessel[i+1,], smth_vessel[i+2,])
  }
  toc()
  toc()

  return(list(tangent = tangent_array, normal = normal_array, binormal = binormal_array, vssl_coords = smth_vessel))
}
