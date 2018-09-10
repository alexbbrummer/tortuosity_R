##### This worksheet presents example code for scaling up tortuosity analysis to be automated on a whole slice of vessels.  Various comments on current bugs/assumptions are provided.


## Initialize list of all vessels in slice and extract each vessel from slice based on Unique ID

vessel_list <- list()
for(i in 1:length(levels(as.factor(vessels_slice$ID)))){
  vessel_list[[i]] <- vessels_slice[which(vessels_slice$ID == i), 1:3]
}

## Initlize smoothed vessel list and perform smoothing on each vessel that containts 20 or more points.  This number is dependent on the number of knots being used for the fitting as well as the degree of the splines.  Since we are varying knot number with different models, the 20 point threshold could be defined dynamically in future versions of this code, but currently it is being hardcoded with the simple "if" statement.
smth_vessel_list <- list()
for(i in 1:length(vessel_list)){
  if(nrow(vessel_list[[i]]) >= 20){
    # print(i)
    smth_vessel_list[[i]] <- vessel_spline_fit(vessel = vessel_list[[i]], number_samples = 1000, spline = "pspline", m = c(4,2), plot = FALSE, subsample_density = 2)
  }else{
    smth_vessel_list[[i]] <- NULL
  }
}

# Initialize tortuosity metric lists
cur_tor_met_list <- list()
dist_met_list <- list()
inflec_met_list <- list()
soam_list <- list()

# Loop through smoothed and interpolated vessels and calculate various tortuosity metrics.  Note that here we are keeping the filtering of torsion spikes turned off.  As coded, the filteirng process is too buggy for scaled up automation.  In future versions of this code and work could improve upon automation either with better error handling, or using analytic methods for prediction and evaluation of spike locations (inflection points).  See comments in the filter_torsion_spikes method within tortuosity_metrics.
for(i in 1:length(smth_vessel_list)){
  if(length(smth_vessel_list[[i]]) == 4){
    print(i)
    cur_tor_met_list[[i]] <- curvature_torsion_calculator(tangent = smth_vessel_list[[i]][[1]], normal = smth_vessel_list[[i]][[2]], binormal = smth_vessel_list[[i]][[3]], vessel_coords = smth_vessel_list[[i]][[4]], filter_torsion_spikes = FALSE)
    dist_met_list[[i]] <- distance_metric(vessel_coords = smth_vessel_list[[i]][[4]])
    inflec_met_list[[i]] <- inflection_count_metric(vessel_coords = smth_vessel_list[[i]][[4]], return_count = TRUE, return_coords = FALSE)
    soam_list[[i]] <- sum_of_all_angles_metric(smth_vessel_list[[i]][[4]])
  }
}
