### This file is meant to serve as an example workflow for analyzing tortuosity of individual vessels.

##  First set your working directory to a location consisting of the .dat files and your directory of R functions

setwd("/Users/alexwork/odrive/Google Drive bruinmail/Research/ucla/savage_lab_vascular_build/tortuosity_calcs_R/code/final/")

source("./dat_file_reader.R")
source("./single_vessel_plotter.R")
source("./vessel_poly_fit.R")
source("./vessel_spline_fit.R")
source("./frenet_vectors.R")
source("./tortuosity_metrics.R")
source("./curve_torse_check_plotter.R")

library("rgl")

## Use the dat_file_reader() function to read in a .dat file for analysis

filename <- "MCA0 day 7 fitc N-13 NIH001_sm_vd_824x824x1000_th_05.dat"

vessels_slice <- dat_file_reader(dat_filename = filename)

## Run the View() function to see the structure of the vessels_slice data.frame.  Note that it contains rescaled coordinates for the vessel backbones (in units of micrometers) and the fourth column is a factor for vessel ID numbers.  Note that these IDs are NOT the same as the vessel nodeids from your other .tsv files.  Integration of the correct vessel ID names happens after performing tortuosity analyses.
View(vessels_slice)

## Quickly use the RGL packages plot3d function to get a coarse view of all vessels in the slice.
plot3d(vessels_slice[,1:3])

## Alternatively, use the custom made single_vessel_plotter() function for viewing individual vessels.  Index via the which function to capture individual vessels based on their ID number.

# Extracting an individual vessel.  Note that we have also indexed so as to remove the ID column from our single vessel.

vessel <- vessels_slice[which(vessels_slice$ID == 29), 1:3]

single_vessel_plotter(vessel_coords = vessel, centered = TRUE, frenet = FALSE)

## Run one of the fitting functions, either vessel_poly_fit or vessel_spline_fit, to interpolate a fitted functional representation of the vessel.  The bulk of this project will be varying the fitting parameters (smoothing for vessel_spline_fit and poly_degree for vessel_poly_fit) to examine how they influence the resulting interpolation and subsequent measures for tortuosity.  Outputs from these functions will be the x, y, and z compoents of the tangent, normal and binormal Frenet vectors, as well as reprinting the x, y, and z coordinates of the voxel backbones.  Also note that the output from these functions are lists of matrices.

vessel_poly_fit(vessel = vessel, number_samples = 10000, poly_degree = 3, plot = TRUE)
vessel_poly_fit(vessel = vessel, number_samples = 10000, poly_degree = 5, plot = TRUE)
vessel_poly_fit(vessel = vessel, number_samples = 10000, poly_degree = 7, plot = TRUE)

# Note that the null option for splining uses a leave-one-out cross-validation method to automatically calculate the ideal pentalizing parameter.  However, this is often an overfitting.
vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = NULL, plot = TRUE)
vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = 0.3, plot = TRUE)
vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = 0.5, plot = TRUE)
vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = 0.7, plot = TRUE)


smth_vessel <- vessel_spline_fit(vessel = vessel, number_samples = 10000, smoothing = 0.5, plot = TRUE)


## Now that we have smoothed interpolated vessels, we can add the Frenet vectors to the visualizaitons if we like when using single_vessel_plotter().
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = TRUE, scale = 1.0, frenet_vectors = seq(from = 100, to = 5000, by = 100))

## Next we'll make graphs of curvature, torsion, and error in curvature and torsion as functions of normalized vessel arc length, using the curve_torse_check_plotter() function.

# To graph curvature and torsion, use plot_type = 1
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 1)

# To graph curvature and error in curvature, use plot_type = 2
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 2)

# To graph torsion and error in torsion, use plot_type = 3
curve_torse_check_plotter(vessel_coords = smth_vessel, plot_type = 3)

## To extract various measures of tortuosity, we can use a combination of the functions curvature_torsion_calculator(), distance_metric(), inflection_count_metric, and sum_of_all_angles_metric().  The former of these four methods measures produces multiple metrics all directly related to curvature and torsion.  These are: total and average and maximum curvature and torsion, as well as total and average combined curvature and torsion.  These metrics are discussed in O'Flynn et al., 2007.  The latter three mtrics of distance, inflection count, and sum of all angles are discussed in Bullitt et al., 2003.

cur_tor_met <- curvature_torsion_calculator(tangent = smth_vessel[[1]], normal = smth_vessel[[2]], binormal = smth_vessel[[3]], vessel_coords = smth_vessel[[4]])

## Distance metric returns both the arclength and the arclength divided by the end-to-end distance
dist_met <- distance_metric(vessel_coords = smth_vessel[[4]])

## Inflection count metric can be toggled to return either the number of inflection points, or the values of the (Delta N)**2 term as a function of vessel arclenth which indicate inflection points when (Delta N)**2 > 2.
inflec_met <- inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = TRUE, return_coords = FALSE)
inflec_met <- inflection_count_metric(vessel_coords = smth_vessel[[4]], return_count = FALSE, return_coords = TRUE)

## Sum of all angles metric returns the sum of all angles as defined in Bullitt et al., 2003.

soam <- sum_of_all_angles_metric(smth_vessel[[4]])

## Finally, we can plot some tortuosity metrics on the vessels them selves using color coding.  Here are two examples using curvature and torsion to color the vessel.


# Coloring by curvature.
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = FALSE, col_metric = cur_tor_met$curvature)


# Coloring by torsion
single_vessel_plotter(vessel_coords = smth_vessel[[4]], centered = TRUE, frenet = FALSE, col_metric = cur_tor_met$torsion)



## Use the above code to experiment with different different fitting parameters for different vessels (and from different slices).  More detailed instruction will come after the weekend.
