# Install the package with devtools::install_github(laduplessis/bdskytools)
# If you cannot install the package source the R files in the folder and make sure the 
# packages boa and RColorBrewer are installed

# Load the package
library(bdskytools)


# Extract the data and HPDs from the logfile
############################################

# Change this to the path of your logfile here!
fname <- "~/Documents/Taming_the_BEAST/Tutorials-Git/Skyline-plots/precooked_runs/hcv_bdsky.log"
lf    <- readLogfile(fname, burnin=0.1)

Re_sky    <- getSkylineSubset(lf, "reproductiveNumber")
Re_hpd    <- getMatrixHPD(Re_sky)
delta_hpd <- getHPD(lf$becomeUninfectiousRate)


# The non-gridded intervals
###########################
# Since here they are equidistant from the origin to the present they should probably not be plotted this way
# Use this when the shift-times in bdsky are fixed (needs to be manually set in the xml at the moment)
# or when the origin is fixed (automatically fixes the shift-times)
# When doing this do NOT plot with type='smooth' as it gives a misleading result!!!
plotSkyline(1:10, Re_hpd, type='step', ylab="R")

# Why is the result below misleading?
plotSkyline(1:10, Re_hpd, type='smooth', ylab="R")


# Extract gridded HPDs  
######################
# The trimegrid specifies the years into the past to grid the skyline to.
# We grid the skyline to every 4th year for the past 400 years
# If the timegrid is too fine it will not be that smooth anymore, 
# but it's fast enough to quickly interpolate a grid of 400 (every year for the past 400 years)
timegrid <- seq(1,400,length.out=100)
#timegrid <- 1:400
Re_gridded     <- gridSkyline(Re_sky,    lf$origin, timegrid)
Re_gridded_hpd <- getMatrixHPD(Re_gridded)


# Plotting the skyline
#####################
# The plotting times, the most recent sample is from 1993
times     <- 1993-timegrid
plotSkyline(times, Re_gridded_hpd, type='smooth', xlab="Time", ylab="R")
plotSkyline(times, Re_gridded_hpd, type='lines', xlab="Time", ylab="R")

Re_rev <- revMatrix(Re_sky)
Re_rev_gridded     <- gridSkyline(Re_rev, lf$origin, timegrid)
Re_rev_gridded_hpd <- getMatrixHPD(Re_rev_gridded)
plotSkyline(times, Re_rev_gridded_hpd, type='smooth', xlab="Time", ylab="R")

Re_rev_gridded2     <- gridSkyline(Re_sky, lf$origin, timegrid, reverse=TRUE)
Re_rev_gridded_hpd2 <- getMatrixHPD(Re_rev_gridded2)
plotSkyline(times, Re_rev_gridded_hpd2, type='smooth', xlab="Time", ylab="R")


# Where do the smooth skyline come from?
########################################
# Plot the traces of the skyline, so you can see how the average gives a smooth hpd
# First 10, then 100, then 1000 of the individual MCMC steps. Now if we take the HPD 
# interval at every fourth year and only plot that we get the smooth skyline in the
# previous step.
plotSkyline(times, Re_gridded, type='steplines', traces=10, col=pal.dark(cblue,0.5),ylims=c(0,5), xlab="Time", ylab="R", main="10 traces")
plotSkyline(times, Re_gridded, type='steplines', traces=100, col=pal.dark(cblue,0.5),ylims=c(0,5), xlab="Time", ylab="R", main="100 traces")
plotSkyline(times, Re_gridded, type='steplines', traces=1000, col=pal.dark(cblue,0.1),ylims=c(0,5), xlab="Time", ylab="R", main="1000 traces")


# Pretty skyline plots
######################
# Use this function for the type of plot you might want to publish
plotSkylinePretty(times, Re_gridded_hpd, axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3))

# Add a line at Re = 1 for comparison
abline(h=1, col=pal.dark(cred), lty=2)


# Plot both Re and delta skylines on one set of axes
# Delta only has a dimension of 1, so the skyline is not really insightful
# Can also use this to compare skylines of the same parameter between different models (eg. changing the priors or number of shifts)
plotSkyline(c(0,1), as.matrix(delta_hpd), type='step', ylab="Delta")
par(mar=c(5,4,4,4)+0.1)
plotSkylinePretty(range(times), as.matrix(delta_hpd), type='step', axispadding=0.0, col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), col.axis=pal.dark(cblue),
ylab=expression(delta), side=4, yline=2, ylims=c(0,1), xaxis=FALSE)
plotSkylinePretty(times, Re_gridded_hpd, type='smooth', axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[e]), side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3), new=TRUE, add=TRUE)
