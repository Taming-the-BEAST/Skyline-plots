#library(bdskytools)

# Requires packages RColorBrewer, boa

fname <- "~/Documents/Repositories/TamingTheBeast/Tutorials/Skyline/Data/hcv_bdsky.log"

lf <- readLogfile(fname, burnin=0.1)

R0_sky    <- getSkylineSubset(lf, "R0")

# Extract the raw HPDs
R0_hpd    <- getMatrixHPD(R0_sky)
delta_hpd <- getHPD(lf$becomeUninfectiousRate)


# Extract gridded HPDs  (if the timegrid is too fine it will not be that smooth anymore, but it's fast enough to quickly interpolate a grid of 500)
# timegrid <- 1:500
timegrid <- seq(1,400,length.out=100)
R0_gridded    <- gridSkyline(R0_sky,    lf$origin, timegrid)

R0_gridded_hpd    <- getMatrixHPD(R0_gridded)


# The plotting times, the most recent sample is 1993
times     <- 1993-timegrid


# Plain skyline plots
plotSkyline(times, R0_gridded_hpd, type='smooth')
plotSkyline(times, R0_gridded_hpd, type='lines')

# Plot the traces of the skyline, so you can see how the average gives a smooth hpd
plotSkyline(times, R0_gridded, type='steplines', traces=10, col=pal.dark(cblue,0.5),ylims=c(0,5))
plotSkyline(times, R0_gridded, type='steplines', traces=100, col=pal.dark(cblue,0.5),ylims=c(0,5))
plotSkyline(times, R0_gridded, type='steplines', traces=1000, col=pal.dark(cblue,0.1),ylims=c(0,5))

# The non-gridded intervals
# Since here they are equidistant from the origin to the present they should probably not be plotted this way
# Use this when the shift-times in bdsky are fixed (needs to be manually set in the xml at the moment)
# or when the origin is fixed (automatically fixes the shift-times)
# When doing this do NOT plot with type='smooth' as it gives a misleading result!!!
plotSkyline(1:10, R0_hpd, type='step')
plotSkyline(1:10, R0_hpd, type='smooth')
plotSkyline(range(times), as.matrix(delta_hpd), type='step')



# Pretty skyline
plotSkylinePretty(times, R0_gridded_hpd, axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[0]), side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3))


# Plot both R0 and delta skylines on one set of axes
# Can also use this to compare skylines of the same parameter between different models (eg. changing the priors or number of shifts)
par(mar=c(5,4,4,4)+0.1)
plotSkylinePretty(range(times), as.matrix(delta_hpd), type='step', axispadding=0.0, col=pal.dark(cblue), fill=pal.dark(cblue, 0.5), col.axis=pal.dark(cblue),
ylab=expression(delta), side=4, yline=2, ylims=c(0,1), xaxis=FALSE)
plotSkylinePretty(times, R0_gridded_hpd, type='smooth', axispadding=0.0, col=pal.dark(corange), fill=pal.dark(corange, 0.5), col.axis=pal.dark(corange),
xlab="Time", ylab=expression("R"[0]), side=2, yline=2.5, xline=2, xgrid=TRUE, ygrid=TRUE, gridcol=pal.dark(cgray), ylims=c(0,3), new=TRUE, add=TRUE)
