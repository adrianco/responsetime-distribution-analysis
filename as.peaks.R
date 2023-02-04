# Take a histogram h$Value h$Count and find the peaks in it, interpolating between bins
# approximate finite mixed model decomposition starting with histograms of data
# Intended for use on response time histograms with logarithmic bins like hdrhistogram
# so that each (approximately) lognormal component is interpreted as a symmetric normal peak
library(pracma) # for findpeaks

# h - is a data frame containing Values and Counts forming a histogram where Values for each bucket are exponential
# plots - optionally shown
# normalize - divides down the counts to a probability density
# peakcount=0 does a quick estimate of all the peaks, otherwise it does up to the specified number of dnorm fits
# epsilon=0.005 - sets a minimum interesting peak size to 0.5% of the largest peak, the ones that are visible in the plot.
# printdebug turns on helpful intermediate print output to see what's going on

# Returns a data frame with one row per peak, sorted biggest peak first
# row names show the order the peaks were discovered in
# PeakDensity is the histogram count or normalized density
# PeakBucket is the center bucket for the peak
# PeakMin and PeakMax are the start and end buckets for the peak
# if a gaussian fit was obtained, then PeakAmplitude is non-zero, along with updated PeakMean and PeakSD
# The Value of the peak is interpolated to estimate PeakLatency
as.peaks <- function(h, plots=FALSE, normalize=FALSE, epsilon=0.005, peakcount=0, printdebug=F) {
  # normalize counts to a probability density
  if (normalize) {
    yl <- "Probability Density"
    h$Count <- h$Count/sum(h$Count)
  } else {
    yl <- "Count"
  }
  # mean of the entire distribution
  hmean <- sum(h$Count*h$Value)/sum(h$Count)
  # add bucket as a column for fitting and use this local copy to subtract peaks from Count
  hb <- cbind(h, Bucket=as.numeric(row.names(h)))
  if (plots) {
    plot(h$Value, h$Count, type="l",main="Response Time Distribution and Mean",
         xlab=paste("Mean X=", format(hmean, digits=4), "Seconds"), ylab=yl)
    points(hmean, 0, col=2, pch='X') # mark the mean value
    barplot(h$Count, main="Response Time Histogram with Logarithmic Buckets", ylab=yl)
    plot(h$Count, type="l", main="Response log(Time) with Peaks", xlab="Bucket Index", ylab=yl)
  }
  # returns the position of the peaks and the tails on each side, sorted biggest first
  p <- findpeaks(h$Count, sortstr=T, minpeakheight = epsilon*max(h$Count))
  # assemble into data frame with extra columns for gaussian peaks and interpolated peak latency
  peaks <- data.frame(PeakDensity=p[,1], PeakBucket=p[,2], PeakMin=p[,3], PeakMax=p[,4],
             PeakMean=0, PeakSD=0, PeakAmplitude=0, PeakLatency=0)
  # do a quick match of all peaks in one shot if peakcount=0, otherwise subtract out and re-find peaks
  if (peakcount >0)
    iters <- peakcount
  else
    iters <- nrow(peaks)
  # iterate over the peaks fitting them more accurately and subtracting them out
  newpeaks <- FALSE # first pass, only process peaks seen before subtraction to reveal more peaks
  for (i in 1:iters) {
    if (i > nrow(peaks)) break # we ran out of peaks before we got to peakcount
    p.bucket <- peaks$PeakBucket[i] # the bucket index at the center of the peak
    p.min <- max(peaks$PeakMin[i], p.bucket-2) # take at most center five points to fit peak to
    p.max <- min(peaks$PeakMax[i], p.bucket+2)
    # mean of the peak in between buckets
    peaks$PeakMean[i] <- sum(p.min:p.max * hb$Count[p.min:p.max])/sum(hb$Count[p.min:p.max])
    # estimate sd based on the peak only, which will be below the actual sd but in the right range
    peaks$PeakSD[i] <- sqrt(sum(hb$Count[p.min:p.max]*((p.min:p.max - peaks$PeakMean[i])^2))/(sum(hb$Count[p.min:p.max])))
    # fit the top few peaks with guassian/normal curves and subtract out the peaks as we go
    # first do all the peaks that were seen at the beginning before looking for new ones
    if (peakcount > 0) {
      if (printdebug) print(peaks[i,])
      if (printdebug) print(hb[p.min:p.max,])
      # fit a normal to it, starting with an underestimate to get the narrowest fit without crashing nls
      # need to use try as it craps out especially on small residual peaks
      fit <- try(nls(Count ~ a * dnorm(Bucket, mean, sd), data=hb[p.min:p.max,],
                 start=list(mean=peaks$PeakMean[i], sd=peaks$PeakSD[i]/2, a=peaks$PeakDensity[i]),
                 control=nls.control(warnOnly=F, printEval=F)), silent=T)
      if (class(fit) == "try-error") {
        if (printdebug) print(fit[1])
      } else {
        cfit <- coef(fit)
        if (printdebug) print(cfit)
        peaks$PeakAmplitude[i] <- cfit[3] # normal amplitude is not the same as bucket density
        peaks$PeakMean[i] <- cfit[1] # overwrite previous estimates
        peaks$PeakSD[i] <- cfit[2]
        p.norm <- predict(fit, newdata=list(Bucket=1:nrow(hb)))
        # Subtract out the current peak from the data
        hb$Count <- hb$Count - p.norm
        hb$Count[hb$Count < 0] <- 0 # truncate any negative residuals
        if (plots) {
          lines(hb$Bucket, p.norm, col=i) # the fitted normal
          lines(hb$Bucket, hb$Count, col=i+1) # the residual that will be fitted next time in the same color
          points(peaks$PeakMean[i], peaks$PeakDensity[i], col=i, pch="x")
        }
        # check to see if we have run out of peaks
        if (i==nrow(peaks)) {
          if (newpeaks==TRUE) {
            break # nothing more to do
          } else {
            # process the remains of the distribution after initial peaks
            newpeaks <- TRUE
            # rescan once for new side peaks revealed by removal
            # that are at least as big as the smallest we've already seen, and at least 5 buckets wide
            p <- findpeaks(hb$Count, minpeakheight=min(peaks$PeakDensity), sortstr=T) # ndowns = 2, nups = 2)
            if (printdebug) {
              print(p)
            }
            if (!is.null(p) && nrow(p) > 0) {
              # CHANGE TO MERGE - need to remove any peaks that are at same buckets as previous run
              # append new peaks to end
              peaks[(i+1):(i+nrow(p)),] <- cbind(p, array(0,c(nrow(p),4))) # pad p to fit the data frame
            } else {
              break # stop trying to process any more peaks
            }
          }
        }
      }
    }
    # interpolate between the values of the buckets each side of the peak to get the latency
    p.logbucketvalue <- log(h$Value[p.bucket]) # the time value at that bucket
    peaks$PeakLatency[i] <- exp((log(hb$Value[p.bucket+1]) - p.logbucketvalue) *   
                                  (peaks$PeakMean[i] - p.bucket) + p.logbucketvalue)
  }
  # final sorted return value
  if (plots && peakcount == 0) points(peaks$PeakMean, peaks$PeakDensity, col=1, pch="x")
  if (iters < nrow(peaks)) peaks <- peaks[1:iters,] # trim any un-used peaks
  peaks[order(peaks$PeakDensity, decreasing=T),]
}