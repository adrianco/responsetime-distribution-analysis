# guillotine.R - chop logfile response time data into sections and process into peaks

# ChatGPT query resulting in the initial version of the code below
# I now have a list of peaks dataframes,  each contains about ten density/latency pairs of data.
# What I want is an iterative clustering algorithm for a stream of data frames over time.
# Each data frame contains a list of points.
# Each point to be clustered is a pair of (x,y) coordinates.
# The first data frame sets the starting points for the clusters.
# Each new data frame that is added is matched to the cluster and points that don't match closely enough create new clusters.
# The current cluster state is returned.
# There is a separate function that adds a new data frame to be clustered. 

addPeaks <- function(peaks, clusters, normalize, threshold=0.05) {
  for (i in 1:nrow(peaks)) {
    point <- peaks[i, c(8, 1)]  # Extract the (x, y) coordinates Latency and Density
    if (point[1] <= 0) next # sometimes latency rounds down to zero, skip to avoid log(0)
    
    point[1] <- log(point[1])/normalize # normalize log latency to consistent max bucket
    
    # Calculate distances to existing clusters
    distances <- sapply(clusters, function(cluster) sqrt((point[1] - cluster$centroid[1])^2 + (point[2] - cluster$centroid[2])^2))
    
    closest_cluster_index <- which.min(distances)
    closest_cluster <- clusters[[closest_cluster_index]]
    
    if (distances[closest_cluster_index] <= threshold) {
      closest_cluster$peaks <- rbind(closest_cluster$peaks, peaks[i, ])
      closest_cluster$points <- rbind(closest_cluster$points, point)
      closest_cluster$centroid <- colMeans(closest_cluster$points)
      clusters[[closest_cluster_index]] <- closest_cluster
    } else {
      new_cluster <- list(peaks = peaks[i, ], points = point, centroid = colMeans(point), normalize=normalize) # has to match the clusters, not a df
      clusters <- c(clusters, list(new_cluster))
    }
  }
  
  return(clusters)
}

# process a log file by chopping into one minute chunks, finding peaks, then clustering them
# return a list of clusters and the corresponding peaks
# initial code structure by ChatGPT
guillotine <- function(df, plot=F, epsilon=0.01, peakcount=10) {
  start_time <- round(min(df$time), "mins")
  end_time <- round(max(df$time), "mins")
  minute_intervals <- seq.POSIXt(start_time, end_time, by = "min")

  #print(minute_intervals)
  # Split the data frame into subsets based on the minute intervals
  msc <- cut(df$time, breaks = minute_intervals, right = FALSE, labels = FALSE)
  
  # Find the last non-NA interval
  last_interval <- max(msc, na.rm = TRUE)
  
  # Replace NAs with an additional interval that captures data before and after the minutes
  # maybe exclude data rather than add an extra interval?
  msc[is.na(msc)] <- last_interval+1
  df$msc <- msc
  
  # use consistent breaks for all the peaks - get breaks for the whole dataset
  hb <- hist(log(df$latency), breaks=40, plot=plot)$breaks
  mhb <- max(hb) # max histogram bucket - needed to normalize latency
  
  # Process each minute section and store results in a list
  results_list <- lapply(unique(df$msc), function(section) {
    subset_df <- df[df$msc == section, , drop = FALSE]
    
    # Find the peaks in each subset
    peaks <- as.peaks(hist(log(subset_df$latency), breaks=hb, plot=F), time=subset_df$time[1], normalize=T, epsilon=epsilon, peakcount=peakcount, plots=plot)
    
  })
  
  # Initialize the clusters using the peaks from the first minute
  first_data_frame <- results_list[[1]]
  first_points <- first_data_frame[, c(8, 1)]  # Extract the (x, y) coordinates Latency and Density
  first_points[,1] <- log(first_points[,1])/mhb   # normalize log latency to max bucket
  
  initial_clusters <- lapply(1:nrow(first_points), function(i) {
    list(peaks = first_data_frame[i, ], points = first_points[i, ], centroid = as.matrix(first_points[i, ]), normalize = mhb)
  })
  
  # Process subsequent data frames
  for (i in 2:length(results_list)) {
    new_data_frame <- results_list[[i]]
    initial_clusters <- addPeaks(new_data_frame, initial_clusters, mhb)
  }
  
  # Return the current cluster state
  initial_clusters
}

library(ggplot2)
# asked ChatGPT "function to plot the current state of the clusters" got this
plotClusters <- function(clusters, xlab="Normalized Log-Latency", ylab="Density") {
  points <- data.frame(x = numeric(), y = numeric(), cluster = factor())
  
  for (i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    cluster_points <- cluster$points
    points <- rbind(points, data.frame(x = cluster_points[, 1], y = cluster_points[, 2], cluster = as.factor(i)))
  }
  
  # Plot the points with different colors for each cluster
  p <- ggplot(points, aes(x = x, y = y, color = cluster)) +
    geom_point() +
    theme_minimal() +
    labs(x = xlab, y = ylab)
  
  # Add black points for the centroids
  centroids <- lapply(clusters, function(cluster) cluster$centroid)
  centroids_df <- data.frame(x = sapply(centroids, "[", 1), y = sapply(centroids, "[", 2))
  #print(centroids_df)
  p <- p + geom_point(data = centroids_df, aes(x = x, y = y), color = "black")
  
  # Return the plot
  return(p)
}

plotClusterDensity <- function(clusters) {
  # Create a data frame to store the cluster density and time information
  cluster_data <- data.frame(Time = numeric(), Cluster = numeric(), Density = numeric())
  
  # Extract cluster density and time information from each cluster
  for (i in 1:length(clusters)) {
    cluster <- clusters[[i]]
    density <- cluster$points$PeakDensity
    time <- cluster$peaks$Time
    cluster_data <- rbind(cluster_data, data.frame(Time = time, Cluster = i, Density = density))
  }
  
  # Plot the cluster densities over time
  p <- ggplot(cluster_data, aes(x = Time, y = Density, group = Cluster, color = as.factor(Cluster))) +
    geom_line() +
    geom_point() +
    labs(x = "Time", y = "Peak Density", color = "Cluster") +
    theme_minimal()
  
  return(p)
}

# haven't quite got this plot figured out yet...
plotClusterPercentile <- function(clusters) {
  # Create an empty data frame to store the combined distributions
  combined_dist <- data.frame(latency = numeric(), density = numeric(), cluster = integer())
  
  # Iterate over each cluster
  for (i in seq_along(clusters)) {
    cluster <- clusters[[i]]
    
    # Calculate the combined standard deviation using the geometric mean
    combined_sd <- exp(mean(log(cluster$peaks$PeakSD)))
    
    # Calculate the combined distribution using dnorm
    latency <- seq(0, cluster$normalize, length.out = 1000)  # Adjust the range and resolution as needed
    density_combined <- dnorm(latency, mean = cluster$centroid[1] * cluster$normalize, sd = combined_sd)
    
    # Add the combined distribution to the data frame
    combined_dist <- rbind(combined_dist, data.frame(latency = latency, density = density_combined, cluster = i))
  }
  
  # Plot the combined distributions
  ggplot(combined_dist, aes(x = latency, y = density, color = factor(cluster))) +
    geom_line() +
    xlab("Latency") +
    ylab("Density") +
    ggtitle("Combined Distributions") +
    theme_minimal()
}
