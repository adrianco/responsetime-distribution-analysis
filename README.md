# Responsetime Distribution Analysis
R code to analyze response time distributions as mixed models, and to cluster seqences of peaks over time as discussed in my Monitorama 2023 talk.

As initially described by [this medium post](https://medium.com/@adrianco/percentiles-dont-work-analyzing-the-distribution-of-response-times-for-web-services-ace36a6a2a19) and extended for Monitorama June 2023 [slides here](https://github.com/adrianco/slides/blob/master/Monitorama%20Histograms.pdf)


Some basic sample data provided as a csv. Will also recognize R histogram class as input directly.
```
> sample <- read.csv("sample.csv")
> as.peaks(sample,plots=T,peakcount=6, normalize = T)
  PeakDensity PeakBucket PeakMin PeakMax  PeakMean    PeakSD PeakAmplitude PeakLatency
4 0.008936170          4       3       5  3.878831 0.3263232    0.00000000    39.47890
3 0.011489362          6       5       8  6.246079 1.2272864    0.03091961    63.38428
5 0.007659574          9       8      10  8.821677 0.3827842    0.00000000   106.09506
2 0.173191489         25      20      26 24.483727 0.7152434    0.00000000  2432.67183
1 0.222553191         28      26      31 27.322500 1.0325453    0.64311539  4291.96511
```

![response](sample-response.png)
![barplot](sample-barplot.png)
![peaks](sample-peaks.png)
