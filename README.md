rSPACE
====

Power analysis for detecting trends in population abundance using occupancy-based monitoring.  A simple, spatially explicit population simulation is used to create sample datasets with which to assess power.

We are currently working to test and develop the code to be more generally accessible.  If you have questions concerning how to implement rSPACE or if there are features that you would like to see added (or if you have any other feedback on the current implementation), please contact us.  Thanks for your patience as we work through the bugs!





##### Download steps
---

1. Make sure you have [Program MARK](http://www.phidot.org/software/mark/downloads/) and the most current version of [R](http://cran.us.r-project.org/) (or at least 3.0.1)
2. Install package dependencies in R

    ````R
    install.packages(c("raster", "RMark", "ggplot2", "tcltk2"))
    ````

3. For 64 bit Windows - download the zip binary from github to someplace handy.
  - This may be a little trickier than it seems.  If you click the handy "Download ZIP" button on the right there, you get the whole github folder.  Within that, the zip binary is the file labeled rSPACE_1.0.zip
4. Install rSPACE package in R
  
    ````R
    install.packages("C:/directory where the zipfile is/rSPACE_1.0.zip", repos=NULL, type='source')
    ````

5. load rSPACE library in R 

    ````R 
    library(rSPACE)
    ````

6. Find the handy package overview

    ````R
    help(rSPACE)
    ````
