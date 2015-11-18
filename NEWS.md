rSPACE 1.0.0
----------------------------------  

* Yay, first release!
* Renaming variables to hopefully be more intuitive (moveDist, moveDistQ, maxDistQ, habitat.cutoff)
* showSteps option in `encounter.history()`

##### rSPACE 1.0.4 (5 Jan 2015)
###### Bug Fixes
* Parameter mismatch in dropN fixed. Had been causing errors for declining populations
* NAs in habitat map were creating errors in encounter.history().  Added repeat map checking steps in encounter.history() when showSteps=T.
* Added warning message for outdated parameter names

##### rSPACE 1.0.5 (15 Jan 2015)
###### Changes
* No longer allows NAs in Habitat or filter maps.  Returns errors rather than trying to fix
* Add saveGrid option to create.landscapes() with default to TRUE.  Will produce a raster layer of the grid used in the simulation in the output folder.

###### Bug fixes
* Error message added for empty grids
* Fixed ratio for creating effective sampling areas
* Error message added for trying to set effective sampling area larger than grid size

##### rSPACE 1.1.0 (27 Feb 2015)
###### Changes
* Cleaning up and adding example scripts
* Updating documentation

###### Bug fixes
* Cleaning up clicky box errors, and setting default detP=1
* Bugs for printing figures in showSteps=T in encounter.history()
* Saving grid from encounter.history() with showSteps=T

##### rSPACE 1.2.0 (17 Nov 2015)
###### Changes
* Deprecating create.landscapes() and test_samples()
* Allow passing cell IDs to analysis function
* Updating check procedure for previous Parameters file.
* Now produces error if trying to overwrite previous output files.
* Adds 'add' argument to createReplicates() to prevent unintending overwriting of encounter history files Defaults to add=F, which returns an error if previous encounter history files are present and restarts N_final.txt counts

###### Bug fixes
* Check for NaNs when dropping individuals in encounter.history()
* Updates labels for RMark results in wolverine_analysis
* Match extent for SamplingFrame and Habitat example layers
