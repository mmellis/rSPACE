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

