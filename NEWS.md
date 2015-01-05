rSPACE 1.0.0
----------------------------------  

* Yay, first release!
* Renaming variables to hopefully be more intuitive (moveDist, moveDistQ, maxDistQ, habitat.cutoff)
* showSteps option in `encounter.history()`

rSPACE 1.0.4 Bug Fixes
* Parameter mismatch in dropN fixed. Had been causing errors for declining populations
* NAs in habitat map were creating errors in encounter.history().  Added repeat map checking steps in encounter.history() when showSteps=T.
* Added warning message for outdated parameter names
