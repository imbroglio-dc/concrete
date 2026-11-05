# concrete 1.1.0

## Bug fixes

* The one-step TMLE hazard update now uses the intervention-evaluated
  g*(a* | W) (`g.star.intervention`) for the update direction. Previously it
  used `g.star.obs`, which left the counterfactual hazards of subjects whose
  observed treatment differed from the target intervention un-updated, so the
  estimator was not doubly robust. The efficient influence curve still uses
  `g.star.obs`. **Point estimates and standard errors will change** relative
  to 1.0.8, particularly when the outcome (hazard) models are misspecified.
* `checkTreatment()` no longer applies `is.nan()` / `is.infinite()` to
  non-numeric treatment columns.
* `getID()` guards against non-character `ID` arguments.

## Improvements

* `getPropScore()` is more robust: it validates that intervention values are
  binary and that each regime has a `g.star` function, and it gives clearer
  errors when a `g.star` function fails.
* Invalid `MaxUpdateIter`, `OneStepEps`, `MinNuisance`, and unrecognized
  `Model` specifications now raise warnings (previously messages) when falling
  back to defaults.
* Logical arguments accept common string and numeric encodings of TRUE/FALSE.
* `CVArg$V` larger than the effective sample size is reset to the default
  folds.
* An informative message is given when a subset of `TargetEvent` is supplied
  (not yet supported).

## Other changes

* License changed from GPL-3 to MIT.
* CI: updated GitHub Actions versions; coverage reporting moved to Codecov.
* Expanded `formatArguments()` test suite.
