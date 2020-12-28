# NetworkExtinction 0.1.1

* Added a `NEWS.md` file to track changes to the package.

# NetworkExtinction 0.1.2

* Fixed parameters in `degree_distribution`, the intercept for power law was not present
* Fixed `degree_distribution` to select one model per family
* Eliminated *Truncated* distribution from `degree_distribution` since it didn't have any theoretical support
* Eliminated the name *argument* `degree_distribution`

# NetworkExtinction 0.1.3

* Added parallel processing for `RandomExtinctions`

# NetworkExtinction 0.2.1

* Change function name from `degree_distribution` to `DegreeDistribution`
* The functions `Mostconnected` and `ExtinctionOrder` are now soft deprecated in favour of `SimulateExtinctions`
* Transformed `ExtinctionOrder` output to unify with `Mostconnected` output so that `SimulateExtinctions` has a common output
