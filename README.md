
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `mosmicrosim`: mosquito microclimate simulator

<!-- badges: start -->
<!-- badges: end -->

A prototype R package for simulating a timeseries of microclimatic
conditions - e.g. local water surface area, water and air temperature
and humidity - in a given location on earth, for a mosquito species,
given a simple set of assumptions about the species’ habitat
preferences.

Most of the computational work is done by the [NicheMapR R
package](https://mrke.github.io/) and its internal Fortran routines for
physical simulation of microclimate conditions. This package provides a
simplified interface for specifying mosquito habitat types, linking
NicheMapR to custom spatial climate/weather data, and including some
additional simulation models for different types of small water bodies.

So far the package provides a single (still in development) function:
`create_micro()` which provides a more user-friendly interface to
creating the config function for the `NicheMapR::microclimate()`
function to compute microclimatic conditions at a given location with
user-specified weather timeseries, for a simplified microclimate model
corresponding to mosquito habitats.

The to do list for the rest of the package is as follows:

- [x] Import NicheMapR and set up package
- [x] Write a shim for functions like `NicheMapR::microclimate()` to
  enable incorporation of other sources of climate timeseries data
- [ ] Streamline simulation of microclimate components *not* including
  local water surface area, and water and air temperature, and humidity
- [ ] Implement simulation of local water surface area, and water and
  air temperature, and humidity using a set of water body models (cone,
  cylinder, with/without leakage)
- [ ] Write a simplified interface to simulating microclimatic
  conditions based on different mosquito habitat assumptions
