
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `mosmicrosim`: mosquito microclimate simulator

<!-- badges: start -->
<!-- badges: end -->

A prototype R package for simulating a timeseries of microclimatic
conditions - e.g. local water surface area, water and air temperature
and humidity - in a given location on earth, given a simple set of
assumptions about the species’ habitat preferences.

Most of the computational work is done by the [NicheMapR R
package](https://mrke.github.io/) and its internal Fortran routines for
physical simulation of microclimate conditions. This package provides a
simplified interface for specifying mosquito habitat types, linking
NicheMapR to custom spatial climate/weather data, and including some
additional simulation models for different types of small water bodies.

So far the package does nothing. The to do list is as follows:

- Import NicheMapR and set up package
- Write a shim for functions like `NicheMapR::microclimate()` to enable
  incorporation of other sources of climate timeseries data
- Streamline simulation of microclimate components not including local
  water surface area, and water and air temperature, and humidity
- Implement simulation of local water surface area, and water and air
  temperature, and humidity using a set of water body models (cone,
  cylinder, with/without leakage)
- Write a simplified interface to simulating microclimatic conditions
  based on different mosquito habitat assumptions
