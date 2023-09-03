# SpectralResampling

[![Build Status](https://github.com/Michael-Reefe/SpectralResampling.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Michael-Reefe/SpectralResampling.jl/actions/workflows/CI.yml?query=branch%3Amain)

## Introduction

This Julia package allows one to resample fluxes and errors onto a new arbitrary wavelength grid while conserving the total flux. It is based off
of the Python package [SpectRes](https://github.com/ACCarnall/SpectRes) by Adam Carnall and implements most of the same functionality, but with
some minor differences (noted below). Please see [Carnall (2017)](https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract) for more information
about the theory behind the spectral resampling procedure.

## Installation

```julia
julia> using Pkg
julia> Pkg.add("SpectralResampling")
```

or alternatively,
```julia
julia> ]
(Environment)> add SpectralResampling
```

## Requirements

This package only has one dependency:

- [LoopVectorization](https://github.com/JuliaSIMD/LoopVectorization.jl)

## Usage

The main function, which is a translation of SpectRes's `spectres` function, is called `resample_conserving_flux` and has the following call signature:

```julia
resample_conserving_flux(new_wave, old_wave, flux[, err, mask]; fill=NaN, verbose=true)
```

### Required arguments:
- `new_wave::AbstractVector`: The new 1D wavelength array that the flux should be resampled onto.
- `old_wave::AbstractVector`: The original 1D wavelength array the the flux vector is currently sampled onto.
- `flux::AbstractArray`: An n-dimensional flux array giving the fluxes at the specified wavelengths. The first axis of
    this array should match the length of `old_wave`, while any additional axes will be treated as independent spectra.

### Optional arguments:
- `err::AbstractArray`: An optional n-dimension error array giving the errors in the fluxes at
    the specified wavelengths. The size of this array must match the size of `flux`.
- `mask::BitArray`: An optional n-dimensional mask array containing boolean values
    indicating whether certain flux values should be masked out. The size of this array must match the size of `flux`.
    Note that the mask is NOT used in the calculation of the new fluxes -- it is instead resampled in such a way that
    the output mask will cover the same pixels that the input mask covered.

### Keyword arguments:
- `fill::Real=NaN`: The fill value that is used to populate bins that are outside of the original wavelength range.
- `verbose::Bool=true`: Whether or not to print a warning message if any bins fall outside of the original wavelength range.

### Returns:
- `new_fluxes::AbstractArray`: The fluxes that have been resampled onto the output wavelength grid. The first axis will
    have the same length as `new_wave` whereas all other axes will be the same as the original `flux` array.

### Optional returns:
- `new_errs::AbstractArray`: Only returned if an input `err` is provided. These are the errors that have been resampled onto
    the output wavelength grid. Will be the same size as `new_fluxes`.
- `new_mask::BitArray`: Only returned if an input `mask` is provided. These are the new mask values that have been resampled
    onto the output wavelength grid. Will be the same size as `new_fluxes`.

A few other utility functions are provided that were not present in the original Python version. Namely, `get_logarithmic_λ` and 
`get_linear_λ`, which convert between linearly-spaced and logarithmically-spaced wavelength vectors. They each take as input a wavelength
vector (for `get_logarithmic_λ` the input must be linear and for `get_linear_λ` the input must be logarithmic), and an optional second
argument specifying the spacing of the new vector in log space or linear space. If no spacing is provided, the output vector is spaced
such that it has the same number of samples as the input vector.

## Differences from the python package

There are a few key differences worth mentioning:

- This code assumes that your flux array's FIRST axis is the wavelength axis, whereas the python version assumes it's the LAST axis.
  Make sure you're aware of the convention used by whichever version you're using! (note: this was changed for performance reasons,
  since Julia uses column-major indexing whereas python uses row-major indexing, so it's more optimal in Julia to make the first axis the one
  that is changing most rapidly).
- There is an additional `mask` argument that allows one to input a mask as a BitArray the same size as `flux` and `err`.
  The mask will be combined by taking an `any` operation along all of the bins in the original flux array that are
  combined into the new bin in the output flux array. That is to say, if any of the points used from the orignal flux
  array to calculate the rebinned flux were masked out, then the new flux value will also be masked out. The mask is NOT used 
  during the calculation of the output fluxes or errors.
- The `get_logarithmic_λ` and `get_linear_λ` functions are new (see above).
- Due to Julia's JIT compiler and optimizations using LoopVectorization.jl, this package is must faster than the original Python version,
  and even beats the numba-jitted Python version as well (after compilation). For very large 70000x40x40 flux and error arrays, during my testing,
  the Julia version completed in ~10s while the numba veresion took ~12s.  For smaller 70000x10x10 arrays, the julia version took ~0.07s versus
  ~0.17s for the numba version.

## Citation

If you use this code in your research, consider citing the original [Carnall (2017)](https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract)
paper.

