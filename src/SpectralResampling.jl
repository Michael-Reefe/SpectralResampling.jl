"""
SpectralResampling.jl - Resample fluxes and errors onto a new arbitrary wavelength grid.
Copyright (C) 2022 Adam Carnall (original Python version)
Copyright (C) 2023 Michael Reefe (Julia version)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""


module SpectralResampling

using LoopVectorization

export get_logarithmic_λ, get_linear_λ, make_bins, resample_conserving_flux


"""
    get_logarithmic_λ(λ[, log_spacing])

Calculate a logarithmically-spaced wavelength grid given a linearly spaced wavelength grid `λ`, optionally
specifying the spacing (in log space) between elements `log_spacing`. For example, if `lnλ` is a logarithmically
spaced wavelength grid, the `log_spacing` would be defined as `log_spacing = log(lnλ[2]/lnλ[1]) = log(lnλ[end]/lnλ[end-1])`.
"""
function get_logarithmic_λ(λ::AbstractVector{<:Real}, log_spacing::Union{Real,Nothing}=nothing)
    diffs = diff(λ)
    @assert diffs[end] ≈ diffs[1] "The input wavelength vector must be linearly spaced!"

    if isnothing(log_spacing)
        # use the 'length' argument to ensure the output is the same length as the input
        return exp.(range(log(λ[1]), log(λ[end]), length=length(λ)))
    end
    # use the 'step' argument to ensure the spacing between elements is given by 'log_spacing'
    return exp.(range(log(λ[1]), log(λ[end]), step=log_spacing))
end


"""
    get_linear_λ(lnλ[, lin_spacing])

Does the opposite of `get_logarithmic_λ`, returning a linearly spaced wavelength vector given a logarithmically
spaced one. Optionally specify the spacing in linear space with the `lin_spacing` argument, defined such that
`lin_spacing = λ[2] - λ[1] = λ[end] - λ[end-1]`.
"""
function get_linear_λ(lnλ::AbstractVector{<:Real}, lin_spacing::Union{Real,Nothing}=nothing)
    @assert lnλ[2]/lnλ[1] ≈ lnλ[end]/lnλ[end-1] "The input wavelength vector must be logarithmically spaced! (natural log)"

    if isnothing(lin_spacing)
        # use the 'length' argument to ensure the output is the same length as the input
        return range(lnλ[1], lnλ[end], length=length(lnλ))
    end
    # use the 'step' argument to ensure the spacing between the elements is given by 'lin_spacing'
    return range(lnλ[1], lnλ[end], step=lin_spacing)
end


"""
    make_bins(vec)

Calculate the bin edges and bin widths for a vector given the distances between each entry

ADDITIONAL NOTES: 

    - This function has been adapted from the SpectRes python package, https://github.com/ACCarnall/SpectRes,
      by Adam Carnall. Please refer to Carnall (2017): https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract
      for an explanation of the theory behind the spectral resampling procedure.  This function is a translation
      of `make_bins` from SpectRes and performs the same functions.
"""
function make_bins(vec::AbstractVector)

    # Get the bin edges 
    edges = zeros(eltype(vec), length(vec)+1)
    edges[1] = vec[1] - (vec[2] - vec[1])/2
    edges[end] = vec[end] + (vec[end] - vec[end-1])/2
    edges[2:end-1] .= (vec[2:end] + vec[1:end-1])/2

    # Get the bin widths
    widths = zeros(eltype(vec), length(vec))
    widths[end] = vec[end] - vec[end-1]
    widths[1:end-1] = edges[2:end-1] - edges[1:end-2]

    # Return the edges and widths
    edges, widths
end


"""
    resample_conserving_flux(new_wave, old_wave, flux, err=nothing, mask=nothing; fill=NaN)

Resample flux (and optionally errors and a mask) onto a new wavelength grid, while convserving the flux.

# Required Arguments
- `new_wave::AbstractVector`: The new 1D wavelength array that the flux should be resampled onto.
- `old_wave::AbstractVector`: The original 1D wavelength array the the flux vector is currently sampled onto.
- `flux::AbstractArray`: An n-dimensional flux array giving the fluxes at the specified wavelengths. The first axis of
    this array should match the length of `old_wave`, while any additional axes will be treated as independent spectra.

# Optional Arguments
- `err::AbstractArray`: An optional n-dimension error array giving the errors in the fluxes at
    the specified wavelengths. The size of this array must match the size of `flux`.
- `mask::BitArray`: An optional n-dimensional mask array containing boolean values
    indicating whether certain flux values should be masked out. The size of this array must match the size of `flux`.
    Note that the mask is NOT used in the calculation of the new fluxes -- it is instead resampled in such a way that
    the output mask will cover the same pixels that the input mask covered.

# Keyword Arguments
- `fill::Real=NaN`: The fill value that is used to populate bins that are outside of the original wavelength range.
- `verbose::Bool=true`: Whether or not to print a warning message if any bins fall outside of the original wavelength range.

# Returns
- `new_fluxes::AbstractArray`: The fluxes that have been resampled onto the output wavelength grid. The first axis will
    have the same length as `new_wave` whereas all other axes will be the same as the original `flux` array.

# Optional returns
- `new_errs::AbstractArray`: Only returned if an input `err` is provided. These are the errors that have been resampled onto
    the output wavelength grid. Will be the same size as `new_fluxes`.
- `new_mask::BitArray`: Only returned if an input `mask` is provided. These are the new mask values that have been resampled
    onto the output wavelength grid. Will be the same size as `new_fluxes`.

ADDITIONAL NOTES: 

    - This module has been adapted from the SpectRes python package, https://github.com/ACCarnall/SpectRes,
      by Adam Carnall. Please refer to Carnall (2017): https://ui.adsabs.harvard.edu/abs/2017arXiv170505165C/abstract
      for an explanation of the theory behind the spectral resampling procedure. This function should perform
      identically to the `spectres.spectres` function from the original python version. However, there are some 
      minor differences enumerated below.

    - This version uses a different convention by assuming the FIRST axis of the flux array is the wavelength axis 
      rather than the last axis in the python version (this was done for performance reasons, since julia is column-major 
      whereas python is row-major)

    - There is an additional `mask` argument that allows one to input a mask as a BitArray the same size as `flux` and `err`.
      The mask will be combined by taking an `any` operation along all of the bins in the original flux array that are
      combined into the new bin in the output flux array. That is to say, if any of the points used from the orignal flux
      array to calculate the rebinned flux were masked out, then the new flux value will also be masked out.

    - Care has been taken to improve the performance over the original python SpectRes package. This version is even faster than
      the numba-jitted version from the original package.
"""
function resample_conserving_flux(new_wave::AbstractVector, old_wave::AbstractVector, flux::AbstractArray; fill::Real=NaN, 
    verbose::Bool=true)::AbstractArray
    new_flux = _resample_conserving_flux_internal(new_wave, old_wave, flux; fill=fill, verbose=verbose)

    new_flux
end

function resample_conserving_flux(new_wave::AbstractVector, old_wave::AbstractVector, flux::AbstractArray, err::AbstractArray;
    fill::Real=NaN, verbose::Bool=true)::Tuple{AbstractArray, AbstractArray}
    new_flux, new_err = _resample_conserving_flux_internal(new_wave, old_wave, flux, err; fill=fill, verbose=verbose)

    new_flux, new_err
end

function resample_conserving_flux(new_wave::AbstractVector, old_wave::AbstractVector, flux::AbstractArray, mask::BitArray;
    fill::Real=NaN, verbose::Bool=true)::Tuple{AbstractArray, AbstractArray}
    new_flux, new_mask = _resample_conserving_flux_internal(new_wave, old_wave, flux, nothing, mask; fill=fill, verbose=verbose)

    new_flux, new_mask
end

function resample_conserving_flux(new_wave::AbstractVector, old_wave::AbstractVector, flux::AbstractArray,
    err::AbstractArray, mask::BitArray; fill::Real=NaN, verbose::Bool=true)::Tuple{AbstractArray,AbstractArray,AbstractArray}
    new_flux, new_err, new_mask = _resample_conserving_flux_internal(new_wave, old_wave, flux, err, mask; fill=fill, verbose=verbose)

    new_flux, new_err, new_mask
end

function _resample_conserving_flux_internal(new_wave::AbstractVector, old_wave::AbstractVector, flux::AbstractArray, 
    err::Union{AbstractArray,Nothing}=nothing, mask::Union{AbstractArray,Nothing}=nothing; fill::Real=NaN, verbose::Bool=true)

    # Do some input checking
    @assert length(old_wave) == size(flux, 1) "The first axis of the flux array must be the same length as the original wavelength array!"

    # Find the edges and widths of the new wavelength bins given the old ones
    old_edges, old_widths = make_bins(old_wave)
    new_edges, new_widths = make_bins(new_wave)

    # Instantiate flux and error arrays with zeros
    new_fluxes = zeros(eltype(flux), (length(new_wave), size(selectdim(flux, 1, 1))...))
    if !isnothing(err)
        @assert size(err) == size(flux) "error must be the same size as flux!"
        new_errs = copy(new_fluxes)
    end
    if !isnothing(mask)
        @assert size(mask) == size(flux) "mask must be the same size as flux!"
        new_mask = falses(size(new_fluxes))
    end

    start = 1
    stop = 1
    warned = false

    # Helper functions
    function accumulate_flux_turbo!(new_fluxes::AbstractArray, flux::AbstractArray, 
        old_widths::AbstractVector, start::Integer, stop::Integer, j::Integer; quad::Bool=false)
        newflux_slice = selectdim(new_fluxes, 1, j)
        for ind ∈ CartesianIndices(newflux_slice)
            num_type = typeof(quad ? (new_fluxes[1]*old_widths[1])^2 : new_fluxes[1]*old_widths[1])
            num = zero(num_type)
            denom = zero(eltype(old_widths))
            f_slice = @view flux[:, ind]
            @turbo for k ∈ start:stop
                q = old_widths[k] * f_slice[k]
                # 'quad' argument combines in quadrature -- to be used if the inputs are errors
                num += quad ? q^2 : q
                denom += old_widths[k]
            end
            num = quad ? √(num) : num
            newflux_slice[ind] = num/denom
        end
    end

    # Calculate new flux and uncertainty values, looping over new bins
    for j ∈ eachindex(new_wave)

        # Add filler values if new_wavs extends outside of spec_wavs
        if (new_edges[j] < old_edges[1]) || (new_edges[j+1] > old_edges[end])
            # since selectdim returns a view, modifying `newflux_slice` also modifies `new_fluxes`
            newflux_slice = selectdim(new_fluxes, 1, j)
            newflux_slice .= fill    
            if !isnothing(err)
                newerr_slice = selectdim(new_errs, 1, j)
                newerr_slice .= fill
            end
            if !isnothing(mask)
                newmask_slice = selectdim(new_mask, 1, j)
                newmask_slice .= 1
            end
            if (j == 1 || j == size(new_wave, 1)) && !warned && verbose
                println("\nnew_wave contains values outside the range " *
                      "in old_wave, new_fluxes and new_errs will be filled " *
                      "with the value set in the 'fill' keyword argument. \n")
                warned = true
            end
            continue
        end
    
        # Find first old bin which is partially covered by the new bin
        while old_edges[start+1] ≤ new_edges[j]
            start += 1
        end

        # Find last old bin which is partially covered by the new bin
        while old_edges[stop+1] < new_edges[j+1]
            stop += 1
        end

        # If new bin is fully inside an old bin start and stop are equal
        if start == stop
            newflux_slice = selectdim(new_fluxes, 1, j)
            newflux_slice .= selectdim(flux, 1, start)
            if !isnothing(err)
                newerr_slice = selectdim(new_errs, 1, j)
                newerr_slice .= selectdim(err, 1, start)
            end
            if !isnothing(mask)
                newmask_slice = selectdim(new_mask, 1, j)
                newmask_slice .= selectdim(mask, 1, start)
            end
        # Otherwise multiply the first and last old bin widths by P_ij
        else 
            start_factor = (old_edges[start+1] - new_edges[j]) / (old_edges[start+1] - old_edges[start])
            stop_factor = (new_edges[j+1] - old_edges[stop]) / (old_edges[stop+1] - old_edges[stop])

            old_widths[start] *= start_factor
            old_widths[stop] *= stop_factor

            # Populate new fluxes and errors
            accumulate_flux_turbo!(new_fluxes, flux, old_widths, start, stop, j)

            # Combine errors in quadrature
            if !isnothing(err)
                accumulate_flux_turbo!(new_errs, err, old_widths, start, stop, j, quad=true)
            end

            # Combine the mask data with an 'any' operation along the wavelength axis
            # such that if any of the pixels in the bins used to calculate the new flux are masked,
            # then the output bin will also be masked
            if !isnothing(mask)
                newmask_slice = selectdim(new_mask, 1, j)
                newmask_slice .= dropdims(any(selectdim(mask, 1, start:stop), dims=1), dims=1)
            end

            # Put the old bin widths back to their initial values
            old_widths[start] /= start_factor
            old_widths[stop] /= stop_factor

        end
    end

    # Only return the quantities that were given in the input
    if !isnothing(err) && !isnothing(mask)
        new_fluxes, new_errs, new_mask
    elseif !isnothing(err)
        new_fluxes, new_errs
    elseif !isnothing(mask)
        new_fluxes, new_mask
    else
        new_fluxes
    end

end


end
