using SpectralResampling

# NOTE: DelimitedFiles and PyPlot are not dependencies of SpectralResampling.jl, make sure you have them installed!
using DelimitedFiles  
using PyPlot
using LaTeXStrings

# Load in the data from the text file
data = readdlm("VST-ATLAS_J025.6821-33.4627.txt", ' ', Float64, '\n', comments=true)

# Specify the new grid to be rebinned to
new_wave = (6510.:5.:9380.) .+ 2.5

# Call the resampling function with the flux and errors
new_flux, new_err = resample_conserving_flux(new_wave, data[:,1], data[:,2], data[:,3])

# Plotting code
f, (ax1, ax2) = plt.subplots(2, figsize=(15,7))
gs = plt.matplotlib.gridspec.GridSpec(2, 1, height_ratios=[3,1])

ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])

ax1.plot(data[:,1], data[:,2] .* 1e19, color="blue", lw=1.5, label=L"1 $\mathrm{\AA}$ Sampling")
ax1.plot(new_wave, new_flux .* 1e19, color="red", lw=1.5, label=L"5 $\mathrm{\AA}$ Sampling")

ax2.plot(data[:,1], data[:,3] .* 1e19, color="blue", lw=1.5)
ax2.plot(new_wave, new_err .* 1e19, color="red", lw=1.5)

ax1.set_ylabel(L"Flux ($10^{-19}\ \mathrm{W/m^2/\AA}$)", size=18)
ax2.set_ylabel("Flux Error", size=18)

ax2.set_xlabel(L"Wavelength $(\mathrm{\AA})$", size=18)

ax1.set_xlim(6800, 9000)
ax1.set_ylim(-0.15, 0.65)
ax2.set_xlim(6800, 9000)
ax2.set_ylim(0, 0.08)
ax1.legend()

plt.show()
