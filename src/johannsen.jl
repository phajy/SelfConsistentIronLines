using Gradus
using Plots

#geometric thin disc
d = ThinDisc(0.0, Inf)

m1 = JohannsenMetric(M = 1.0, a = 0.95, α22 = 0.0)
m2 = JohannsenMetric(M = 1.0, a = 0.85, α22 = 0.97902)

# define custom bins for g
bins = collect(range(0.1, 2.0, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)

model = LampPostModel(h = 10.0)
spectrum = PowerLawSpectrum(3.0)

em_prof1 = Gradus.emissivity_profile(
    m1,
    d,
    model,
    spectrum,
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

em_prof2 = Gradus.emissivity_profile(
    m2,
    d,
    model,
    spectrum,
    n_samples = 10_000,
    sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator()),
)

prof1 = Gradus.RadialDiscProfile(em_prof1)
prof2 = Gradus.RadialDiscProfile(em_prof2)

function calculate_line_profile(m, x, d, prof, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d,
        prof;
        method = TransferFunctionMethod(),
        minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        #minrₑ = 6.0,
        maxrₑ = 100.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        # abstol = 1e-10,
        # reltol = 1e-10,
        kwargs...,
    )
    return f
end

# self-consistently calculated disc illumination profile
# edrat1 = calculate_line_profile(m1, x, d, prof1, bins)
# edrat2 = calculate_line_profile(m2, x, d, prof2, bins)

function calculate_line_profile(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        #minrₑ = 6.0,
        maxrₑ = 100.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        # abstol = 1e-10,
        # reltol = 1e-10,
        kwargs...,
    )
    return f
end

# r^-3 emissivity profile used by Johannsen (2014) - see their description just beofore their equation (38) 
# note that lineprofile assumes emissivity goes as r^-3 by default
edrat1 = calculate_line_profile(m1, x, d, bins)
edrat2 = calculate_line_profile(m2, x, d, bins)


# function min_max_scaling(data::Vector{T}) where T
#     min_val = minimum(data)
#     max_val = maximum(data)
#     scaled_data = [(x - min_val) / (max_val - min_val) for x in data]
#     return scaled_data
# end


# scaled_edrat1 = min_max_scaling(edrat1)
# scaled_edrat2 = min_max_scaling(edrat2)

# plot(bins, scaled_edrat1, xlabel = "energy (a.u.)", ylabel = "number flux density (a.u.)", label = "a = 0.95, α13 = 0.0", title = "Scaled")
# plot!(bins, scaled_edrat2, xlabel ="energy (a.u.)", ylabel = "number flux density (a.u.)", label = "a = 0.85, α13 = 0.97902")

plot(bins, edrat1, xlabel = "energy (a.u.)", ylabel = "number flux density (a.u.)", label = "a = 0.95, α22 = 0.0", title = "Johannsen 2014 fig. 9 centre panel recipe", xlims = (0.2, 1.1))
plot!(bins, edrat2, xlabel ="energy (a.u.)", ylabel = "number flux density (a.u.)", label = "a = 0.85, α22 = 0.97902")
