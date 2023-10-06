using Gradus
using Plots

# define custom bins for g
bins = collect(range(0.1, 1.4, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x1 = SVector(0.0, 1000.0, deg2rad(5), 0.0)
x2 = SVector(0.0, 1000.0, deg2rad(30), 0.0)
x3 = SVector(0.0, 1000.0, deg2rad(45), 0.0)
x4 = SVector(0.0, 1000.0, deg2rad(85), 0.0)

#geometric thin disc
d = ThinDisc(0.0, Inf)

m = KerrMetric(M = 1.0, a = 0.998)

function calculate_line_profile(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        #minrₑ = Gradus.isco(m) + 1e-2,
        verbose = true,
        bins = bins,
        minrₑ = 1.23,
        maxrₑ = 20.0,
        # resolution
        numrₑ = 350,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end


edrat1 = calculate_line_profile(m, x1, d, bins)
edrat2 = calculate_line_profile(m, x2, d, bins)
edrat3 = calculate_line_profile(m, x3, d, bins)
edrat4 = calculate_line_profile(m, x4, d, bins)

plot(bins, edrat1, label = "5 degrees", title = "a = 0.998, r_in = 1.23, r_out = 20.0")
plot!(bins, edrat2, label = "30 degrees")
plot!(bins, edrat3, label = "45 degrees")
plot!(bins, edrat4, label = "85 degrees")

