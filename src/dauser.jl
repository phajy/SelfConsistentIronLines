using Gradus
using Plots


# define custom bins for g
bins = collect(range(0.1, 1.25, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x1 = SVector(0.0, 1000.0, deg2rad(40), 0.0)


#geometric thin disc
d = ThinDisc(0.0, Inf)

m = KerrMetric(M=1.0, a = 0.998)

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
        maxrₑ = 50.0,
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


plot(bins, edrat1, label = "40 degrees", title = "a = 0.998, r_in = 1.23, r_out = 50.0")


