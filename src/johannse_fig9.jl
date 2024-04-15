using Gradus, Plots, LaTeXStrings
using Plots.PlotMeasures

#geometric thin disc
d = ThinDisc(0.0, Inf)

m1 = JohannsenMetric(M = 1.0, a = 0.95, α13 = 0.0)
m2 = JohannsenMetric(M = 1.0, a = 0.85, α13 = -1.07115)
m3=JohannsenMetric(M=1.0,a=0.95,α22=0.0)
m4=JohannsenMetric(M=1.0,a=0.85,α22=0.97902)
m5=JohannsenMetric(M=1.0,a=0.95,α22=0.0)
m6=JohannsenMetric(M=1.0,a=0.8,α22=2.12516)


# define custom bins for g
bins = collect(range(0.1, 2.0, 300))

# define the plane to perform the binning over
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)

# inclination
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)



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


edrat1 = calculate_line_profile(m1, x, d, bins)
edrat2 = calculate_line_profile(m2, x, d, bins)

edrat3 = calculate_line_profile(m3, x, d, bins)
edrat4 = calculate_line_profile(m4, x, d, bins)

edrat5 = calculate_line_profile(m5, x, d, bins)
edrat6 = calculate_line_profile(m6, x, d, bins)



p1=plot(bins,edrat1,color=:black,xlims=(0.1,1.2),xlabel = L"\textrm{Observed~Frequency~Shift~}(\nu_o/\nu_c)", ylabel = L"\textrm{Flux~(arbitrary~units)}",label=L"\textrm{a = 0.95, α13 = 0.0}")
p1=plot!(bins,edrat2,color=:tan1,label=L"\textrm{a = 0.85, α13 = -1.07115}")

p2=plot(bins,edrat3,color=:black,xlims=(0.1,1.2),xlabel = L"\textrm{Observed~Frequency~Shift~}(\nu_o/\nu_c)", ylabel = L"\textrm{Flux~(arbitrary~units)}",label=L"\textrm{a=0.95,α22=0.0}")
p2=plot!(bins,edrat4,color=:hotpink2,label=L"\textrm{a=0.85,α22=0.97902}")

p3=plot(bins,edrat5,color=:black,xlims=(0.1,1.2),xlabel = L"\textrm{Observed~Frequency~Shift~}(\nu_o/\nu_c)", ylabel = L"\textrm{Flux~(arbitrary~units)}",label=L"\textrm{a=0.95,α22=0.0}")
p3=plot!(bins,edrat6,color=:darkorchid1,label=L"\textrm{a=0.8,α22=2.12516}")





plot(p1,p2,p3,layout=(3,1),size=(600,1100),left_margin = [10mm 0mm])
