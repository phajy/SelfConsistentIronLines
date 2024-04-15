using Gradus
using Plots, LaTeXStrings

#geometric thin disc
d = ThinDisc(0.0, Inf)

m1 = JohannsenMetric(M = 1.0, a = 0.4, α22=-2.0)
m2 = JohannsenMetric(M = 1.0, a = 0.4, α22 = 0.0)
m3=JohannsenMetric(M=1.0,a=0.4,α22=2.0)

m4=JohannsenMetric(M=1.0,a=0.8,α22=-2.0)
m5=JohannsenMetric(M=1.0,a=0.8,α22=0.0)
m6=JohannsenMetric(M=1.0,a=0.8,α22=2.0)

m7=JohannsenMetric(M=1.0,a=0.8,ϵ3=-2.0)
m8=JohannsenMetric(M=1.0,a=0.8,ϵ3=0.0)
m9=JohannsenMetric(M=1.0,a=0.8,ϵ3=2.0)

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


# self-consistently calculated disc illumination profile
edrat1 = calculate_line_profile(m1, x, d, bins)
edrat2 = calculate_line_profile(m2, x, d, bins)
edrat3 = calculate_line_profile(m3, x, d, bins)

edrat4 = calculate_line_profile(m4, x, d, bins)
edrat5 = calculate_line_profile(m5, x, d, bins)
edrat6 = calculate_line_profile(m6, x, d, bins)

#change inclination to 60
x1 = SVector(0.0, 1000.0, deg2rad(60), 0.0)
edrat7 = calculate_line_profile(m7, x1, d, bins)
edrat8 = calculate_line_profile(m8, x1, d, bins)
edrat9 = calculate_line_profile(m9, x1, d, bins)



p1=plot(bins,edrat1,color=:darkorchid1,xlims=(0.2,1.2),xlabel = L"\textrm{Observed~ Frequency~ Shift~}(\nu_o / \nu_e)", ylabel = L"\textrm{Flux~ (arbitrary~ units)}",label=L"\textrm{α22 = -2.0}",title=L"\textrm{a=0.4M}")
p1=plot!(bins,edrat2,color=:hotpink2,label=L"\textrm{α22 = 0.0}")
p1=plot!(bins,edrat3,color=:tan1,label=L"\textrm{α22=2.0}")

p2=plot(bins,edrat4,color=:darkorchid1,xlims=(0.2,1.2),xlabel = L"\textrm{Observed~ Frequency~ Shift}(\nu_o / \nu_e)", ylabel = L"\textrm{Flux~ (arbitrary~ units)}",label=L"\textrm{α22 = -2.0}",title=L"\textrm{a=0.8M}")
p2=plot!(bins,edrat5,color=:hotpink2,label=L"\textrm{α22=0.0}")
p2=plot!(bins,edrat6,color=:tan1,label=L"\textrm{α22=2.0}")

p3=plot(bins,edrat7,color=:darkorchid1,xlims=(0.2,1.4),xlabel = L"\textrm{Observed~ Frequency~ Shift}(\nu_o / \nu_e)", ylabel = L"\textrm{Flux~ (arbitrary ~units)}",label=L"\textrm{ϵ3=-2.0}",title=L"\textrm{a=0.8M}")
p3=plot!(bins,edrat8,color=:hotpink2,label=L"\textrm{ϵ3=0.0}")
p3=plot!(bins,edrat9,color=:tan1,label=L"\textrm{ϵ3=2.0}")





plot(p1,p2)
plot(p3)
