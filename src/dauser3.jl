using Gradus
using Plots

bins = collect(range(0.2, 1.15, 300))
plane = PolarPlane(GeometricGrid(); Nr = 3000, Nθ = 1000, r_max = 1000.0)


x1=SVector(0.0, 1000.0, deg2rad(20),0.0)
x2=SVector(0.0, 1000.0, deg2rad(30),0.0)
x3=SVector(0.0, 1000.0, deg2rad(40),0.0)
x4=SVector(0.0, 1000.0, deg2rad(60),0.0)

d=ThinDisc(0.0, Inf)


m1=KerrMetric(M=1.0, a=0.998)

m2=KerrMetric(M=1.0, a=0.000)
m3=KerrMetric(M=1.0, a=-0.998)


function calculate_line_profile(m, x, d, bins; kwargs...)
    _, f = lineprofile(
        m,
        x,
        d;
        method = TransferFunctionMethod(),
        minrₑ = Gradus.isco(m) + 1e-2,
        
        verbose = true,
        bins = bins,
        #minrₑ = 9*1.23,
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
function plot_data(p1,p2,p3)
    plot(bins, p1, label = "a=0.998", title = "Inclination=40, R_in=1.23, r_out=50")
    plot!(bins, p2, label = "a=0")
    plot!(bins, p3, label = "a=-0.998")
    
end
b1=calculate_line_profile(m1,x3,d,bins)
b2=calculate_line_profile(m2,x3,d,bins)
b3=calculate_line_profile(m3,x3,d,bins)

data_40=plot_data(b1,b2,b3)
display(data_40)