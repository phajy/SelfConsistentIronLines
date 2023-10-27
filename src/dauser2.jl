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
        #minrₑ = Gradus.isco(m) + 1e-2,
        
        verbose = true,
        bins = bins,
        minrₑ = 9,
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
    plot(bins, p1, label = "a=0.998")
    plot!(bins, p2, label = "a=0")
    plot!(bins, p3, label = "a=-0.998")
    
end

p1=calculate_line_profile(m1,x1,d,bins)
p2=calculate_line_profile(m2,x1,d,bins)
p3=calculate_line_profile(m3,x1,d,bins)

data_20=plot_data(p1,p2,p3)


display(data_20)

a1=calculate_line_profile(m1,x2,d,bins)
a2=calculate_line_profile(m2,x2,d,bins)
a3=calculate_line_profile(m3,x2,d,bins)

data_30=plot_data(a1,a2,a3)


b1=calculate_line_profile(m1,x3,d,bins)
b2=calculate_line_profile(m2,x3,d,bins)
b3=calculate_line_profile(m3,x3,d,bins)

data_40=plot_data(b1,b2,b3)


c1=calculate_line_profile(m1,x4,d,bins)
c2=calculate_line_profile(m2,x4,d,bins)
c3=calculate_line_profile(m3,x4,d,bins)

data_60=plot_data(c1,c2,c3)


plot(data_20,data_30,data_40,data_60,layout=grid(2,2),size=(1100,1100))
