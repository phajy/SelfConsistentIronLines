# Simple line profile based on the Gradus documentation
using Gradus
using Plots

d = ThinDisc(0.0, 400.0)
x = SVector(0.0, 1000.0, deg2rad(30), 0.0)
m = KerrMetric(1.0, 0.998)

# maximal integration radius
maxrₑ = 50.0

# emissivity function
ε(r) = r^(-3)

# g grid to do flux integration over
gs = range(0.0, 1.2, 500)
_, flux = lineprofile(gs, ε, m, x, d, maxrₑ = maxrₑ, verbose = true)

# plot flux as a function of energy
plot(gs, flux, legend=false)
