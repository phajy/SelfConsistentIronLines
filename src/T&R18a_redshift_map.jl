using Gradus, Plots

# metric and metric parameters
m = KerrMetric(M=1.0, a=0.0)
# observer position
x = SVector(0.0, 1000.0, deg2rad(60), 0.0)
# accretion disc
d = ShakuraSunyaev(m, eddington_ratio = 0.1)
#d = ThinDisc(Gradus.isco(m), 30.0)

# define point function which filters geodesics that intersected the accretion disc
# and use those to calculate redshift
# filter the outer radius -- specify by adjusting value by gp.x[2]

r_filter = Gradus.FilterPointFunction((m, gp, t) -> gp.x[2] < 30, NaN)
pf = ConstPointFunctions.redshift(m, x) ∘ r_filter ∘ ConstPointFunctions.filter_intersected()

α, β, img = rendergeodesics(
    m,
    x,
    d,
    # maximum integration time
    2000.0,
    αlims = (-35, 35), 
    βlims = (-35, 35),
    image_width = 800,
    image_height = 800,
    verbose = true,
    pf = pf,
    # filter out the phantom image inside the shadow
    callback = domain_upper_hemisphere()
)

heatmap(α, β, img, aspect_ratio = 1, colour = :gist_rainbow, title = "Eddington ratio 10%")
