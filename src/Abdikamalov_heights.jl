using Gradus
using Plots
using Printf
using LaTeXStrings
using Plots.PlotMeasures

_format_metric(m::JohannsenMetric) = Printf.@sprintf "a=%.2f, α13=%.2f" m.a m.α13
_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.2f" m.a

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
        maxrₑ = 500.0,
        # resolution
        numrₑ = 200,
        Nr = 3000,
        abstol = 1e-10,
        reltol = 1e-10,
        kwargs...,
    )
    return f
end

function em_prof(m, d, model, spectrum; kwargs...)
    em_prof = Gradus.emissivity_profile(
        m,
        d,
        model,
        spectrum,
        n_samples = 10_000,
    )
    return em_prof
end



function run_all_parameter_combinations(m, θ, model, bins; kwargs...)

    x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)

    spectrum = PowerLawSpectrum(2.0)

    @info "m = $m θ = $(θ)"

    # discs
    dthin = ThinDisc(0.0, Inf)
    d1 = ShakuraSunyaev(m, eddington_ratio = 0.1)
    d2 = ShakuraSunyaev(m, eddington_ratio = 0.2)
    d3 = ShakuraSunyaev(m, eddington_ratio = 0.3)

    @info "0%"
    em_prof0 = em_prof(m, dthin, model, spectrum)
    prof0 = Gradus.RadialDiscProfile(em_prof0)
    edrat0 = @time calculate_line_profile(m, x, dthin, prof0, bins; kwargs...)
    @info "10%"
    em_prof10 = em_prof(m, d1, model, spectrum)
    prof10 = Gradus.RadialDiscProfile(em_prof10)
    edrat10 = @time calculate_line_profile(m, x, d1, prof10, bins; kwargs...)
    @info "20%"
    em_prof20 = em_prof(m, d2, model, spectrum)
    prof20 = Gradus.RadialDiscProfile(em_prof20)
    edrat20 = @time calculate_line_profile(m, x, d2, prof20, bins; kwargs...)
    @info "30%"
    em_prof30 = em_prof(m, d3, model, spectrum)
    prof30 = Gradus.RadialDiscProfile(em_prof30)
    edrat30 = @time calculate_line_profile(m, x, d3, prof30, bins; kwargs...)

    return (; metric = m, f = [edrat0, edrat10, edrat20, edrat30], θ = θ, bins = bins, h = model.h)

end

function plot_all(data)
    incl_text = Printf.@sprintf " θ=%0.f h=%.1f" data.θ data.h
    p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft, xlabel = L"\textrm{Observed~frequency~shift~} (\nu_o / \nu_e)", ylabel = L"\testrm{Flux~(arbitrary~units)")
    for (edd, f) in zip((0, 10, 20, 30), data.f)
        label_string = L"\dot{M} / \dot{N}_\textrm{edd} = " * string(edd / 100)
        plot!(p, data.bins, f, label = label_string)
    end
    p
end

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

INCLINATION = 70.0
# check for β = 2
KWARGS = (; β₀ = 2)

m_n1 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model5 = LampPostModel(h = 9.5)
data_n1 = run_all_parameter_combinations(m_n1, INCLINATION, model5, bins; KWARGS...)

pn1 = plot_all(data_n1)
display(pn1)

m_n2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model10 = LampPostModel(h = 9.6)
data_n2 = run_all_parameter_combinations(m_n2, INCLINATION, model10, bins; KWARGS...)


pn2 = plot_all(data_n2)
display(pn2)

m_n3 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
model15 = LampPostModel(h = 9.7)
data_n3 = run_all_parameter_combinations(m_n3, INCLINATION, model15, bins; KWARGS...)

pn3 = plot_all(data_n3)
display(pn3)


# put everything together
plot(pn1, pn2, pn3, layout = grid(1, 3, heights=[0.3]), size = (1200, 1200), left_margin = [5mm 0mm], right_margin = [5mm 0mm])
