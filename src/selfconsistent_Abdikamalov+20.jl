using Gradus
using Plots
using Printf

_format_metric(m::JohannsenMetric) = Printf.@sprintf "a=%.3f, α13=%.2f" m.a m.α13
_format_metric(m::KerrMetric) = Printf.@sprintf "a=%.3f" m.a

_format_label(edd) = Printf.@sprintf "Ṁ / Ṁedd = %.1f" (edd / 100)

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
        sampler = EvenSampler(BothHemispheres(), GoldenSpiralGenerator())
    )
    return em_prof
end



function run_all_parameter_combinations(m, θ, bins; kwargs...)

    x = SVector(0.0, 1000.0, deg2rad(θ), 0.0)

    model = LampPostModel(h = 10.0)
    spectrum = PowerLawSpectrum(3.0)

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

    return (; metric = m, f = [edrat0, edrat10, edrat20, edrat30], θ = θ, bins = bins)

end

function plot_all(data)
    incl_text = Printf.@sprintf " θ=%0.f" data.θ
    p = plot(title = _format_metric(data.metric) * incl_text, legend = :topleft)
    for (edd, f) in zip((0, 10, 20, 30), data.f)
        plot!(p, data.bins, f, label = _format_label(edd))
    end
    p
end

# define custom bins for g
bins = collect(range(0.1, 1.4, 200))

INCLINATION = 70.0
KWARGS = (; β₀ = 2)

################## negative α13


m_n1 = JohannsenMetric(M = 1.0, a = 0.0, α13 = -0.35, ϵ3 = 0.0)
data_n1 = run_all_parameter_combinations(m_n1, INCLINATION, bins; KWARGS...)

pn1 = plot_all(data_n1)
display(pn1)

m_n2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = -0.35, ϵ3 = 0.0)
data_n2 = run_all_parameter_combinations(m_n2, INCLINATION, bins; KWARGS...)


pn2 = plot_all(data_n2)
display(pn2)

m_n3 = JohannsenMetric(M = 1.0, a = 0.998, α13 = -0.35, ϵ3 = 0.0)
data_n3 = run_all_parameter_combinations(m_n3, INCLINATION, bins; KWARGS...)

pn3 = plot_all(data_n3)
display(pn3)

################## zero α13

m_k1 = JohannsenMetric(M = 1.0, a = 0.0, α13 = 0.0, ϵ3 = 0.0)
data_k1 = run_all_parameter_combinations(m_k1, INCLINATION, bins; KWARGS...)

pk1 = plot_all(data_k1) 
display(pk1)

m_k2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.0, ϵ3 = 0.0)
data_k2 = run_all_parameter_combinations(m_k2, INCLINATION, bins; KWARGS...)

pk2 = plot_all(data_k2) 
display(pk2)

m_k3 = JohannsenMetric(M = 1.0, a = 0.998, α13 = 0.0, ϵ3 = 0.0) # spike at 30% Eddington ratio
data_k3 = run_all_parameter_combinations(m_k3, INCLINATION, bins; KWARGS...)

pk3 = plot_all(data_k3) 
display(pk3)

################## positive α13

m_p1 = JohannsenMetric(M = 1.0, a = 0.0, α13 = 0.35, ϵ3 = 0.0)
data_p1 = run_all_parameter_combinations(m_p1, INCLINATION, bins; KWARGS...)

pp1 = plot_all(data_p1) 
display(pp1)

m_p2 = JohannsenMetric(M = 1.0, a = 0.9, α13 = 0.35, ϵ3 = 0.0)
data_p2 = run_all_parameter_combinations(m_p2, INCLINATION, bins; KWARGS...)

pp2 = plot_all(data_p2) 
display(pp2)

m_p3 = JohannsenMetric(M = 1.0, a = 0.998, α13 = 0.35, ϵ3 = 0.0)
data_p3 = run_all_parameter_combinations(m_p3, INCLINATION, bins; KWARGS...)

pp3 = plot_all(data_p3) 
display(pp3)

# put everything together
plot(pn1, pk1, pp1, pn2, pk2, pp2, pn3, pk3, pp3, layout = grid(3, 3), size = (1100, 1100))




# begin
#     import Pkg; Pkg.add("JLD2")

#     using JLD2

#     all_data_names = [:data_n1, :data_n2, :data_n3, :data_k1, :data_k2, :data_k3, :data_p1, :data_p2, :data_p3]
#     data_to_save = [data_n1, data_n2, data_n3, data_k1, data_k2, data_k3, data_p1, data_p2, data_p3]

#     data = Dict(i => j for (i, j) in zip(all_data_names, data_to_save))

#     save("binning-method-full.jld2", 
#     data
#     )
# end