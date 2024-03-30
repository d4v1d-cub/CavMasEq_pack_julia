include("graphs.jl")
include("obs.jl")

# This is the rate of FMS algorithm for KSAT
function rate_FMS_KSAT(Ep::Int64, Em::Int64, T::Float64, avE::Float64, K::Number)
    dE = Em - Ep
    if dE > 0
        return Ep / K / avE * exp(-dE / T)
    else
        return Ep / K / avE
    end
end

function build_args_rate_FMS(graph::HGraph, p_cav::Array{Float64, 4}, probi::Vector{Float64},
                             pu::Array{Float64, 3}, ch_u::Matrix{Int64}, T::Float64)
    avE = ener(graph, probi, pu, ch_u)
    return T, avE, graph.K
end