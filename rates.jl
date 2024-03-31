# This script contains some default rates for different models that could be useful for the user
# It contains also some functions that compute global observables that might be needed inside the rates
# or in the data post-processing by the user

mutable struct State
    p_cav::Array{Float64, 4}
    probi::Vector{Float64}
    pu::Array{Float64, 3}
end


# This function computes the energy in a K-SAT formula, where there is only one unsatisfied
# configuration of the variables in a clause
function ener(graph::HGraph, probi::Vector{Float64}, pu::Array{Float64, 3}, ch_u::Vector{Int64})
    e = 0.0
    for he in 1:graph.M
        node = graph.he_2_var[he, 1]
        bit_node = ch_u & 1
        e += (bit_node + (1 - 2 * bit_node) * probi[node]) * pu[he, 1, 2]
    end
    return e
end

# This is the rate of FMS algorithm for KSAT
function rate_FMS_KSAT(Ep::Int64, Em::Int64, T::Float64, avE::Float64, K::Number)
    dE = Em - Ep
    if dE > 0
        return Ep / K / avE * exp(-dE / T)
    else
        return Ep / K / avE
    end
end


# Each rate function needs some specific arguments. The general form of these functions 
function build_args_rate_FMS(graph::HGraph, st::State, ch_u::Matrix{Int64}, T::Float64)
    avE = ener(graph, st.probi, st.pu, ch_u)
    return T, avE, graph.K
end