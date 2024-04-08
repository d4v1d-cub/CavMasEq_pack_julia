# This script contains some default rates for different models that could be useful for the user
# It contains also some functions that compute global observables that might be needed inside the rates
# or in the data post-processing by the user

# This struct is passed to the build_rate function at every step of the integration of the CME
mutable struct State_CME
    p_cav::Array{Float64, 4}
    probi::Vector{Float64}
    pu::Array{Float64, 3}
end


# This struct is passed to the build_rate function at every step of the integration of the CDA
mutable struct State_CDA
    p_joint::Matrix{Float64}
    pu::Array{Float64, 3}
end


# This function computes the energy in a K-SAT formula, where there is only one unsatisfied
# configuration of the variables in a clause
function ener(graph::HGraph, probi::Vector{Float64}, pu::Array{Float64, 3}, ch_u::Vector{Int64})
    e = 0.0
    for he in 1:graph.M
        node = graph.he_2_var[he, 1]
        bit_node = ch_u[he] & 1
        e += (bit_node + (1 - 2 * bit_node) * probi[node]) * pu[he, 1, 2]
    end
    return e
end


# This function computes the energy in a K-SAT formula, where there is only one unsatisfied
# configuration of the variables in a clause
function ener(p_joint::Matrix{Float64}, ch_u::Vector{Int64})
    e = 0.0
    for he in 1:size(p_joint, 1)
        e += p_joint[he, ch_u[he] + 1]
    end
    return e
end

# This is the rate of FMS algorithm for KSAT
function rate_FMS_KSAT(Ep::Int64, Em::Int64, eta::Float64, avE::Float64, K::Number, N::Int64)
    dE = Em - Ep
    if dE > 0
        return Ep * N / K / avE * eta ^ dE
    else
        return Ep * N / K / avE
    end
end


# Each rate function needs some specific arguments. The general form of these functions 
function build_args_rate_FMS(graph::HGraph, st::State_CME, ch_u::Vector{Int64}, eta::Float64)
    avE = ener(graph, st.probi, st.pu, ch_u)
    return eta, avE, graph.K, graph.N
end


# Each rate function needs some specific arguments. The general form of these functions 
function build_args_rate_FMS(graph::HGraph, st::State_CDA, ch_u::Vector{Int64}, eta::Float64)
    avE = ener(st.p_joint, ch_u)
    return eta, avE, graph.K, graph.N
end