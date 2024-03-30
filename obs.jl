include("graphs.jl")


function ener(graph::HGraph, probi::Vector{Float64}, pu::Array{Float64, 3}, ch_u::Vector{Int64})
    e = 0.0
    for he in 1:graph.M
        node = graph.he_2_var[he, 1]
        bit_node = ch_u & 1
        e += (bit_node + (1 - 2 * bit_node) * probi[node]) * pu[he, 1, 2]
    end
    return e
end