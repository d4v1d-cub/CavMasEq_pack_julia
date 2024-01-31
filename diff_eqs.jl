include("graphs.jl")

# This function initializes the cavity conditional probabilities with a fully 
# independent initial distribution given by the vector 'p0'
function init_p_cav(graph::HGraph, p0::Vector{Float64})
    ch_exc = graph.chains_he ÷ 2
    p_cav = zeros(Float64, (graph.M, graph.K, ch_exc))
    for he in 1:graph.M
        for i in 1:graph.K
            for ch in 0:ch_exc - 1
                bits = digits(ch, base=2, pad=graph.K-1)
                p_cav[he, i, ch + 1] = prod(bits + (1 .- 2 * bits) .* p0[graph.nodes_except[he, i, :]])
            end
        end
    end
    return p_cav
end


# The user can just pass a single float 'p0' and the vector of initial conditions 
# is assumed to be homogeneous
function init_p_cav(graph::HGraph, p0::Float64)
    p0 = fill(p0, graph.N)
    init_p_cav(graph, p0)
end


function sum_product_KSAT()
    
end