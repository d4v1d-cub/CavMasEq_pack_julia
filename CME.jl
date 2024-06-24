# This function initializes the cavity conditional probabilities with a fully 
# independent initial distribution given by the vector 'p0'.
# It produces an array p_cav[he, i, val_cond, chain_others]
function init_p_cav(graph::HGraph, p0::Vector{Float64})
    ch_exc = graph.chains_he ÷ 2
    p_cav = zeros(Float64, (graph.M, graph.K, 2, ch_exc))
    for he in 1:graph.M
        for i in 1:graph.K
            for s in 1:2
                for ch in 0:ch_exc - 1
                    bits = digits(ch, base=2, pad=graph.K-1)
                    p_cav[he, i, s, ch + 1] = prod(bits + (1 .- 2 * bits) .* p0[graph.nodes_except[he, i, :]])
                end
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


# This function computes one term in the derivative of pcav. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_cav. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch_exc is an integer that codes the combination of the variables in the argument of the
# probability 'p_cav'. In this case all p_cav are conditioned to one already unsatisfied link.
function der_pcav_contr_node_unsat(p_cav::Vector{Float64}, all_sums::Array{Float64, 3}, 
    ch_exc::Int64, ch_exc_unsat::Int64, 
    index_in::Int64)

    ch_exc_flip = (ch_exc ⊻ (2 ^ (index_in - 1)))  # The ⊻ (xor) operation flips the variable
    val = ((ch_exc >> (index_in - 1)) & 1)         # Takes the value of the variable
    Ea = (ch_exc == ch_exc_unsat)                  # As the variable in the conditional is set 
    Ea_flip = (ch_exc_flip == ch_exc_unsat)        # to always unsatisfy its link, the clause
    # will be unsatisfied only if the rest of the variables are in their unsatisfied configuration

    return -all_sums[val + 1, Ea + 1, Ea_flip + 1] * p_cav[ch_exc + 1] + 
        all_sums[2 - val, Ea_flip + 1, Ea + 1] * p_cav[ch_exc_flip + 1]
end



# This function computes one term in the derivative of pcav. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_cav. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch_exc is an integer that codes the combination of the variables in the argument of the
# probability 'p_cav'. In this case all p_cav are conditioned to one satisfied link.
function der_pcav_contr_node_sat(p_cav::Vector{Float64}, all_sums::Array{Float64, 3}, 
    ch_exc::Int64, index_in::Int64)

    ch_exc_flip = (ch_exc ⊻ (2 ^ (index_in - 1)))  # The ⊻ (xor) operation flips the variable
    val = ((ch_exc >> (index_in - 1)) & 1)         # Takes the value of the variable

    return -all_sums[val + 1, 1, 1] * p_cav[ch_exc + 1] + 
    all_sums[2 - val, 1, 1] * p_cav[ch_exc_flip + 1]
end



# This function computes the part of the derivatives of the p_cav that correspond to
# the flipping of variable 'node' inside a clause 'he' conditioned on another node that is
# defined in the function all_ders_node_KSAT (see below)
# The result is cumulated in the matrix 'ders'
function der_pcav_KSAT(p_cav::Matrix{Float64}, pu::Array{Float64, 3}, he::Int64, 
    node::Int64, place_node::Int64, ch_exc_unsat::Int64, graph::HGraph, 
    all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ders::Matrix{Float64}, nch::Int64)

    all_sums = compute_all_sums(pu, node, he, graph, all_lp, all_lm, ratefunc, rate_args)
    for ch_exc in 0:nch-1
        ders[1, ch_exc + 1] += der_pcav_contr_node_sat(p_cav[1, :], all_sums, ch_exc, place_node)
        ders[2, ch_exc + 1] += der_pcav_contr_node_unsat(p_cav[2, :], all_sums, ch_exc, ch_exc_unsat, 
                                                         place_node)   
    end
    return all_sums, ders
    # It returns the sums corresponding to flipping node with fixed values of 'he' and the 
    # computed derivatives
end


# Given the sum-products corresponding to the flipping variable 'i', summed over all the clauses
# containing 'i' different from a given 'he', this function computes the derivative of the local
# probability 'probi'. It receives the cavity conditional probability of having 'he' unsatisfied
# (pu_he) and the value of the link between 'i' and 'he' (li_he)
function der_pi_KSAT(probi::Float64, pu_he::Vector{Float64}, all_sums::Array{Float64, 3}, 
             li_he::Int8)
    cumul = 0.0
    for u_neigh in 0:1
        pr = 1 - u_neigh - (1 - 2 * u_neigh) * pu_he[li_he + 1]
        pr_flip = 1 - u_neigh - (1 - 2 * u_neigh) * pu_he[2 - li_he]
        Ea = u_neigh * (li_he == 1)             # 'u_neigh' controls the state of the rest of the variables in 'he'
        # If u_neigh = 1, the rest of the variables unsatisfy their links and u=0 otherwise 
        # However, the clause 'he' is unsat only is 'i' also unsatisfies its link
        # As 'probi' is defined for si = 0, this happens only if 'li_he=1'
        Ea_flip = u_neigh * (li_he == 0)       # On the other hand, if li_he=0 the clause will be unsatisfied
        # after flipping 'si'
        cumul += -all_sums[1, Ea + 1, Ea_flip + 1] * pr * probi + 
        all_sums[2, Ea_flip + 1, Ea + 1] * pr_flip * (1 - probi)
    end

    return cumul
end


# This function computes the all the derivatives related to a hyperedge 'he'
function all_ders_he_KSAT(p_cav::Array{Float64, 4}, pu::Array{Float64, 3}, he::Int64, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ders_pcav::Array{Float64, 4}, nch::Int64, 
    ch_u_cond::Matrix{Int64}, all_sums::Array{Float64, 4})

    inner_sums = zeros(Float64, (2, 2, 2))
    for i in 1:graph.K                      # Computing all the derivatives of cavity probabilities
        node = graph.he_2_var[he, i]        # in the hyperedge 'he' when 'node' is flipped
        for j in 1:graph.K - 1        
            node_neigh = graph.nodes_except[he, i, j]
            place_neigh = graph.nodes_in[he][node_neigh] 
            ch_exc_unsat = ch_u_cond[he, place_neigh]  # gives the configuration that unsatisfies every link
                                                    # converting it to an integer
            place_in_exc = graph.place_there[he, place_neigh][node] # Gets the index of 'node' in the 
                                                    # array  graph.nodes_except[he, place_neigh, :]
            inner_sums, ders_pcav[he, place_neigh, :, :] = 
            der_pcav_KSAT(p_cav[he, place_neigh, :, :], pu, he, node, place_in_exc, ch_exc_unsat, 
            graph, all_lp, all_lm, ratefunc, rate_args, ders_pcav[he, place_neigh, :, :], nch) 

            if he == graph.var_2_he[node][end]
                all_sums[node, :, :, :] .= inner_sums
            end
        end
    end
end


function all_ders_node_KSAT(probi::Float64, pu::Array{Float64, 3}, node::Int64, graph::HGraph, 
    links::Matrix{Int8}, all_sums::Array{Float64, 4})

    if length(graph.var_2_he[node]) > 0
        he = graph.var_2_he[node][end]
        place_in = graph.nodes_in[he][node]
        return der_pi_KSAT(probi, pu[he, place_in, :], all_sums[node, :, :, :], links[he, place_in])
    else
        return 0
    end
end


# This function computes all the derivatives for the CME in the K-SAT
function all_ders_CME_KSAT(p_cav::Array{Float64, 4}, probi::Vector{Float64}, pu::Array{Float64, 3}, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, links::Matrix{Int8}, ch_u_cond::Matrix{Int64})

    d_pcav = zeros(Float64, size(p_cav))
    d_node = zeros(Float64, size(probi))
    all_sums = zeros(Float64, (graph.N, 2, 2, 2))
    nch = graph.chains_he ÷ 2

    Threads.@threads for he in 1:graph.M
    all_ders_he_KSAT(p_cav, pu, he, graph, all_lp, all_lm, ratefunc, rate_args, d_pcav, 
                    nch, ch_u_cond, all_sums)
    end

    # for he in 1:graph.M
    #     all_ders_he_KSAT(p_cav, pu, he, graph, all_lp, all_lm, ratefunc, rate_args, d_pcav, 
    #                     nch, ch_u_cond, all_sums)
    # end

    for node in 1:graph.N
    d_node[node] = all_ders_node_KSAT(probi[node], pu, node, graph, links, all_sums)
    end

    GC.gc()
    return d_pcav, d_node
end


# This function computes the derivative of the probabilities and feeds the julia's ODE integrator
function fder_KSAT_CME(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build, len_cav, efinal = p
    # These are the parameters of the integration
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.chains_he ÷ 2))
    probi = u[len_cav + 1:len_cav + graph.N]
    # The probabilities are reshaped in their original forms

    pu_cond = comp_pu_KSAT(p_cav, graph, ch_u_cond)
    p_joint_u = get_pju_CME(graph, probi, pu_cond, ch_u)

    st = State_CME(p_cav, probi, pu_cond, p_joint_u)  
    # A struct of type state is created just to pass it as an argument
    # for the builder of the rates' function arguments
    # This way, the user can choose what information to use inside the rate

    rates_arg = rarg_build(graph, st, rarg_cst...)

    d_pc, d_pi = all_ders_CME_KSAT(p_cav, probi, pu_cond, graph, all_lp, all_lm, rfunc, rates_arg, links, 
                                   ch_u_cond)

    GC.gc()
    du .= vcat(reshape(d_pc, len_cav), d_pi)
end


function get_pju_CME(graph::HGraph, probi::Vector{Float64}, pu_cond::Array{Float64, 3}, 
                     ch_u::Vector{Int64})
    p_joint_u = zeros(Float64, graph.M)
    for he in 1:graph.M
        node = graph.he_2_var[he, 1]
        bit_node = ch_u[he] & 1
        p_joint_u[he] = (bit_node + (1 - 2 * bit_node) * probi[node]) * pu_cond[he, 1, 2]
    end
    return p_joint_u
end


function reshape_u_to_probs_CME(u::Vector{Float64}, integrator)
    graph = get_graph(integrator)
    len_cav = integrator.p[10]
    p_cav = reshape(u[1:len_cav], (graph.M, graph.K, 2, graph.chains_he ÷ 2))
    probi = u[len_cav + 1:len_cav + graph.N]
    return p_cav, probi
end


function get_pju_CME(u::Vector{Float64}, integrator)
    graph = get_graph(integrator)
    ch_u = get_ch_u(integrator)
    ch_u_cond = get_ch_u_cond(integrator)
    p_cav, probi = reshape_u_to_probs_CME(u, integrator)
    pu_cond = comp_pu_KSAT(p_cav, graph, ch_u_cond)
    return get_pju_CME(graph, probi, pu_cond, ch_u)
end

function save_ener_CME(u, t, integrator)
    p_joint_u = get_pju_CME(u, integrator)
    e = ener(p_joint_u)
    println(t, "\t", e)
    return e
end


function stopcond_CME(u, t, integrator)
    efinal = get_efinal(integrator)
    p_joint_u = get_pju_CME(u, integrator)
    e = ener(p_joint_u)
    return e - efinal
end


# This function integrates the CME's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CME_KSAT_base(ratefunc::Function, rargs_cst, rarg_build::Function, 
                       graph::HGraph, links::Matrix{Int8}, tspan::Vector{Float64}, p0::Float64, 
                       method, eth::Float64, cbs_save::CallbackSet, dt_s::Float64, 
                       abstol::Float64, reltol::Float64)
    efinal = eth * graph.N
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_cav = init_p_cav(graph, p0)
    probi = fill(p0, graph.N)
    len_cav = length(p_cav)
    params = graph, all_lp, all_lm, links, ch_u, ch_u_cond, ratefunc, rargs_cst, rarg_build, len_cav, 
             efinal
    u0 = vcat(reshape(p_cav, len_cav), probi)
    prob = ODEProblem(fder_KSAT_CME, u0, tspan, params)

    affect!(integrator) = terminate!(integrator)
    cb_stop = ContinuousCallback(stopcond_CME, affect!)

    cbs = CallbackSet(cbs_save, cb_stop)

    sol = solve(prob, method, progress=true, callback=cbs, saveat=dt_s, 
                abstol=abstol, reltol=reltol)
    return sol
end


function fder_CME_custom(p_cav::Array{Float64, 4}, probi::Vector{Float64}, pu_cond::Array{Float64, 3}, 
    p_joint_u::Vector{Float64}, graph::HGraph, ch_u_cond::Matrix{Int64}, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    rfunc::Function, rarg_cst, rarg_build::Function, links::Matrix{Int8})

    st = State_CME(p_cav, probi, pu_cond, p_joint_u)  
    # A struct of type state is created just to pass it as an argument
    # for the builder of the rates' function arguments
    # This way, the user can choose what information to use inside the rate

    rates_arg = rarg_build(graph, st, rarg_cst...)

    return all_ders_CME_KSAT(p_cav, probi, pu_cond, graph, all_lp, all_lm, rfunc, rates_arg, links, 
                             ch_u_cond)
end


# Euler step to intialize Adams Moulton method of the 2nd order
function init_Euler_CME(p_cav::Array{Float64, 4}, probi::Vector{Float64}, pu_cond::Array{Float64, 3}, 
    p_joint_u::Vector{Float64}, graph::HGraph, ch_u_cond::Matrix{Int64}, 
    all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    rfunc::Function, rarg_cst, rarg_build::Function, links::Matrix{Int8}, dt0::Float64)


    d_pc, d_pi = fder_CME_custom(p_cav, probi, pu_cond, p_joint_u, graph, ch_u_cond, all_lp, all_lm, 
                             rfunc, rarg_cst, rarg_build, links)

    p_cav_1 = p_cav .+ dt0 * d_pc
    probi_1 = probi .+ dt0 * d_pi

    return p_cav_1, probi_1, d_pc, d_pi
end


# Adams Moulton predictor-corrector method of the second order
function AM2_CME(p_cav_0::Array{Float64, 4}, probi_0::Vector{Float64}, graph::HGraph, 
    ch_u::Vector{Int64}, ch_u_cond::Matrix{Int64}, all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, rfunc::Function, rarg_cst, rarg_build::Function, 
    links::Matrix{Int8}, t0::Float64, dt0::Float64, tl::Float64, tol::Float64, efinal::Float64, 
    dt_min::Float64, max_iter::Int64)

    t_list = Vector{Float64}()
    e_list = Vector{Float64}()
    
    pu_cond = comp_pu_KSAT(p_cav_0, graph, ch_u_cond)
    p_joint_u = get_pju_CME(graph, probi_0, pu_cond, ch_u)

    e = ener(p_joint_u)
    push!(t_list, t0)
    push!(e_list, e)

    println(t0, "\t", e)

    p_cav_1, probi_1, d_pc_0, d_pi_0 = init_Euler_CME(p_cav_0, probi_0, pu_cond, p_joint_u, graph, ch_u_cond, 
                                                  all_lp, all_lm, rfunc, rarg_cst, rarg_build, links, dt0)

    pu_cond .= comp_pu_KSAT(p_cav_1, graph, ch_u_cond)
    p_joint_u .= get_pju_CME(graph, probi_1, pu_cond, ch_u)
    t = t0 + dt0

    e = ener(p_joint_u)
    push!(t_list, t)
    push!(e_list, e)

    println(t, "\t", e)

    p_cav_2 = zeros(Float64, size(p_cav_1))
    probi_2 = zeros(Float64, size(probi_1))
    p_cav_pred = zeros(Float64, size(p_cav_1))
    probi_pred = zeros(Float64, size(probi_1))
    d_pc_1 = zeros(Float64, size(d_pc_0))
    d_pi_1 = zeros(Float64, size(d_pi_0))
    d_pc_pred = zeros(Float64, size(d_pc_0))
    d_pi_pred = zeros(Float64, size(d_pi_0))
    dif_pc = zeros(Float64, size(d_pc_0))
    dif_pi = zeros(Float64, size(d_pi_0))

    dt1 = dt0

    while t < tl
        if e < efinal
            break
        end

        d_pc_1, d_pi_1 = fder_CME_custom(p_cav_1, probi_1, pu_cond, p_joint_u, graph, ch_u_cond, all_lp, 
                                     all_lm, rfunc, rarg_cst, rarg_build, links)

        cond = true
        counter = 0
        while cond && counter < max_iter
            dif_pc .= dt1 / dt0 * (d_pc_1 .- d_pc_0)
            dif_pi .= dt1 / dt0 * (d_pi_1 .- d_pi_0)

            p_cav_pred .= p_cav_1 .+ dt1 * (d_pc_1 .+ 0.5 * dif_pc)
            probi_pred .= probi_1 .+ dt1 * (d_pi_1 .+ 0.5 * dif_pi) 

            pu_cond .= comp_pu_KSAT(p_cav_pred, graph, ch_u_cond)
            p_joint_u .= get_pju_CME(graph, probi_pred, pu_cond, ch_u)

            d_pc_pred, d_pi_pred = fder_CME_custom(p_cav_pred, probi_pred, pu_cond, p_joint_u, graph, 
                                               ch_u_cond, all_lp, all_lm, rfunc, rarg_cst, rarg_build, 
                                               links)

            p_cav_2 .= p_cav_pred .+ dt1 * (0.5 - dt1 / 6 / (dt1 + dt0)) * (d_pc_pred .- d_pc_1 .- dif_pc)
            probi_2 .= probi_pred .+ dt1 * (0.5 - dt1 / 6 / (dt1 + dt0)) * (d_pi_pred .- d_pi_1 .- dif_pi)
            
            err_cav = maximum(abs.(p_cav_2 .- p_cav_pred) ./ p_cav_2)
            err_pi = maximum(abs.(probi_2 .- probi_pred) ./ probi_2)
            err = max(err_cav, err_pi)

            if minimum(p_cav_2) < 0 || minimum(probi_2) < 0
                dt1 /= 2
                if dt1 < dt_min
                    dt_min /= 2
                end
            elseif tol * 0.9 < err < tol * 1.001
                cond = false
                dt1 = max(dt1 * sqrt(tol / err), dt_min)
            else
                dt1 = max(dt1 * sqrt(tol / err), dt_min)
            end

            counter += 1
            GC.gc()
        end

        if counter == max_iter
            println("Maximum number of iterations excedeed inside ABM_CME")
        end

        p_cav_0 .= p_cav_1
        probi_0 .= probi_1

        p_cav_1 .= p_cav_2
        probi_1 .= probi_2

        pu_cond .= comp_pu_KSAT(p_cav_1, graph, ch_u_cond)
        p_joint_u .= get_pju_CME(graph, probi_1, pu_cond, ch_u)
        t += dt1
    
        e = ener(p_joint_u)
        push!(t_list, t)
        push!(e_list, e)
    
        println(t, "\t", e)

        dt0 = dt1
    end

    return t_list, e_list
end


# ODE integration by hand
function CME_KSAT_custom(ratefunc::Function, rargs_cst, rarg_build::Function, 
    graph::HGraph, links::Matrix{Int8}, tspan::Vector{Float64}, p0::Float64, 
    method::Function, eth::Float64, tol::Float64, dt0::Float64, dt_min::Float64, 
    max_iter::Int64)


    efinal = eth * graph.N
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_cav = init_p_cav(graph, p0)
    probi = fill(p0, graph.N)

    t0 = tspan[1]
    tl = tspan[2]


    t_list, e_list = method(p_cav, probi, graph, ch_u, ch_u_cond, all_lp, all_lm, ratefunc, rargs_cst, 
                            rarg_build, links, t0, dt0, tl, tol, efinal, dt_min, max_iter)
    return t_list, e_list
end


# Decorator for the function integrates the CME's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CME_KSAT(ratefunc::Function, rargs_cst, rarg_build::Function;
                  graph::HGraph=build_empty_graph(), 
                  N::Int64=0, K::Int64=0, alpha::Union{Float64, Int64}=0.0, seed_g::Int64=rand(1:typemax(Int64)),
                  links::Matrix{Int8}=Matrix{Int8}(undef, 0, 0), seed_l::Int64=rand(1:typemax(Int64)), 
                  tspan::Vector{Float64}=[0.0, 1.0], 
                  p0::Float64=0.5, method=VCABM(), eth::Float64=1e-6, cbs_save::CallbackSet=CallbackSet(),
                  dt_s::Float64=0.1, abstol::Float64=1e-6, reltol::Float64=1e-3, custom=false,
                  dt0::Float64=0.01, dt_min::Float64=1e-7, max_iter=100)
    if N > 0
        if graph.N == 0
            c = K * alpha
            graph = build_ER_HGraph(N, c, K, seed_g)
        end

        if length(links) == 0
            links = gen_links(graph, seed_l)
        end

        if custom
            return CME_KSAT_custom(ratefunc, rargs_cst, rarg_build, graph, links, tspan, p0, method, eth, 
                                   reltol, dt0, dt_min, max_iter)
        else
            return CME_KSAT_base(ratefunc, rargs_cst, rarg_build, graph, links, tspan, p0,
            method, eth, cbs_save, dt_s, abstol, reltol)
        end

        
    else
        throw("In CME_KSAT function: The user should provide either a graph::HGraph or valid values for N, K and alpha")
    end  
end

