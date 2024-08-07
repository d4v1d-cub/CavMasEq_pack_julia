# This function initializes the joint probabilities in a hyperedge with a fully 
# independent initial distribution given by the vector 'p0'.
# It produces an array p_joint[he, chain]
function init_p_joint(graph::HGraph, p0::Vector{Float64})
    p_joint = zeros(Float64, (graph.M, graph.chains_he))
    for he in 1:graph.M
        for ch in 0:graph.chains_he - 1
            bits = digits(ch, base=2, pad=graph.K)
            p_joint[he, ch + 1] = prod(bits + (1 .- 2 * bits) .* p0[graph.he_2_var[he, :]])
        end
    end
    return p_joint
end


# The user can just pass a single float 'p0' and the vector of initial conditions 
# is assumed to be homogeneous
function init_p_joint(graph::HGraph, p0::Float64)
    p0 = fill(p0, graph.N)
    init_p_joint(graph, p0)
end


function compute_p_cond(p_joint::Matrix{Float64}, graph::HGraph, links::Matrix{Int8})
    nch_exc = graph.chains_he ÷ 2
    p_cond = Array{Float64, 4}(undef, graph.M, graph.K, 2, nch_exc)
    for he in 1:graph.M
        for i in 1:graph.K
            li = links[he, i]
            for s in 1:2
                sat_i = (s - 1) == li
                for ch_exc in 0:nch_exc - 1
                    ch = convert_chain(ch_exc, s, i)
                    p_cond[he, i, 2 - sat_i, ch_exc + 1] = p_joint[he, ch + 1] 
                end                
                p_cond[he, i, 2 - sat_i, :] .= normalize(p_cond[he, i, 2 - sat_i, :])
            end
        end
    end
    return p_cond
end


# This function computes one term in the derivative of p_joint. It corresponds to
# the flipping of the variable indexed as 'index_in' among the variables in the argument
# of p_joint. The array 'all_sums' is computed using the function 'compute_all_sums'
# ch is an integer that codes the combination of the variables in the argument of the
# probability 'p_joint'.
function der_pjoint_contr_node(p_joint::Vector{Float64}, all_sums::Array{Float64, 3}, 
    ch::Int64, ch_unsat::Int64, index_in::Int64)

    ch_flip = (ch ⊻ (2 ^ (index_in - 1)))      # The ⊻ (xor) operation flips the variable
    val = ((ch >> (index_in - 1)) & 1)         # Takes the value of the variable
    Ea = (ch == ch_unsat)                      # Compares with the unsatisfied combination 
    Ea_flip = (ch_flip == ch_unsat)            # Compares the state with the variable flipped with
    # the unsat combination inside the clause

    return -all_sums[val + 1, Ea + 1, Ea_flip + 1] * p_joint[ch + 1] + 
        all_sums[2 - val, Ea_flip + 1, Ea + 1] * p_joint[ch_flip + 1]
end



# This function computes the part of the derivatives of the p_cav that correspond to
# the flipping of variable 'node' inside a clause 'he' conditioned on another node that is
# defined in the function all_ders_node_KSAT (see below)
# The result is cumulated in the matrix 'ders'
function der_pjoint_KSAT(p_joint::Vector{Float64}, pu::Array{Float64, 3}, he::Int64, 
    node::Int64, place_node::Int64, ch_unsat::Int64, graph::HGraph, 
    all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ders::Vector{Float64})

    all_sums = compute_all_sums(pu, node, he, graph, all_lp, all_lm, ratefunc, rate_args)
    for ch in 0:graph.chains_he-1
        ders[ch + 1] += der_pjoint_contr_node(p_joint, all_sums, ch, ch_unsat, place_node)   
    end
    return ders
    # It returns the sums corresponding to flipping node with fixed values of 'he' and the 
    # computed derivatives
end


# This function computes the all the derivatives related to a hyperedge 'he'
function all_ders_he_KSAT(p_joint::Matrix{Float64}, pu::Array{Float64, 3}, he::Int64, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ders::Matrix{Float64}, ch_u::Vector{Int64})

    for i in 1:graph.K                      # Computing all the derivatives of cavity probabilities
        node = graph.he_2_var[he, i]        # in the hyperedge 'he' when 'node' is flipped
        ders[he, :] .= der_pjoint_KSAT(p_joint[he, :], pu, he, node, i, ch_u[he], 
                                       graph, all_lp, all_lm, ratefunc, rate_args, ders[he, :]) 
    end 
end


# This function computes all the derivatives for the CDA in the K-SAT
function all_ders_CDA_KSAT(p_joint::Matrix{Float64}, pu::Array{Float64, 3}, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    ratefunc::Function, rate_args, ch_u::Vector{Int64})

    ders = zeros(Float64, size(p_joint))

    Threads.@threads for he in 1:graph.M
    all_ders_he_KSAT(p_joint, pu, he, graph, all_lp, all_lm, ratefunc, rate_args, ders, 
                     ch_u)
    end

    # for he in 1:graph.M
    #     all_ders_he_KSAT(p_joint, pu, he, graph, all_lp, all_lm, ratefunc, rate_args, ders, 
    #                      ch_u)
    # end

    GC.gc()
    return ders
end


# This function computes the derivative of the probabilities and feeds the julia's ODE integrator
function fder_KSAT_CDA(du::Vector{Float64}, u::Vector{Float64}, p, t::Float64)
    graph, all_lp, all_lm, links, ch_u, ch_u_cond, rfunc, rarg_cst, rarg_build, efinal = p
    # These are the parameters of the integration
    p_joint = reshape(u, (graph.M, graph.chains_he))
    # The probabilities are reshaped in their original forms

    p_cond = compute_p_cond(p_joint, graph, links)
    pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
    p_joint_u = get_pju_CDA(graph, p_joint, ch_u)

    st = State_CDA(p_joint, pu_cond, p_joint_u)  
    # A struct of type state is created just to pass it as an argument
    # for the builder of the rates' function arguments
    # This way, the user can choose what information to use inside the rate

    rates_arg = rarg_build(graph, st, rarg_cst...)

    dp = all_ders_CDA_KSAT(p_joint, pu_cond, graph, all_lp, all_lm, rfunc, rates_arg, ch_u)

    GC.gc()
    du .= reshape(dp, length(du))
end


function get_pju_CDA(graph::HGraph, p_joint::Matrix{Float64}, ch_u::Vector{Int64})
    p_joint_u = zeros(Float64, graph.M)
    for he in 1:graph.M
        p_joint_u[he] = p_joint[he, ch_u[he] + 1]
    end
    return p_joint_u
end


function reshape_u_to_probs_CDA(u::Vector{Float64}, integrator)
    graph = get_graph(integrator)
    p_joint = reshape(u, (graph.M, graph.chains_he))
    return p_joint
end


function get_pju_CDA(u::Vector{Float64}, integrator)
    graph = get_graph(integrator)
    ch_u = get_ch_u(integrator)
    p_joint = reshape_u_to_probs_CDA(u, integrator)
    return get_pju_CDA(graph, p_joint, ch_u)
end


function save_ener_CDA(u, t, integrator)
    p_joint_u = get_pju_CDA(u, integrator)
    e = ener(p_joint_u)
    println(t, "\t", e)
    return e
end

function stopcond_CDA(u, t, integrator)
    efinal = get_efinal(integrator)
    p_joint_u = get_pju_CDA(u, integrator)
    e = ener(p_joint_u)
    return e - efinal
end

# This function integrates the CDA's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CDA_KSAT_base(ratefunc::Function, rargs_cst, rarg_build::Function, 
                       graph::HGraph, links::Matrix{Int8}, tspan::Vector{Float64}, p0::Float64, 
                       method, eth::Float64, cbs_save::CallbackSet, dt_s::Float64, 
                       abstol::Float64, reltol::Float64)
    efinal = eth * graph.N
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_joint = init_p_joint(graph, p0)
    params = graph, all_lp, all_lm, links, ch_u, ch_u_cond, ratefunc, rargs_cst, rarg_build, efinal
    u0 = reshape(p_joint, graph.M * graph.chains_he)
    prob = ODEProblem(fder_KSAT_CDA, u0, tspan, params)

    affect!(integrator) = terminate!(integrator)
    cb_stop = ContinuousCallback(stopcond_CDA, affect!)

    cbs = CallbackSet(cbs_save, cb_stop)

    sol = solve(prob, method, progress=true, callback=cbs, saveat=dt_s, 
                abstol=abstol, reltol=reltol)
    return sol
end


function fder_CDA_custom(p_joint::Matrix{Float64}, pu_cond::Array{Float64, 3}, p_joint_u::Vector{Float64}, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    rfunc::Function, rarg_cst, rarg_build::Function, ch_u::Vector{Int64})
    

    st = State_CDA(p_joint, pu_cond, p_joint_u)  
    # A struct of type state is created just to pass it as an argument
    # for the builder of the rates' function arguments
    # This way, the user can choose what information to use inside the rate

    rates_arg = rarg_build(graph, st, rarg_cst...)

    dp = all_ders_CDA_KSAT(p_joint, pu_cond, graph, all_lp, all_lm, rfunc, rates_arg, ch_u)

    return dp
end


function init_Euler_CDA(p_joint::Matrix{Float64}, pu_cond::Array{Float64, 3}, p_joint_u::Vector{Float64}, 
    graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, all_lm::Vector{Vector{Vector{Int64}}}, 
    rfunc::Function, rarg_cst, rarg_build::Function, ch_u::Vector{Int64}, dt0::Float64)

    dp = fder_CDA_custom(p_joint, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, rarg_cst, rarg_build,
                         ch_u)
    p_joint_1 = p_joint .+ dt0 * dp

    return p_joint_1, dp
end



# Adams Moulton predictor-corrector method of the second order
function AM2_CDA(p_joint_0::Matrix{Float64}, graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, rfunc::Function, rarg_cst, rarg_build::Function, 
    ch_u::Vector{Int64}, ch_u_cond::Matrix{Int64}, links::Matrix{Int8}, t0::Float64, dt0::Float64, tl::Float64, 
    efinal::Float64, fileener::String, tol::Float64=1e-2, dt_min::Float64=1e-7, max_iter::Int64=100)

    fe = open(fileener, "w")
    
    t_list = Vector{Float64}()
    e_list = Vector{Float64}()
    
    p_cond = compute_p_cond(p_joint_0, graph, links)
    pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
    p_joint_u = get_pju_CDA(graph, p_joint_0, ch_u)

    e = ener(p_joint_u)
    push!(t_list, t0)
    push!(e_list, e)

    println("t=", t0)
    write(fe, string(t0) * "\t" * string(e) * "\n")
    flush(fe)

    p_joint_1, dp_0 = init_Euler_CDA(p_joint_0, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, 
                                   rarg_cst, rarg_build, ch_u, dt0)

    p_cond = compute_p_cond(p_joint_1, graph, links)
    pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
    p_joint_u = get_pju_CDA(graph, p_joint_1, ch_u)
    t = t0 + dt0

    e = ener(p_joint_u)
    push!(t_list, t)
    push!(e_list, e)

    println("t=", t)
    write(fe, string(t) * "\t" * string(e) * "\n")
    flush(fe)

    p_joint_2 = zeros(Float64, size(p_joint_1))

    dt1 = dt0

    while t < tl
        if e < efinal
            println("Final energy reached")
            break
        end

        p_cond .= compute_p_cond(p_joint_1, graph, links)
        pu_cond .= comp_pu_KSAT(p_cond, graph, ch_u_cond)
        p_joint_u .= get_pju_CDA(graph, p_joint_1, ch_u)
        dp_1 = fder_CDA_custom(p_joint_1, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, rarg_cst, rarg_build,
                               ch_u)

        counter = 0
        cond = true
        while cond && counter < max_iter
            dif_p = dt1 / dt0 * (dp_1 .- dp_0)

            p_pred = p_joint_1 .+ dt1 * (dp_1 .+ 0.5 * dif_p) 

            p_cond = compute_p_cond(p_pred, graph, links)
            pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
            p_joint_u = get_pju_CDA(graph, p_pred, ch_u)

            dp_pred = fder_CDA_custom(p_pred, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, rarg_cst, rarg_build,
                                      ch_u)

            p_joint_2 = p_pred .+ dt1 * (0.5 - dt1 / 6 / (dt1 + dt0)) * (dp_pred .- dp_1 .- dif_p)
            
            err = maximum(abs.(p_joint_2 .- p_pred) ./ p_joint_2)

            if minimum(p_joint_2) < 0
                println("Some probabilities are negative")
                dt1 /= 2
                println("step divided by half")
                if dt1 < dt_min
                    dt_min /= 2
                    println("dt_min also halfed")
                end
            elseif err < tol * 1.01
                cond = false
                println("step dt=", dt1, "  accepted")
                dt1 = dt1 * sqrt(tol / err)
            else
                println("step dt=", dt1, "  rejected")
                println("iter=", counter)
                dt1 = dt1 * sqrt(tol / err)
            end

            counter += 1
            GC.gc()
        end

        if counter == max_iter
            println("Maximum number of iterations excedeed inside ABM_CME")
        end

        p_joint_0 .= p_joint_1

        p_joint_1 .= p_joint_2

        p_cond = compute_p_cond(p_joint_1, graph, links)
        pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
        p_joint_u = get_pju_CDA(graph, p_joint_1, ch_u)
        t += dt1
    
        e = ener(p_joint_u)
        push!(t_list, t)
        push!(e_list, e)
    
        println("t=", t)
        write(fe, string(t) * "\t" * string(e) * "\n")
        flush(fe)

        dt0 = dt1
    end

    close(fe)
    return t_list, e_list
end


# Adams Moulton predictor-corrector method of the second order
function RK2_CDA(p_joint::Matrix{Float64}, graph::HGraph, all_lp::Vector{Vector{Vector{Int64}}}, 
    all_lm::Vector{Vector{Vector{Int64}}}, rfunc::Function, rarg_cst, rarg_build::Function, 
    ch_u::Vector{Int64}, ch_u_cond::Matrix{Int64}, links::Matrix{Int8}, t0::Float64, dt0::Float64, tl::Float64, 
    efinal::Float64, fileener::String; tol::Float64=1e-2, dt_min::Float64=1e-7)

    fe = open(fileener, "w")
    
    t_list = Vector{Float64}()
    e_list = Vector{Float64}()
    
    p_cond = compute_p_cond(p_joint, graph, links)
    pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
    p_joint_u = get_pju_CDA(graph, p_joint, ch_u)

    e = ener(p_joint_u)
    push!(t_list, t0)
    push!(e_list, e)

    println("t=", t0)
    write(fe, string(t0) * "\t" * string(e) * "\n")
    flush(fe)

    k1 = zeros(Float64, size(p_joint))
    k2 = zeros(Float64, size(p_joint))
    p_joint_1 = zeros(Float64, size(p_joint))
    dp = zeros(Float64, size(p_joint))

    dt1 = dt0
    t = t0

    while t < tl
        if e < efinal
            println("Final energy reached")
            break
        end

        p_cond .= compute_p_cond(p_joint, graph, links)
        pu_cond .= comp_pu_KSAT(p_cond, graph, ch_u_cond)
        p_joint_u .= get_pju_CDA(graph, p_joint, ch_u)

        dp .= fder_CDA_custom(p_joint, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, rarg_cst, rarg_build,
                              ch_u)

        k1 .= dt1 * dp
        p_joint_1 .= (p_joint .+ k1)

        p_cond .= compute_p_cond(p_joint_1, graph, links)
        pu_cond .= comp_pu_KSAT(p_cond, graph, ch_u_cond)
        p_joint_u .= get_pju_CDA(graph, p_joint_1, ch_u)

        dp .= fder_CDA_custom(p_joint_1, pu_cond, p_joint_u, graph, all_lp, all_lm, rfunc, rarg_cst, rarg_build,
                              ch_u)

        k2 .= dt1 * dp

        valid = minimum(p_joint .+ (k1 .+ k2) / 2) >= 0

        if !valid
            println("Some probabilities would be negative if dt=", dt1, " is taken")
            dt1 /= 2
            println("step divided by half")
            if dt1 < dt_min
                dt_min /= 2
                println("dt_min also halfed")
            end
        else
            error = maximum(abs.(k1 .- k2))
            if error < 2 * tol
                println("step dt=", dt1, "  accepted")
                t += dt1
                p_joint .= (p_joint .+ (k1 .+ k2) / 2)

                p_cond = compute_p_cond(p_joint, graph, links)
                pu_cond = comp_pu_KSAT(p_cond, graph, ch_u_cond)
                p_joint_u = get_pju_CDA(graph, p_joint, ch_u)

                e = ener(p_joint_u)
                push!(t_list, t)
                push!(e_list, e)

                println("t=", t)
                write(fe, string(t) * "\t" * string(e) * "\n")
                flush(fe)
                
            else
                println("step dt=", dt1, "  rejected  new step will be attempted")
                println("error=", error)
            end

            dt1 = 4 * dt1 * sqrt(2 * tol / error) / 5

            if dt1 > graph.M
                dt1 = M
            elseif dt1 < dt_min
                dt1 = dt_min
            end
            
            println("Recommended step is dt=", dt1)
        end
        GC.gc()
    end

    close(fe)
    return t_list, e_list
end


# ODE integration by hand
function CDA_KSAT_custom(ratefunc::Function, rargs_cst, rarg_build::Function, 
    graph::HGraph, links::Matrix{Int8}, tspan::Vector{Float64}, p0::Float64, 
    method::Function, eth::Float64, tol::Float64, dt0::Float64, dt_min::Float64, 
    fileener::String)


    efinal = eth * graph.N
    all_lp, all_lm = all_lpm(graph, links)
    ch_u, ch_u_cond = unsat_ch(graph, links)
    p_joint = init_p_joint(graph, p0)

    t0 = tspan[1]
    tl = tspan[2]


    t_list, e_list = method(p_joint, graph, all_lp, all_lm, ratefunc, rargs_cst, rarg_build, 
                            ch_u, ch_u_cond, links, t0, dt0, tl, efinal, fileener, tol=tol, 
                            dt_min=dt_min)
    return t_list, e_list
end


# Decorator for the function integrates the CDA's equations for a specific algorithm (given by ratefunc)
# and some boolean formula (given by graph and links)
function CDA_KSAT(ratefunc::Function, rargs_cst, rarg_build::Function;
                  graph::HGraph=build_empty_graph(), 
                  N::Int64=0, K::Int64=0, alpha::Union{Float64, Int64}=0.0, seed_g::Int64=rand(1:typemax(Int64)),
                  links::Matrix{Int8}=Matrix{Int8}(undef, 0, 0), seed_l::Int64=rand(1:typemax(Int64)), 
                  tspan::Vector{Float64}=[0.0, 1.0], p0::Float64=0.5, method=VCABM(), 
                  eth::Float64=1e-6, cbs_save::CallbackSet=CallbackSet(), dt_s::Float64=0.1, 
                  abstol::Float64=1e-6, reltol::Float64=1e-3, custom=false,
                  dt0::Float64=0.1, dt_min::Float64=1e-7, 
                  fileener::String="ener.txt")
    if N > 0
        if graph.N == 0
            c = K * alpha
            graph = build_ER_HGraph(N, c, K, seed_g)
        end

        if length(links) == 0
            links = gen_links(graph, seed_l)
        end

        if custom
            return CDA_KSAT_custom(ratefunc, rargs_cst, rarg_build, graph, links, tspan, p0, method, eth, 
                                   reltol, dt0, dt_min, fileener)
        else
            return CDA_KSAT_base(ratefunc, rargs_cst, rarg_build, graph, links, tspan, p0,
                                 method, eth, cbs_save, dt_s, abstol, reltol)
        end
    else
        throw("In CDA_KSAT function: The user should provide either a graph::HGraph or valid values for N, K and alpha")
    end  
end

