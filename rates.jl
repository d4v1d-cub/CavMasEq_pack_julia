include("graphs.jl")

function rate_FMS_KSAT(Ep::Int64, Em::Int64, T::Float64, avE::Float64, K::Int64)
    dE = Em - Ep
    if dE > 0
        return Ep / K / avE * exp(-dE / T)
    else
        return Ep / K / avE
    end
end