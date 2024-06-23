# Script that contains the functions to export the results

function print_ener(saved_e::SavedValues, filename::String)
    fout = open(filename, "w")
    for i in eachindex(saved_e.t)
        write(fout, string(saved_e.t[i]) * "\t" * string(saved_e.saveval[i]) * "\n")
    end
    close(fout)
end


function print_ener(t_list::Vector{Float64}, e_list::Vector{Float64}, filename::String)
    fout = open(filename, "w")
    for i in eachindex(t_list)
        write(fout, string(t_list[i]) * "\t" * string(e_list[i]) * "\n")
    end
    close(fout)
end