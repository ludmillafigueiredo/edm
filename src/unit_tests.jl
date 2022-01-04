"""
This module contains functions for unit tests.
"""

function check_overwrite(path_results)

    if isdir(abspath(path_results))
       #error("Verify the simulation name: $path_results exists already.")
    end

end

function check_duplicates(plants::Array{Plant,1})

    if length(unique(getfield.(plants, :id))) != length(getfield.(plants, :id))
       error("Duplicated individuals: $(length(unique(getfield.(plants, :id)))) ids and $(length(getfield.(plants, :id))) plants")
    end

end

#Wrapper to parallelize check_ages function
function check_ages_matrix!(plants_matrix::Matrix{Vector{Plant}}, chunks::Array{Tuple{Int,Int,Int,Int}})
    Threads.@threads for chunk in chunks
        for x in chunk[1]:chunk[2]
            for y in chunk[3]:chunk[4]
                check_ages(plants_matrix[x,y])
            end
        end
    end
end

function check_ages(plants::Array{Plant,1})

    #if true in (map(x -> x.age > x.span, plants))
    #   error("Individual older than its span")
    #else
    juvs = filter(x -> x.stage == "j", plants)
    if true in map(x -> x.age > x.firstflower, juvs)
       	log_age()
        error("Immature juvenile")
    end

end
