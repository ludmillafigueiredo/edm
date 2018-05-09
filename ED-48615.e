ERROR: LoadError: BoundsError: attempt to access ""
  at index [1]
Stacktrace:
 [1] getindex(::String, ::UnitRange{Int64}) at ./strings/string.jl:235
 [2] reproduce!(::Array{Setworld.WorldCell,2}, ::Array{Organisms.Organism,1}, ::Int64, ::Dict{String,Any}) at /home/luf74xx/Dokumente/model/Organisms.jl:347
 [3] simulate() at /home/luf74xx/Dokumente/model/main.jl:186
 [4] include_from_node1(::String) at ./loading.jl:576
 [5] include(::String) at ./sysimg.jl:14
 [6] process_options(::Base.JLOptions) at ./client.jl:305
 [7] _start() at ./client.jl:371
while loading /home/luf74xx/Dokumente/model/main.jl, in expression starting on line 216
