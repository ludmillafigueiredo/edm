module Fileprep

"""
This module contains functions related to converting between input, model and output units and formats

# Example:
A fragment of 5 ha area: 5.10⁴ m² = 5.10⁴.10⁴ = 5.10⁸ cm² (Mühlenberg area): nb of x and y cells (3cm side): √5.10⁸ / 3 = 74500
grid
areatocell(50000)

  gridsizes(realxlen, realylen)
Convert fragments real length values give in `realxlen` and `realylen` from m to the equivalent number of grid cell lengths.
1 m = 100 cm, 3 cm of cell side size
"""

function gridsizes(realxlen::Array{Float64,1}, realylen::Array{Float64,1})
  xlength = round(Int64,((realxlen.*100)./5), RoundNearestTiesAway)
  ylength = round(Int64,((realylen.*100)./5), RoundNearestTiesAway)
  return xlength, ylength
end

"""
  lengthtocell(d)
Converts distance values from m to cell size.
"""
function lengthtocell(d::Float64)
  celldist = round(Int64,((d.*100)./5), RoundNearestTiesAway)
  return celldist
end

"""
  areatocell(area)
Calculate the side length of a square grid of a fragment of `area` m².
1 m² = 10000 cm², 3 cm of cell side size
# 2 methods:
- single fragment
- multiple fragments
"""
function areatocell(area::Float64)
  side = round(Int64,((sqrt(area*10000))/5), RoundNearestTiesAway)
  return side
end

function areatocell(area::Array{T,1} where {T<:Number})
  side = round.(Int64,(sqrt.(area*10000))/5, RoundNearestTiesAway)
  return side
end


"""
  connectivity()
Calculate connectivity between ...
"""
function connectivity()
end

end
