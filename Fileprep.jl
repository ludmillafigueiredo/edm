module Fileprep

# 10 fragments of 5 ha area: 5.10⁴ m² = 5.10⁴.10⁴ = 5.10⁸ cm² (Max. Goettigen fragment area)
# nb of x and y cells (3cm side): √5.10⁸ / 3 = 74500

"""
  gridsizes(realxlen, realylen)
Convert fragments real length values give in `realxlen` and `realylen` from m to the equivalent number of grid cell lengths.
1 m = 100 cm, 3 cm of cell side size
"""

function lengthstocell(realylen::Array{Float64,2}, realylen::Array{Float64,2})
  xlength = round(Int64,()(realxlen.*100)./3), RoundNearestTiesAway)
  ylength = round(Int64,((realylen.*100)./3, RoundNearestTiesAway)
  return xlength, ylength
end

"""
  areatocell(area)
Calculate the side length of a square grid of a fragment of `area` m².
1 m² = 10000 cm², 3 cm of cell side size
"""
function areatocell(area::Array{Float64,2})
  side = round(Int64,((sqrt.(area.*10000))./3), RoundNearestTiesAway)
  return side
end

"""
  connectivity()
Calculate connectivity between ...
"""
function connectivity()
end

end
