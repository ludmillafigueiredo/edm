mutable struct plant
    species::AbstractString
    radius_FON::Float64
    a_biomass::Float64
    b_biomass::Float64
end
sp_table[:max_temp]

function read_plantslist()
    countlines("splist_test.txt")
    splist = open("splist_test.txt", "r")
    for name in 1:countlines("splist_test.txt")
        readline(splist)
        # dinamically set strings as variables names: use metaprogramming or Dict instead!
    end
end

function read_ptfs()
    countlines("PTFlist_test.txt")
    ptfs = open("PTFlist_test.txt", "r")
    for name in 1:countlines("PTFlist_test.txt")
        readline(ptfs)
        # dinamically

dummy = plant("dumm",9.2,34.5,23)
print(dummy)

for i in 1:20
    dummy.radius_FON = +(1,dummy.radius_FON)
    println("radius is now $(dummy.radius_FON)")
end

sp = "sp_name"
splist = open("splist_test.txt", "r")

readlines(splist)
readline(splist)

sp_table = readtable("splist_test.csv", header = true)
sp_table[:sps]
sp_table[:sps][8]
sp_table[:max_temp]


function read_ptfs()
    countlines("PTFlist_test.txt")
    ptfs = open("PTFlist_test.txt", "r")
    for name in 1:countlines("PTFlist_test.txt")
        readline(ptfs)
        # dinamically set strings as variables names: use metaprogramming or Dict instead!
    end
end

# strings into variable names -------------------------------------
function string_as_name(s::String,v::Any)
    s = Symbol(s)
	@eval (($s) = ($v))
    return v
    #ex: string_as_name("newvar", 2)
end

# more: https://stackoverflow.com/questions/38986764/save-variable-name-as-string-in-julia
mutable struct plant
    species::AbstractString
    radius_FON::Float64
    a_biomass::Float64
    b_biomass::Float64
end
mutable struct landcell
    temperature::Int64
    plant::plant
end

plant("sp",3,9,8)

land = Matrix{landcell}(3,3)

for i in collect(1:3),  j in collect(1:3)
    if i == j
    land[i,j] = landcell(39,
    plant("sp",2,3,5))
    end
end
land[1,1].plant.species
#------------------------------------------------------------------
# Initialize landscape
landscape = Array{WordCell}(4,5)
landscape[1,2]
typeof(landscape)
methods(Patch)
# So far, all elements of landscape are undefined, so they won't have any fields. I have to set the type - only know how to set individually
landscape[1,1]=Patch(false, true, 0.5, 0.5)
landscape[1,1].pollinators

landinit_teste =  landscape_init(10,10)

mutable struct type2
    something::Float64
end


mutable struct type1
    something::Float64
    locations::Tuple
end

ex1 = type1(23,(8,7))
ex2 = type2(44)
ex3 = type1(3,(7,7))

plantlist = Set([ex1,ex2,ex3])

ex4 = type1(88,(6,8))
push!(plantlist,ex4)

mutable struct simpleplant
    sp::String
    value::Float64
end

list = Set()

for i in 1:countlines("splist_test.txt")
    x = simpleplant(
    readline("splist_test.txt"),
    rand()
    )
    push!(list,x)
end
list

realplantlist = Set()

for i in 1:countlines("splist_test.txt")
    global::realplantlist = Set()
    x = Plant(
    readline("splist_test.txt"),
    rand(), #biomass
    rand(),
    "some_genotype",
    rand(), #seed trait
    (rand(),rand()), #location: make it from tha landscape size
    rand() #FON
    )
    push!(realplantlist,x)
end
realplantlist # getting fields is still a problem

for l in countlines("splist_test.csv")
    if ((sp_list[:min_temp][l] < 34 < sp_list[:max_temp][l])&&   (sp_list[:min_precip][l] < 10 < sp_list[:max_precip][l]))
    possible_sps = sp_list[l,:sps]
else
    continue
end
end

l=5
(sp_list[:min_temp][l] < 34 < sp_list[:max_temp][l])&&
(sp_list[:min_precip][l] < 10 < sp_list[:max_precip][l])

findin((sp_list))

sp_list[.&(sp_list[:min_temp].<34,sp_list[:max_temp].>34,
sp_list[:min_precip].<10,sp_list[:max_precip].>10),:sps]

possible =sp_list[.&(sp_list[:min_temp].<34,sp_list[:max_temp].>34,
sp_list[:min_precip].<10,sp_list[:max_precip].>10),:sps]

rand(possible)

land = Matrix{Int64}(3,3)

for i in collect(1:3),  j in collect(1:3)
    if i == j
        land[i,j] = i+j
    end
end
land

t = (8,9)
typeof(t)


# Radial neighborhood
R = 3 #Rmax of ZOI
ZOI_matrix = zeros(Int64,2*(R+1)+1,2*(R+1)+1)
center_x = R+2 #lands in the middle of the region, since the extends for 2R
center_y = R+2

# Trying for a fixed value of r, which is the distance from the
# NOT GREAT, BUT SOMETHING:
r = 1
#for a in collect(0:90:360) # includes 360, and sets 90 where 0 was
a = 0
while a < 360
    ZOI_matrix[Int(floor(center_x + r*cosd(a))),Int(floor(center_y + r*sind(a)))] = a
    a += 45
#trunc provavelmente nao eh a melhor solucao
end
atan2(4,3) #atan2(y,x) gives the angle to x axis in radians docs.juliaorg

radians_matrix = Array{Float64}(R+2,R+2)
coord_matrix = Array{Tuple,2}(R+2,R+2) # Tuple matrix because it will hold the coordinates in a cartesian plan
distance_matrix = Array{Float64}(R+2,R+2)

for y in collect(1:1:size(radians_matrix,2)), x in collect(1:1:size(radians_matrix,1))
    @inbounds begin
        coord_matrix[x,y] = (y - 1, size(radians_matrix,1) - x) #coordinates to origin - x and y in a cartesian plan are j and i in a matrix, respectly
        # angles in radian to positive x axis
        radian_angle = atan2(coord_matrix[x,y][2], coord_matrix[x,y][1])
        #ZOI_matrix[x,y] = atan2(coord_matrix[x,y]/π)
        radians_matrix[x,y] = Float32(radian_angle) #em algum momento os eixos se inverteram
        #distance_matrix[x,y] = sqrt(
        #(coord_matrix[x,y][1] * sin(ZOI_matrix[x,y]))^2 + #coord_matrix[x,y][1] dah o x com relacao a origem
        #(coord_matrix[x,y][2] * cos(ZOI_matrix[x,y]))^2
        #) # esta faltando alguma correcao para a geometria:
        distance_matrix[x,y] = sqrt(
        (coord_matrix[x,y][1] - 0)^2 + (coord_matrix[x,y][2] - 0)^2
        )
        # coord_matrix[x,y][1] is the first element of the "coordinate pair" (i.e. x) and [2], the second (i.e. y)
        #Float32(r*sqrt(
        #  (sin(atan2(abs(y-center_y), abs(x-center_x))))^2 +
        #  (cos(atan2(abs(y-center_y), abs(x-center_x))))^2
        #  ))
end
end

# distance_matrix for stem > grid:
radians_matrix = Array{Float64}(R+2,R+2)
coord_matrix = Array{Tuple,2}(R+2,R+2)
distance_matrix = Array{Float64}(R+2,R+2)

for y in collect(1:1:size(radians_matrix,2)), x in collect(1:1:size(radians_matrix,1))
    @inbounds begin
        coord_matrix[x,y] = (y - 1, size(radians_matrix,1) - x)
end
end

origin = (coord_matrix[,][],origin = (coord_matrix[,][])
for y in collect(1:1:size(radians_matrix,2)), x in collect(1:1:size(radians_matrix,1))
    @inbounds begin
        distance_matrix[x,y] = sqrt(
        (coord_matrix[x,y][1] - 0)^2 + (coord_matrix[x,y][2] - 0)^2
        )
end
end

#coord_matrix[size(coord_matrix,1),19 is the origin!
atan2(1,0)
atan2(3,1)
Bool(0)

convert(Bool,distance_matrix[:,:] .> R) && distance_matrix[:,:] .= 0

convert(::Type{BoolArrax}, x::BitArray) = x==(true) ? true : x==false ? false : throw(InexactError())

exp.(-distance_matrix)

# putting individual FON is FON matrix
A = Matrix{Real}(8,9)
B = ones(4,4)

# Projecting FONs in all quadrants:
res_FONs_matrix = zeros(12,12)

px = 7
py = 6

function projection_FON(individualResFON,res_FONs_matrix, px,py)
    # Projecting to 1st quadrant:
    for j in 1:(size(individualResFON)[2]), i in 1:(size(individualResFON)[1])
        res_FONs_matrix[(px-size(individualResFON)[1]+i-1), (py-j+1)] =
        individualResFON[i,j]
    end
    # Projecting to 2nd quadrant
    for j in 1:(size(individualResFON)[2]), i in 1:(size(individualResFON)[1])
        res_FONs_matrix[(px-size(individualResFON)[1]+i-1),(py+j-1)] =
        individualResFON[i,j]
    end
    # Projecting to 3rd quadrant:
    for j in 1:(size(individualResFON)[2]), i in 1:(size(individualResFON)[1])
        res_FONs_matrix[(px+size(individualResFON)[1]-i -1), (py + j -1)] =
        individualResFON[i,j]
    end
    # Projecting to 4th quadrant: Do what was done from 2nd to 1st, but from the 3rd
    for j in 1:(size(individualResFON)[2]), i in 1:(size(individualResFON)[1])
        res_FONs_matrix[(px+size(individualResFON)[1]-i -1), (py-j+1)] =
        individualResFON[i,j]
    end
    return res_FONs_matrix
end

# Storing values of the resouce FON in a Dictionnary matrix, associate with the individual that produced it
A = Matrix{Dict}(8,8)
A[1,1] = Dict("p1" => 0.98)
A[2,1] = Dict("p1" => 0.80)
# add value
A[2,1]["key1"] = 0.66
push!(A[1,1],"p1"=>5)
# change
A[2,1]["key1"]  = 4
# When Matrix{Dict}(size) is created, it actually has undefined elements and therefore, those must be defined before adding anything
#no good:
#A = Array{Dict}(size(res_FONs_matrix))
#id ="plantid"
#for j in (py-R):(py+R), i in (px-R-1):(px+R-1)
#    if res_FONs_matrix[i,j] > 0
#        if isassigned(A,i) #problematic because only uses one dimension
#            push!(A[i,j], [id] => res_FONs_matrix[i,j])
#        else
#            A[i,j] = Dict()
#            A[i,j][id] = res_FONs_matrix[i,j]
#        end
#    end
#end
#isassigned(A,(7,7))

res_FONs_idDict = [Dict() for i = 1:size(res_FONs_matrix)[1], j = 1:size(res_FONs_matrix)[2]]
for j in (py-R):(py+R), i in (px-R-1):(px+R-1)
    if res_FONs_matrix[i,j] > 0
        push!(res_FONs_idDict[i,j], [id] => res_FONs_matrix[i,j])
    else
        continue
    end
end
A[7,7]
keys(A[7,7])
values(A[7,7])
A[6:7,7:9]
push!(A[7,7], "plant2id" => 0.4567)
list_sps =Array{String}( keys(A[7,7]))
list_sps

#Ressource competition & Plant growth
#Dummy id → FON
idDict_matrix = [Dict() for i=1:8, j=1:8]
#fill it:
for j in 1:8, i in 1:8
    for x in 1:20
        if (rand() > 0.5)
            push!(idDict_matrix[i,j], string("id",x) => rand())
        end
    end
end
location = (6,6)
R_ZOI = 2
# Scan: scanear a area, mas ficar dentro da FON CIRCULAR
Fa = 0
occupied_cells = 0
for j in (location[2]-R_ZOI):(location[2]+R_ZOI), i in (location[1]-R_ZOI-1):(location[1]+R_ZOI-1)
    if in(string("id",3), collect(keys(idDict_matrix[i,j])))
        occupied_cells += 1
        Fa += sum(values(idDict_matrix[i,j])) - idDict_matrix[i,j][string("id",3)]
        println(i,j)
    end
end
Fa = Fa/occupied_cells
#sqrt((i-px)+(j-py) < R_ZOI could also check whether the distance between the point in the square area that is being "scanned" is inside the circular area of the FON or not.

# Calculating the biomass
using DifferentialEquations
#1 Defining a problem:
f(t,u) = 1.01*u
u0=1/2
tspan = (0.0,1.0)
prob = ODEProblem(f,u0,tspan)
#2 Solving a problem
sol = solve(prob,Tsit5(),reltol=1e-8,abstol=1e-8)
using Plots
using GR
plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
     xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")

#New model tests
# Create 3D array with fragments
function fragment_init(xlength, ylength) #TODO read x and y length external file
	# Create basic fragment "layout"
	fragment = Array{WordCell}(xlength,ylength)
	fill!(fragment,WordCell())

	for y in 1:ylength, x in 1:xlength
		@inbounds begin
			fragment[x,y].temperature = 25
			fragment[x,y].precipitation = 500
		end
	end

    return fragment
end

fragment = fragment_init(8,8)

function landscape_init(n_frags, xlength::Int64, ylength::Int64)
	landscape = Array{WordCell}(xlength,ylength,n_frags)
	#Create landscape 3D array to store fragments
	for f in 1:n_frags
		landscape[:,:,f] = fragment_init(xlength,xlength)
	end
    return landscape
end

landscape_test1 = landscape_init(3,5,5)


for y in 1:5,  x in 1:5
    landscape[x,y].suitability = true
    celltemp = rand(Normal(5,1),1)[1]
    landscape[x,y].temperature =  celltemp #TODO whole fragment is getting the same values
    cellprec = rand(Normal(100,10),1)[1]
    landscape[x,y].precipitation = cellprec
end


randomarray = zeros(Float64, 2,3)

for y in 1:3, x in 1:2
    randomarray[x,y] = rand(Normal(4,1),1)[1]
end

randomarray


# Cartesian indexing
@generated function mysum(A::Array{T,N}) where {T,N}
    quote #creates an expression object
        s = zero(T)
        @nloops $N i A begin
            s += @nref $N A i
        end
        s
    end
end

array = ones(10,10)

for cell in eachindex(array)
    array[cell] += cell
end

array[1:2,2:3]

array[[1,7],[2,5]]

locals = hcat([1,7],[2,5])

for row in 1:size(locals,1)
   array[locals[row,1],locals[row,2]] = 10
end
