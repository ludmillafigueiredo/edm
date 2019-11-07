const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const aE = 0.63 #-eV Brown & Sibly MTE book chap 2
const tK = 273.15 # °C to K converter
# These factors adapt the base mortality rates to the developmental stage of the individual. They were defined through a study of the rates yielded by the metabolic theory for species in the Goettingen list
const SEED_MFACTOR = 15

const B0_GERM = 141363714.221475
const B0_MORT = 1113239249.49412

const EVOLVABLE_TRAITS = [:compartsize, :span, :firstflower, :floron, :floroff,
      		       	  :seednumber, :seedon, :seedoff, :bankduraiton]

# SPP_REFERENCE and TRAIT RANGES are also contants, but are input-dependent.
# Therefore, there are initialized in initialisation.jl.

# minimal values to avoid errors in sum() and Normal()
NOT_0 = 1e-7
# proportion of reproductive mass allocated to seed production
const ALLOC_SEED = 0.05
# default proportion of visited plants
const VISITED_DEFAULT = 1

# Dispersal kernels
const dispersal_pars = Dict("short" => (mu = 1, lambda = 0.2, factor = 4),
	       	      "medium" => (mu = 0.2, lambda = 3, factor = 1000),
		      "long" => (mu = 1000, lambda = 100, factor = 1))

const Q = 5
const cK = 350.0

# Probability of management
const manage_prob = 1.0