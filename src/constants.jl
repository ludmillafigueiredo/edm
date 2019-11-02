const Boltz = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const aE = 0.63 #-eV Brown & Sibly MTE book chap 2
const tK = 273.15 # °C to K converter
# These factors adapt the base mortality rates to the developmental stage of the individual. They were defined through a study of the rates yielded by the metabolic theory for species in the Goettingen list
const seed_mfactor = 15
const juv_mfactor = 1
const adult_mfactor = 1

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