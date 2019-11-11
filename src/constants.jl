const BOLTZ = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const A_E = 0.63 #-eV Brown & Sibly MTE book chap 2
const T_K = 273.15 # °C to K converter
# These factors adapt the base mortality rates to the developmental stage of the individual. They were defined through a study of the rates yielded by the metabolic theory for species in the Goettingen list
const SEED_MFACTOR = 1

# Normalization constants of metabolic rates
# ------------------------------------------
const B0_GERM = 141363714.221475/100
const B0_MORT = 1113239249.49412
# Growth constant is species-specific

# List of evolvable traits
# ------------------------
const EVOLVABLE_TRAITS = [:compartsize, :span, :firstflower, :floron, :floroff,
      		       	  :seednumber, :seedon, :seedoff, :bankduration]

# Pollination parameters
# ----------------------
# proportion of reproductive mass allocated to seed production
const ALLOC_SEED = 0.05
# default proportion of visited plants
const VISITED_DEFAULT = 1

# Parameters of dispersal kernels
# -------------------------------
const DISPERSAL_PARS = Dict("short" => (mu = 1, lambda = 0.2, factor = 4),
	       	      "medium" => (mu = 0.2, lambda = 3, factor = 1000),
		      "long" => (mu = 1000, lambda = 100, factor = 1))

# Growth parameters
# -----------------
const Q = 5
const C_K = 1000.0

# Probability of management
# -------------------------
const MANAGE_PROB = 1.0

# Various
# -------
# SPP_REF and TRAIT RANGES are also contants, but are input-dependent.
# Therefore, they are initialized in initialisation.jl.

# Minimal value to avoid errors in sum() and Normal()
NOT_0 = 1e-7
