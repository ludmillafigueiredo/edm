const BOLTZ = 8.62e-5 #- eV/K Brown & Sibly MTE book chap 2
const A_E = 0.63 #-eV Brown & Sibly MTE book chap 2
const T_K = 273.15 # °C to K converter
# These factors adapt the base mortality rates to the developmental stage of the individual.
# They were defined through a study of the rates yielded by the metabolic theory for
# species in the listed for the calcareous grasslands around Goettingen.
const SEED_MFACTOR = 1

# Normalization constants of metabolic rates
# ------------------------------------------
const B0_GERM = 141363714.221475
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
const VST_DFLT = 1
# Efficiency of insect visit in pollination
const INSCT_EFFC = 0.6
# default proportion of plants pollinated by wind
const WIND_DFLT = 1
# Efficiency of wind pollination
const WIND_EFFC = 0.6
const SELFING_PROBA = 0.5
# Probability of clonal distribution
const CLONAL_PROBA

# Parameters of dispersal kernels
# -------------------------------
const DISPERSAL_PARS = Dict("short" => (mu = 1, lambda = 0.2, factor = 4),
	       	      "medium" => (mu = 0.2, lambda = 3, factor = 1000),
		      "long" => (mu = 1000, lambda = 100, factor = 1))

# Growth parameters
# -----------------
const Q = 5
const C_K = 500.0

# Probability of management
# -------------------------
const MANAGE_PROB = 1.0

# Various
# -------
# SPP_REF and TRAIT RANGES are also contants, but are input-dependent.
# Therefore, they are initialized in initialisation.jl.

# Minimal value to avoid errors in sum() and Normal()
NOT_0 = 1e-7
