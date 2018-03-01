# 10 fragments of 5 ha area: 5.10⁴ m² = 5.10⁴.10⁴ = 5.10⁸ cm² (Max. Goettigen fragment area)
# nb of x and y cells (3cm side): √5.10⁸ / 3 = 74500
# 1 fragments
@time landscape_init(false,74500,74500,1) #OutOfMemory error in ThinkPad alone
# 10 fragments

@time newOrgs!( mylandscape, fgroups, sps, init_stage, init_abund, biomassμ, biomasssd, genotypes, dispμ, dispsd, radius, IDcounter)
