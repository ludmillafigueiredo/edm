# From Kubisch's SISP.jl:
using ProfileView
@profile simulate(false)
ProfileView.view()
