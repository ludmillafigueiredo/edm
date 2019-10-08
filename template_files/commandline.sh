#!/bin/bash
# "simID" is the name the simulation should get

/home/ubuntu/build/julia-9d11f62bcb/bin/julia src/run_EDM.jl --simID "simID" --rseed "100" --outputat "/home/ubuntu/model/EDoutputs" --spinput "/home/ubuntu/model/inputs/weisssingles/weissplants_sppinput.csv" --insect "/home/ubuntu/model/inputs/insects_indep.csv" --initialland "/home/ubuntu/model/inputs/control_4m2.grd" --disturbtype "none" --landmode "artif" --timesteps "5200" --tout "4" --temp_ts "/home/ubuntu/model/inputs/temp1917_2017.csv" --timemsg "true"
