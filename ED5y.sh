#!/bin/bash
#SBATCH -n 12 #number of cores
#SBATCH --mem-per-cpu=8G
#SBATCH --output=ED-%j.o
#SBATCH --error=ED-%j.e

/home/luf74xx/builds/julia-d386e40c17/bin/julia main.jl --simID "singlesp1-5y" --landconfig "/home/luf74xx/Dokumente/model/inputs/landpars1m.csv" --spinput "/home/luf74xx/Dokumente/model/inputs/singlespinput1.csv" --temp_ts "/home/luf74xx/Dokumente/model/inputs/envtimeseries_19995y.csv" --timesteps "260"