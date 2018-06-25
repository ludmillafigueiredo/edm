#!/bin/bash
#SBATCH -n 12 #number of cores
#SBATCH --mem-per-cpu=8G
#SBATCH --output=ED-%j.o
#SBATCH --error=ED-%j.e

/home/luf74xx/builds/julia-d386e40c17/bin/julia main.jl --simID "singlesp1-25yexpgrowth" --landconfig "/home/luf74xx/Dokumente/model/inputs/landpars5m.csv" --spinput "/home/luf74xx/Dokumente/model/inputs/singlespinputj1g.csv" --temp_ts "/home/luf74xx/Dokumente/model/inputs/envtimeseries_199925y.csv" --timesteps "1300"