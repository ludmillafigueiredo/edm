#!/bin/bash
#SBATCH -n 1 #number of cores
#SBATCH --mem-per-cpu=8G
#SBATCH --workdir="/home/luf74xx/Dokumente/model/EDoutputs"
#SBATCH --output=ED-%j.o
#SBATCH --error=ED-%j.e

/home/luf74xx/builds/julia-d386e40c17/bin/julia main.jl
