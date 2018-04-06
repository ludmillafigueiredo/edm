#!/bin/bash
#SBATCH -n 12 #number of cores
#SBATCH --mem-per-cpu=8G
#SBATCH --output=ED-%j.o
#SBATCH --error=ED-%j.e

/home/luf74xx/builds/julia-d386e40c17/bin/julia main.jl
