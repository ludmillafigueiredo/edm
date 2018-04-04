#!/bin/bash
#SBATCH -n 1 #number of cores
#SBATCH --mem-per-cpu=2G

/home/luf74xx/builds/julia-d386e40c17/bin/julia main.jl
