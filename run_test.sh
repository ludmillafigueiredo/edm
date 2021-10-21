#!/bin/bash
/your_path_to/julia-1.0.5/bin/julia src/run_EDM.jl --simID "test" --rseed "666" --sppinput "examples/test/sppinput_example.csv" --temp_ts "examples/test/temperature_example.csv" --initial_land "examples/test/4m2.grd" --disturb_type "none" --pollination "examples/test/pollination_example.jl" --output_freq "13" --outputat "."
