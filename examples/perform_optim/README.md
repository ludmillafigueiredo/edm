# Files to test the performance of the code

There are three pontentially critical "conditions" of a simulation that I understand can slow it down:

- Area (`--initial_land` argument): the bigger the area, the more individuals can occupy it, and thus, the simulation slows down as a consequence of the number of individuals to simulate.
- Number of species  (`--sppinput` argument): the higher the number of species (accompanied by bigger area, that can theoretically sustain all of them), the more likely it is that some species produce more individuals than others, and the simulation is slowed down in the same manner to what is described for the area. In fact, what is important here it to evaluate the species coexistance, i.e. are there species that dominate the community to the point where they exclude many others by producing much more seeds, or much bigger adults, for example?
- Period of time to be simulated (`--temp_ts` argument): the objective of the model is to analyze extinction processes, which usually take up to decades or a couple of centuries to complete. Therefore, after equilibrium is reached, a simulation should simulated at least 200 years, to provide a realistic estimation of species responses.

I have also decreased the frequency of outputs (`--output_freq`): from every 13 time steps, to every 26.

With that in mind, here are some scenarios for which performance should be tested:

100 spp, in 100 m2, over 160 years
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test_100_100_100" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_100m2.grd" 
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 1024 m2 (the closest a square grid can get to 1000), over 100 years: increasing area, which is the most critical

```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test_100_1000_100" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_1024m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 1ha (100 x 100 m2), over 100 years
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_10000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

194 spp, in 1ha, over 100 years
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/194spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_10000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 1ha, over 300 years: realistic period of simulation, but smaller area that what would be expected in .
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_10000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 4ha (200 x 200 m2), over 300 years: this is would be a minimally realistic simulation.

```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_40000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

194 spp, in 4ha, over 300 years: same as before, but increased number of species to evaluate coexistence

```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/194spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_40000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 9ha (300 x 300 m2), over 100 years: increased only size, to see the difference in performance
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_90000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

194 spp, in 9ha, over 100 years: same as before, but increased number of species to evaluate coexistence
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/194spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_160years.csv" 
--initial_land "examples/test/control_90000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

100 spp, in 9ha, over 300 years: it would be great if we could get this to run on a feasible time
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/100spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_90000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

194 spp, in 9ha, over 300 years: same as before, but increased number of species to evaluate coexistence
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/194spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_90000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

194 spp, in 16ha (400 x 400 m2), over 300 years: it would be **awesome** if we could get this to run on a feasible time
```
#!/bin/bash
/your_path_to/julia-X.X.X/bin/julia src/run_EDM.jl 
--simID "test" 
--rseed "777" 
--sppinput "examples/perform_optim/194spp_sppinput.csv" 
--temp_ts "examples/perform_optim/temp_300years.csv" 
--initial_land "examples/test/control_160000m2.grd"
--disturb_type "none" 
--pollination "examples/test/pollination_example.jl" 
--output_freq "26" 
--outputat "."
```

Code used to generate land files and species lists:
```
# landscape files
creatEDnlm(100, 0,"examples/perform_optim")
creatEDnlm(10000, 0,"examples/perform_optim")
creatEDnlm(40000, 0,"examples/perform_optim")
creatEDnlm(90000, 0,"examples/perform_optim")
creatEDnlm(160000, 0,"examples/perform_optim")

# 100 spp
create_sppfile("100spp", rseed = 777 , mode = "random", richp = 100, traitsmode = "unif", traitsdir = "", spplist = NULL, inputsdir = "examples/perform_optim")

# 194
create_sppfile("194spp", rseed = 777 , mode = "random", richp = 194, traitsmode = "unif", traitsdir = "", spplist = NULL, inputsdir = getwd())

# temperature file of 320 years (160 years the same used in simulations)
rbind(readr::read_csv("temp_160years.csv"), 
      readr::read_csv("temp_160years.csv")) %>% 
	  mutate(week = row_number()) %>%
	  write_csv("temp_320years.csv")
```
