# File structure
`src`: contains the source code of the model and all scripts related to creating inputs or analysing outputs
- `run_EDM.jl`: set up of simulation and scheduling of biological processes
- `submodels.jl`: functions to simulate the processes of the model; entities are also in this script
- `landconfig`: functions used to derive the size of the simulation arena
- `constants_globalpars.jl`: this file should contain all constants and global parameters, but some might still be scathered around in the code
- `auxfunctions.jl`: contain functions involved in conversion of formats or measures (distances in m to number od grid-cells, for example)
- `createnlm.R`: creates neutral landscapes as one of the possible inputs of the model.
- `analysED.R`: functions containing the main analysis of the outputs of the model

`template_files`: examples of how input files should look like

`EDM_ODD.pdf`: model description
---
# Style guide
Blocks of code Julia are separated by: `# -------------------`
Sections of R code code are names with:
```
#### Title ####
# -------------# # this just makes it look cute
```
`pseudo:` or `#TODO` identifies sections of pseudo code. Not supposed to be infromative messages on the current status of the code
