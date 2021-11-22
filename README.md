# File structure

The source code (`src`) is structured to fit into the ODD framework (Overview, Design concepts and Details) devised by Grimm et al.(2010) for the description of individual-based models. The full description is available in `docs/EDM_ODD.pdf`.
- `inputs.jl`: function for reading in inputs
- `initialisation.jl`: functions and custom-types for initalisation of model data structures
- `submodels_lifecycle.jl`: functions to simulate the biological processes of the model
- `submodels_disturbances.jl`: functions to simulate habitat loss and/or pollination loss. Climate change is also a possible disturbance, but it is defined by an input of temperature time-serie (details in the documentation).
- `entities.jl`: entities are also in this script
- `landconfig`: functions used to derive the size of the simulation arena
- `constants.jl`: this file should contain all constants used in the model, but some might still be scathered around in the code (except the path to source code, which is defined in `run_EDM.jl`
- `auxiliary.jl`: functions and custom-types responsible for operations in the model that do directly realted to the points covered by the ODD structure (inputs, initialisation, submodels, scheduling, outputs).
- `outputs.jl`: functions responsible for organizing and write outputs
- `createnlm.R`: creates neutral landscapes as one of the possible inputs of the model.
- `analysED.R`: functions containing the main analysis of the outputs of the model


`examples`: examples of how input files should look like

To run the model, use the command line:

```
your_path_to/julia-1.6.3/bin/julia run_EDM.jl --arguments
```

`run_test.sh` shows how to run the `test` example.

`docs/TRACE/EDM_TRACE.html`: Detailed model description, including reasoning and decisions behind development, following the TRACE framework (TRAnsparent and Comprehensive Ecological modelling documentation), devised by Grimm et al.(2014).

*References*
Grimm, V., Berger, U., DeAngelis, D. L., Polhill, J. G., Giske, J., & Railsback, S. F. (2010). The ODD protocol: A review and first update. Ecological Modelling, 221(23), 2760–2768. https://doi.org/10.1016/j.ecolmodel.2010.08.019

Grimm, V., Augusiak, J., Focks, A., Frank, B. M., Gabsi, F., Johnston, A. S. A., … Railsback, S. F. (2014). Towards better modelling and decision support: Documenting model development, testing, and analysis using TRACE. _Ecological Modelling_, 280, 129–139. https://doi.org/10.1016/j.ecolmodel.2014.01.018

