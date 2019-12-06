# Analysis code for paper "Soil characteristics and land-use drive bacterial community assembly patterns"

This repository contains the code for analysis and figure generatation for the paper:

Barnett SE, Youngblut ND, Buckley DH, *In press* Soil characteristics and land-use drive bacterial community assembly patterns, *FEMS Microbiology Ecology*

## Contents
Here is a brief description of each of the files:

[Diversity_measures.md](./Diversity_measures.md): Examines alpha and beta diversity of the soil microbial communities including running CAP analysis to examine the how different soil properties explain variation in community composition.
  
  * Figure 1
  * Supplementary figure S2
  * Supplementary figure S5
  

[Neutral_model.md](./Neutral_model.md): Fits the Sloan neutral model to the regional dataset to assess how well this model explains OTU occupancy.

  * Figure 2
  
[Phylogenetic_signal.md](./Phylogenetic_signal.md): Assess the phylogenetic signal in niche preference for pH, percent SOM, and C:N ratio of the bacteria across the region.

  * Supplementary figure S3
  
[bNTI_calculation.md](./bNTI_calculation.md): Calculation of the βNTI metric (Deterministic proceses).

[bNTI_soil_properties.md](./bNTI_soil_properties.md): Looking for correlations between βNTI and variation in soil parameters across all soil samples.

  * Supplementary figure S4

[bNTI_across_landuse.md](./bNTI_across_landuse.md): Testing if βNTI differes across land use regimes. Also looking for correlations between βNTI and variation in soil parameters within each land use regime. 

  * Figure 4

[bNTI_across_pH_groups.md](./bNTI_across_pH_groups.md): Testing if βNTI differes between acidic or pH neutral soils.

[Comm_assembly_processes.md](./Comm_assembly_processes.md): Calculating RCbray metric. Quantifying how different community assembly processes differ across land use regimes and pH groups (Deterministic: Homogeneous selection and variable selection, Stochastic: limited dispersal, homogenizing dispersal, and drift alone). 

  * Figure 3
  * Figure 5
  
## License

This is provided under the MIT license. See [LICENSE](./LICENSE) for more information.


