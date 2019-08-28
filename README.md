# Development of ecophysiology module in GEOS-Chem 

## Summary

Ground-level ozone is a major air pollutant that adversely affects human health and agricultural productivity. Removal of ozone from the atmosphere by vegetation is controlled mostly by the process of dry deposition in the form of plant stomatal uptake, which in turn causes damage to plant tissues with ramifications for ecosystem and crop health. The openness of plant stomata is generally represented by a bulk stomatal conductance, which is often semi-empirically parameterized in many atmospheric and land surface models, highly fitted to historical observations. A lack of mechanistic linkage to ecophysiological processes such as photosynthesis may render models insufficient to represent plant-mediated responses of atmospheric chemistry to long-term changes in CO2, climate and short-lived air pollutant concentrations. We therefore developed a new ecophysiology module to mechanistically simulate land-atmosphere exchange of important gas species in GEOS-Chem, based on version 12.2.0. We adopted the formulations from the Joint UK Land Environmental Simulator (JULES) to couple photosynthesis rate and bulk stomatal conductance dynamically. The implementation not only allows dry deposition to be coupled with plant ecophysiology, but also enables plant/crop productivity and functions to respond dynamically to atmospheric chemical changes. 

## Code updates

Most of the new stuffs are inside GeosCore/ecophy_mod.F90. I also incoporated several changes in drydep_mod.F, state_diag_mod.F90 and diagnostics_mod.F90 to deal with diagnostics. Updates in other files are mostly about enabling my module to work in GEOS-Chem. When the module is turned on, the dry deposition over vegetated land will be modified. The resulting ozone concentration over land becomes higher.

### Model structure
Currently, my module is called inside the drydep_mod.F when the bulk canopy stomatal resistance is calculated. If ecophysiology is turned on, it replaces the default parameterizations in Wesely scheme.

### Clarifications on my changes in drydep_mod.F
I found that the differene provided by git is not very informative for this file, so I would like to clarify that nost of those changes are involved to interchange the order of calculation of aerodynamic resistance RA and surface conductance RSURFC (which involves removing a DO loop), and include an IF-statement to call the ecophysiology module.

## Run directory updates 

I used some new and edited input files for input.geos and HEMCO_Config.rc, but since they are not included in the source code directory, I did not put them here, and probably I will send them directly to the GEOS-Chem Support Team soon.

## Data file update

To enable the parameterizations in the ecophysiology module, I included a map of some soil parameters. It is read by HEMCO.

# References:
1. Clark, D. B., et al., _The Joint UK Land Environment Simulator (JULES), model description – Part 2: Carbon fluxes and vegetation dynamics_, Geosci. Model Dev., 4, 701–722, [https://doi.org/10.5194/gmd-4-701-2011](https://doi.org/10.5194/gmd-4-701-2011), 2011.
1. Best, M. J., et al. _The Joint UK Land Environment Simulator (JULES), model description – Part 1: Energy and water fluxes_, Geosci. Model Dev., 4, 677–699, [https://doi.org/10.5194/gmd-4-677-2011](https://doi.org/10.5194/gmd-4-677-2011), 2011.
1. P.M Cox, C Huntingford, R.J Harding, _A canopy conductance and photosynthesis model for use in a GCM land surface scheme_, Journal of Hydrology, Volumes 212–213, Pages 79-94, ISSN 0022-1694, [https://doi.org/10.1016/S0022-1694(98)00203-0](https://doi.org/10.1016/S0022-1694(98)00203-0), 1998.
1. Raoult, N. M., _Land-surface parameter optimisation using data assimilation techniques: the adJULES system V1.0_, Geosci. Model Dev., 9, 2833–2852, [https://doi.org/10.5194/gmd-9-2833-2016](https://doi.org/10.5194/gmd-9-2833-2016), 2016.

# Others

Also check out:

- GEOS-Chem: http://www.geos-chem.org
- GEOS-Chem Wiki: http://wiki-geos.chem.org
- TGABI: http://www.cuhk.edu.hk/sci/essc/tgabi/
