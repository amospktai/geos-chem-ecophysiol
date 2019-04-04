This branch is where I test how I can work with diagnostics.

Copied from the ecophysiology branch:

# Development of ecophysiology module in GEOS-Chem 

## Summary

I am working on developing an ecophysiology module in GEOS-Chem. This module aims at better representing atmosphere-biosphere fluxes by mechanistically simulating vegetation processes, e.g. photosynthesis and dry deposition. For the moment, I am still using some of the most simple parameterizations and schemes from other land surface models (Mostly from [JULES](https://jules.jchmr.org), [Clark et. al. (2011)](https://doi.org/10.5194/gmd-4-701-2011) summarizes the science of this model), but having the framework of this module enables easier future modifications.

## Code updates

Most of the new stuffs are inside GeosCore/ecophy_mod.F90. Updates in other files are mostly about enabling my module to work in GEOS-Chem. I have already tested that when the ecophysiology module is turned off, the other modules work exactly the same as in version 12.2.0. When the module is turned on, the dry deposition over vegetated land will be modified. The resulting ozone concentration over land becomes higher after a 1-month simulation in summer, due to a decrease in ozone deposition flux. This is not expected, so I will keep checking and updating this module.

### Model structure
Currently, my module is called inside the drydep_mod.F when the bulk canopy stomatal resistance is calculated. If ecophysiology is turned on, it replaces the default parameterizations in Wesely scheme.

### Clarifications on my changes in drydep_mod.F
I found that the differene provided by git is not very informative for this file, so I would like to clarify that all I did is just interchanging the order of calculation of aerodynamic resistance RA and surface conductance RSURFC (which involves removing a DO loop), and include an IF-statement to call the ecophysiology module.

## Run directory updates 

I used some new and edited input files for input.geos and HEMCO_Config.rc, but since they are not included in the source code directory, I did not put them here, and probably I will send them directly to the GEOS-Chem Support Team soon.

## Data file update

To enable the parameterizations in the ecophysiology module, I included a map of some soil parameters. It is read by HEMCO.

# Documentation of the ecophysiology module

This section is under development.

# Future work

ecophy_mod.F90:
1. Deal with diagnostics
1. Move away from module variables? Especially the PFT-specific parameters.
1. Tidy up the codes: formatting, comments and enable the use of Protex.
1. Add error handlings and floating point exception checkings.
1. Include other parameterizations into the module.

HEMCO_Config.rc:
1. Try using the original soil parameters map with finer resolution, insteaed of the regridded one.

# Others

Also check out:

- GEOS-Chem: http://www.geos-chem.org
- GEOS-Chem Wiki: http://wiki-geos.chem.org
- TGABI: http://www.cuhk.edu.hk/sci/essc/tgabi/
