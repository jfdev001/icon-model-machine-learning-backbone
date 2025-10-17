# On building ICON with FTorch

The ICON build system uses GNU autotools. You should really only have to 
modify the `configure.ac`.

FTorch is built with cmake. If you add FTorch as a submodule, then you can 
end up just building it by itself. 

Once built with FTorch, just configure and build and then attempt to import 
the FTorch library into some module in the upper atmosphere directory 
(or elsewhere...).

`upper_atmosphere/mo_machine_learning_parametrization.f90` then 
`upper_atmosphere/mo_upatmo_phy_iondrag.f90`
