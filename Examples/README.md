# Examples

This section of the repository contains three examples to familiarise users with Kariba, the CompPS spectra generated to compute 
the radiative transfer corrections in Kariba, and a set of Python scripts to replicate figures 3,4,6 and 7 of Lucchini et al. 2021. 

The three example codes are as follows:

----------------------------------------------------------------------------------------------------------------------------------
## kariba_particles_examples.cpp

This code highlights how to set up each particle distribution needs to be set up, and also highlights how the three most complex
particle distributions in Kariba (Bknpower, Kappa, Mixed) compare to each other for a near-identical set of parameters. It can 
also be used to replicate figure 1 of the model paper. The classes used are Thermal, Mixed, Kappa and Bknpowerlaw.

----------------------------------------------------------------------------------------------------------------------------------
## kariba_corona_examples.cpp

This code highlights how to set up the an accretion disk+spherical corona model in Kariba, and compares the thermal Comptonisation 
for a range of temperatures and optical depths, and compares the output in Kariba and CompPS. This is also shown in figure 2 of 
the model paper. The classes used are Thermal, ShSDisk and Compton.

----------------------------------------------------------------------------------------------------------------------------------
## kariba_singlezone_examples.cpp

This code highlights how to set up the simplest iteration possible of a single zone blazar emission model, neglecting cooling as
well as external photon fields. Only synchrotron and single-scattering synchrotron-self Compton emission from a spherical region 
are calculated. It replicates the resultes of the EHT MWL campaign of M87, and the default parameters are the same as Model 2 in 
that paper. The classes used are Powerlaw, Cyclosyn and Compton.

----------------------------------------------------------------------------------------------------------------------------------
## Use

All the example codes can be compiled with the ./Makecode command, which will generate all three executables. Unlike BHJet, no
additional wrappers are provided, so integrating these examples in ISIS or Xspec is left to  users.
